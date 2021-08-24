# Import the modules of interest
import pandas as pd
import numpy as np
from time import sleep
import geopandas as gpd
import subprocess
import tqdm
import time
import datetime
from pathlib import Path
import ee
from functools import partial
from contextlib import contextmanager
import multiprocessing

ee.Initialize()

###########################################################
# Load data
taxonList = list(pd.read_csv('~/data/EMF_listToModel.csv')['unite_sh'])

# If script doesn't run at once due to multiprocessing issues; rerun - lines below filter out jobs currently in queue
taskList = [str(i) for i in ee.batch.Task.list()]
subList = [s for s in taskList if any(xs in s for xs in ['probability'])]
subList = [s for s in subList if any(xs in s for xs in ['READY'])]
runningList = list(map(lambda s: s.strip("<Task EXPORT_IMAGE: ").strip("_probability (READY)>"), subList))

taxonList = [x for x in taxonList if x not in runningList]

###########################################################

# General settings
# Input the name of the username that serves as the home folder for asset storage
usernameFolderString = ''

# Input the name of the classification property
classProperty = 'presence'

# Input the normal wait time (in seconds) for "wait and break" cells
normalWaitTime = 5

# Input a longer wait time (in seconds) for "wait and break" cells
longWaitTime = 10

# Input the Cloud Storage Bucket that will hold the bootstrap collections when uploading them to Earth Engine
# !! This bucket should be pre-created before running this script
bucketOfInterest = ''

# Specify the column names where the latitude and longitude information is stored
latString = 'Pixel_Lat'
longString = 'Pixel_Long'

# Input the name of the property that holds the CV fold assignment
cvFoldString = 'CV_Fold'

# Input a list of the covariates being used
covariateList = [
'CGIAR_PET',
'CHELSA_BIO_Annual_Mean_Temperature',
'CHELSA_BIO_Annual_Precipitation',
'CHELSA_BIO_Max_Temperature_of_Warmest_Month',
'CHELSA_BIO_Precipitation_Seasonality',
'ConsensusLandCover_Human_Development_Percentage',
'ConsensusLandCoverClass_Barren',
'ConsensusLandCoverClass_Deciduous_Broadleaf_Trees',
'ConsensusLandCoverClass_Evergreen_Broadleaf_Trees',
'ConsensusLandCoverClass_Evergreen_Deciduous_Needleleaf_Trees',
'ConsensusLandCoverClass_Herbaceous_Vegetation',
'ConsensusLandCoverClass_Mixed_Other_Trees',
'ConsensusLandCoverClass_Shrubs',
'EarthEnvTexture_CoOfVar_EVI',
'EarthEnvTexture_Correlation_EVI',
'EarthEnvTexture_Homogeneity_EVI',
'EarthEnvTopoMed_AspectCosine',
'EarthEnvTopoMed_AspectSine',
'EarthEnvTopoMed_Elevation',
'EarthEnvTopoMed_Slope',
'EarthEnvTopoMed_TopoPositionIndex',
'EsaCci_BurntAreasProbability',
'GHS_Population_Density',
'GlobBiomass_AboveGroundBiomass',
'GlobPermafrost_PermafrostExtent',
'MODIS_NPP',
'PelletierEtAl_SoilAndSedimentaryDepositThicknesses',
'SG_Depth_to_bedrock',
'SG_Sand_Content_005cm',
'SG_SOC_Content_005cm',
'SG_Soil_pH_H2O_005cm',
]

# Load the composite on which to perform the mapping, and subselect the bands of interest
compositeToUse = ee.Image("").select(covariateList)

# Load a geometry to use for the export
exportingGeometry = ee.Geometry.Polygon([[[-180, 88], [180, 88], [180, -88], [-180, -88]]], None, False);

# Set k for k-fold CV
k = 10

# Make a list of the k-fold CV assignments to use
kList = list(range(1,k+1))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# SDM inputs
# original resolution = 30 arc sec;
# *2 for 1min; *60 for 1 deg; *2 for 2 deg
# dist = compositeToUse.projection().nominalScale().getInfo()*2*60*2

# Number of iteratitions (= re-runs with randomly selected PAs)
n_it = 10

# Make a list of the k-fold CV assignments to use
nList = list(range(1,n_it+1))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Random forest inputs
nTrees = 250

def pipeline(taxonOfInterest):
    success = False
    idle = 0

    while not success:
            try:
                # Write the name of a local staging area folder for outputted CSV's
                holdingFolder = '~/data/training_data/'+taxonOfInterest

                # Create directory to hold training data
                # Path(holdingFolder).mkdir(parents=True, exist_ok=True)

                # Input the name of the project folder inside which all of the assets will be stored
                # This folder will be generated automatically below, if it isn't yet present
                projectFolder = '2021_EMF/'+taxonOfInterest

                ####################################################################################################################################################################
                # Start of modeling
                ####################################################################################################################################################################

                obs_points = ee.FeatureCollection('20210204_EMF_sampled')\
                                    .filterMetadata('unite_sh', 'equals', taxonOfInterest)\
                                    .distinct('.geo')\
                                    .map(lambda f: f.set('presence', 1)).select(covariateList + ['presence'])

                dist = compositeToUse.projection().nominalScale().getInfo()*2*2*60

                buffer = obs_points.map(lambda f: f.buffer(dist))

                PA_full = ee.FeatureCollection('20210204_random_points_sampled').select(covariateList + ['presence'])
                PA_2dfar = PA_full.filter(ee.Filter.geometry(buffer).Not())

                fcList = []
                for n in nList:
                    pa_points = PA_2dfar.randomColumn(seed=n)\
                                        .filterMetadata('random', 'less_than', 0.05)\
                                        .distinct('.geo')\
                                        .limit(obs_points.size())\
                                        .map(lambda f: f.set('presence', 0))

                    training_data = obs_points.merge(pa_points).set('iteration',n);

                    fcList.append(training_data)

                '''
                for n in nList:
                        fc = fcList[n-1]
                        fc = fc.map(lambda f: f.set('run', n).set('unite_sh',taxonOfInterest)).toList(50000).getInfo()

                        result = []

                        for item in fc:
                            values = item['properties']
                            row = [str(values[key]) for key in covariateList+['run','unite_sh','presence']]
                            row = ",".join(row)
                            result.append(row)

                        df = pd.DataFrame([item.split(",") for item in result], columns = covariateList+['run','unite_sh','presence'])
                        df.replace('None', np.nan, inplace = True)
                        with open('data/training_data/'+taxonOfInterest+'_training_data.csv', 'a') as f:
                            df.to_csv(f, mode='a', header=f.tell()==0)
                '''

                # Instantiate RF classifier with default settings
                classifier = ee.Classifier.smileRandomForest(
                                        numberOfTrees = nTrees,
            #                             variablesPerSplit=2,
                                        bagFraction = 0.632
                                        ).setOutputMode('PROBABILITY')



                # Make a feature collection from the k-fold assignment list
                kFoldAssignmentFC = ee.FeatureCollection(ee.List(kList).map(lambda n: ee.Feature(ee.Geometry.Point([0,0])).set('Fold',n)))

                # Define a function to take a feature with a classifier of interest
                def computeAccuracy(fcOI):
                    # Create a function to map through the fold assignments and compute the overall accuracy
                    # for all validation folds
                    def computeAccuracyForFold(foldFeature):
                        # Organize the training and validation data
                        foldNumber = ee.Number(ee.Feature(foldFeature).get('Fold'))
                        trainingData = fcOI.filterMetadata(cvFoldString,'not_equals',foldNumber)
                        validationData = fcOI.filterMetadata(cvFoldString,'equals',foldNumber)

                        # Train the classifier and classify the validation dataset
                        trainedClassifier = classifier.train(trainingData,classProperty,covariateList)
                        outputtedPropName = classProperty+'_probability'

                        # Compute accuracy metrics
                        accuracyToSet = validationData.classify(trainedClassifier,outputtedPropName).errorMatrix(classProperty,outputtedPropName).accuracy()
                        kappaAccuracyToSet = validationData.classify(trainedClassifier,outputtedPropName).errorMatrix(classProperty,outputtedPropName).kappa()
                        return foldFeature.set('accuracy',accuracyToSet).set('kappaAccuracy', kappaAccuracyToSet)

                    # Compute the accuracy values of the classifier across all folds
                    accuracyFC = kFoldAssignmentFC.map(computeAccuracyForFold)
                    meanAccuracy = accuracyFC.aggregate_mean('accuracy')
                    sdAccuracy = accuracyFC.aggregate_total_sd('accuracy')
                    meankappaAccuracy = accuracyFC.aggregate_mean('kappaAccuracy')
                    sdkappaAccuracy = accuracyFC.aggregate_total_sd('kappaAccuracy')

                    # Compute the feature to return
                    featureToReturn = fcOI.set('Mean_Accuracy',meanAccuracy,'StDev_Accuracy',sdAccuracy).set('Mean_KappaAccuracy',meankappaAccuracy,'StDev_KappaAccuracy',sdkappaAccuracy)
                    return featureToReturn

                accuracy_fc = ee.FeatureCollection(ee.List(nList).map(lambda n: ee.Feature(ee.Geometry.Point([0,0])).set('iteration',n)))

                # Perform accuracy k-fold test
                # list(map(computeAccuracy, fcList))

                def varImpFunc(fc):

                    fcOI = ee.FeatureCollection(fc)

                    iteration = fcOI.get('iteration').getInfo()

                    # Train the classifier with the collection
                    trainedClassifer = classifier.train(fcOI,classProperty,covariateList)

                    # Variable importance metrics
                    classifierDict = trainedClassifer.explain().get('importance')
                    featureImportances = classifierDict.getInfo()
                    featureImportances = pd.DataFrame(featureImportances.items(),
                                                      columns=['Covariates', 'Feature_Importance']).sort_values(by='Feature_Importance',
                                                                                                                ascending=False)
                    featureImportances['iteration'] = iteration
                    Path('output/feature_importances/EMF').mkdir(parents=True, exist_ok=True)

                    with open('output/feature_importances/EMF/'+taxonOfInterest+'_featureImportances.csv', 'a') as f:
                        featureImportances.to_csv(f, mode='a', header=f.tell()==0)
                    print(taxonOfInterest, 'variable importance #',iteration)

                    # print('Feature Importances: ', '\n', featureImportances)
                    # plt = featureImportances[:10].plot(x='Covariates', y='Feature_Importance', kind='bar', legend=False,
                    #                               title='Feature Importances')
                    # fig = plt.get_figure()
                    # fig.savefig('output/'+classProperty+'_FeatureImportances.png', bbox_inches='tight')

                list(map(varImpFunc, fcList))

                def bootstrapFunc(fc):
                    # Load the collection with the pre-assigned K-Fold assignments
                    fcOI = ee.FeatureCollection(fc)

                    # Train the classifier with the collection
                    trainedClassifer = classifier.train(fc,classProperty,covariateList)

                    # Classify the image
                    classifiedImage = compositeToUse.classify(trainedClassifer,classProperty+'_probability')

                    return classifiedImage

                # Construct final image to export
                meanImage = ee.ImageCollection.fromImages(list(map(bootstrapFunc, fcList))).reduce(
                    reducer = ee.Reducer.mean()
                ).rename('mean_probability')

                stdDevImage = ee.ImageCollection.fromImages(list(map(bootstrapFunc, fcList))).reduce(
                    reducer = ee.Reducer.stdDev()
                ).rename('stdDev_probability')

                finalImage = ee.Image.cat(meanImage, stdDevImage)

                export_task = ee.batch.Export.image.toAsset(
                    image = finalImage.toFloat(),
                    description = taxonOfInterest+'_probability',
                    assetId = 'EMF_probabilities/'+taxonOfInterest,
                    crs = 'EPSG:4326',
                    crsTransform = '[0.008333333333333333,0,-180,0,-0.008333333333333333,90]',
                    region = exportingGeometry,
                    maxPixels = int(1e13),
                    pyramidingPolicy = {'.default': 'mode'}
                );

                export_task.start()
                print(taxonOfInterest, "started!")

                success = True

            except Exception as e:
                        print(e)
                        success = False
                        idle = (1 if idle > 5 else idle + 1)
                        print("idling for %d" % idle)
                        sleep(idle)

number_of_processes = 50

@contextmanager
def poolcontext(*args, **kwargs):
    pool = multiprocessing.Pool(*args, **kwargs)
    yield pool
    pool.terminate()

if __name__ == '__main__':

    with poolcontext(number_of_processes) as pool:
        results = pool.map(pipeline, taxonList)


def pipeline(taxonOfInterest):
    success = False
    idle = 0

    while not success:
            try:
                # Write the name of a local staging area folder for outputted CSV's
                holdingFolder = '~/data/training_data/'+taxonOfInterest

                # Create directory to hold training data
                # Path(holdingFolder).mkdir(parents=True, exist_ok=True)

                # Input the name of the project folder inside which all of the assets will be stored
                # This folder will be generated automatically below, if it isn't yet present
                projectFolder = '2021_EMF/'+taxonOfInterest

                ####################################################################################################################################################################
                # Start of modeling
                ####################################################################################################################################################################

                obs_points = ee.FeatureCollection('20210108_EMF_presences_to_map')\
                                    .filterMetadata('unite_sh', 'equals', taxonOfInterest)\
                                    .distinct('.geo')\
                                    .map(lambda f: f.set('presence', 1))

                dist = compositeToUse.projection().nominalScale().getInfo()*2*2*60

                buffer = obs_points.map(lambda f: f.buffer(dist))

                PA_full = ee.FeatureCollection('20210204_random_points_sampled').select(covariateList + ['presence'])
                PA_2dfar = PA_full.filter(ee.Filter.geometry(buffer).Not())


                fcList = []
                for n in nList:
                    pa_points = PA_2dfar.randomColumn(seed=n)\
                                        .filterMetadata('random', 'less_than', 0.05)\
                                        .distinct('.geo')\
                                        .limit(obs_points.size())\
                                        .map(lambda f: f.set('presence', 0))

                    training_data = obs_points.merge(pa_points)

                    sampled_training_data = compositeToUse.select(covariateList).sampleRegions(
                            collection = training_data,
                            scale = compositeToUse.projection().nominalScale().getInfo(),
                            tileScale = 16,
                            geometries = True
                        ).set('iteration',n);

                    fcList.append(sampled_training_data)
                '''
                for n in nList:
                        fc = fcList[n-1]
                        fc = fc.map(lambda f: f.set('run', n).set('unite_sh',taxonOfInterest)).toList(50000).getInfo()

                        result = []

                        for item in fc:
                            values = item['properties']
                            row = [str(values[key]) for key in covariateList+['run','unite_sh','presence']]
                            row = ",".join(row)
                            result.append(row)

                        df = pd.DataFrame([item.split(",") for item in result], columns = covariateList+['run','unite_sh','presence'])
                        df.replace('None', np.nan, inplace = True)
                        with open('data/training_data/'+taxonOfInterest+'_training_data.csv', 'a') as f:
                            df.to_csv(f, mode='a', header=f.tell()==0)
                '''

                # Instantiate RF classifier with default settings
                classifier = ee.Classifier.smileRandomForest(
                                        numberOfTrees = nTrees,
            #                             variablesPerSplit=2,
                                        bagFraction = 0.632
                                        ).setOutputMode('PROBABILITY')



                # Make a feature collection from the k-fold assignment list
                kFoldAssignmentFC = ee.FeatureCollection(ee.List(kList).map(lambda n: ee.Feature(ee.Geometry.Point([0,0])).set('Fold',n)))

                # Define a function to take a feature with a classifier of interest
                def computeAccuracy(fcOI):
                    # Create a function to map through the fold assignments and compute the overall accuracy
                    # for all validation folds
                    def computeAccuracyForFold(foldFeature):
                        # Organize the training and validation data
                        foldNumber = ee.Number(ee.Feature(foldFeature).get('Fold'))
                        trainingData = fcOI.filterMetadata(cvFoldString,'not_equals',foldNumber)
                        validationData = fcOI.filterMetadata(cvFoldString,'equals',foldNumber)

                        # Train the classifier and classify the validation dataset
                        trainedClassifier = classifier.train(trainingData,classProperty,covariateList)
                        outputtedPropName = classProperty+'_probability'

                        # Compute accuracy metrics
                        accuracyToSet = validationData.classify(trainedClassifier,outputtedPropName).errorMatrix(classProperty,outputtedPropName).accuracy()
                        kappaAccuracyToSet = validationData.classify(trainedClassifier,outputtedPropName).errorMatrix(classProperty,outputtedPropName).kappa()
                        return foldFeature.set('accuracy',accuracyToSet).set('kappaAccuracy', kappaAccuracyToSet)

                    # Compute the accuracy values of the classifier across all folds
                    accuracyFC = kFoldAssignmentFC.map(computeAccuracyForFold)
                    meanAccuracy = accuracyFC.aggregate_mean('accuracy')
                    sdAccuracy = accuracyFC.aggregate_total_sd('accuracy')
                    meankappaAccuracy = accuracyFC.aggregate_mean('kappaAccuracy')
                    sdkappaAccuracy = accuracyFC.aggregate_total_sd('kappaAccuracy')

                    # Compute the feature to return
                    featureToReturn = fcOI.set('Mean_Accuracy',meanAccuracy,'StDev_Accuracy',sdAccuracy).set('Mean_KappaAccuracy',meankappaAccuracy,'StDev_KappaAccuracy',sdkappaAccuracy)
                    return featureToReturn

                accuracy_fc = ee.FeatureCollection(ee.List(nList).map(lambda n: ee.Feature(ee.Geometry.Point([0,0])).set('iteration',n)))

                # Perform accuracy k-fold test
                # list(map(computeAccuracy, fcList))

                def varImpFunc(fc):
                    fcOI = ee.FeatureCollection(fc)

                    iteration = fcOI.get('iteration').getInfo()

                    # Train the classifier with the collection
                    trainedClassifer = classifier.train(fc,classProperty,covariateList)

                    # Variable importance metrics
                    classifierDict = trainedClassifer.explain().get('importance')
                    featureImportances = classifierDict.getInfo()
                    featureImportances = pd.DataFrame(featureImportances.items(),
                                                      columns=['Covariates', 'Feature_Importance']).sort_values(by='Feature_Importance',
                                                                                                                ascending=False)
                    Path('output/feature_importances/'+taxonOfInterest).mkdir(parents=True, exist_ok=True)
                    featureImportances.to_csv('output/feature_importances/'+taxonOfInterest+'/'+taxonOfInterest+'_featureImportances_'+str(iteration)+'.csv')
                    # print('Feature Importances: ', '\n', featureImportances)
                    # plt = featureImportances[:10].plot(x='Covariates', y='Feature_Importance', kind='bar', legend=False,
                    #                               title='Feature Importances')
                    # fig = plt.get_figure()
                    # fig.savefig('output/'+classProperty+'_FeatureImportances.png', bbox_inches='tight')

                list(map(varImpFunc, fcList))

                def bootstrapFunc(fc):
                    # Load the collection with the pre-assigned K-Fold assignments
                    fcOI = ee.FeatureCollection(fc)

                    # Train the classifier with the collection
                    trainedClassifer = classifier.train(fc,classProperty,covariateList)

                    # Classify the image
                    classifiedImage = compositeToUse.classify(trainedClassifer,classProperty+'_probability')

                    return classifiedImage

                # Construct final image to export
                meanImage = ee.ImageCollection.fromImages(list(map(bootstrapFunc, fcList))).reduce(
                    reducer = ee.Reducer.mean()
                ).rename('mean_probability')

                export_task = ee.batch.Export.image.toAsset(
                    image = meanImage.toFloat(),
                    description = taxonOfInterest+'_probability',
                    assetId = 'users/johanvandenhoogen/2021_EMF/EMF_probabilities/'+taxonOfInterest,
                    crs = 'EPSG:4326',
                    crsTransform = '[0.008333333333333333,0,-180,0,-0.008333333333333333,90]',
                    region = exportingGeometry,
                    maxPixels = int(1e13),
                    pyramidingPolicy = {'.default': 'mode'}
                );

                # export_task.start()
                print("Mean image", taxonOfInterest, "started!")

            except Exception as e:
                        print(e)
                        success = False
                        idle = (1 if idle > 5 else idle + 1)
                        print("idling for %d" % idle)
                        sleep(idle)

number_of_processes = 12

@contextmanager
def poolcontext(*args, **kwargs):
    pool = multiprocessing.Pool(*args, **kwargs)
    yield pool
    pool.terminate()

if __name__ == '__main__':

    with poolcontext(number_of_processes) as pool:
        results = pool.map(pipeline, taxonList)
