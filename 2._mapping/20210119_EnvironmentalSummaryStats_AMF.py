import pandas as pd
from time import sleep
from itertools import repeat
from functools import partial
from contextlib import contextmanager
import multiprocessing
import ee
ee.Initialize()

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
compositeToUse = ee.Image("projects/crowtherlab/Composite/CrowtherLab_Composite_30ArcSec").select(covariateList)

predictedProbs = ee.ImageCollection('projects/crowtherlab/johan/2021_fungi_SDMs/AMF_probabilities');

# list of AMF species
AMF_list = list(pd.read_csv('/Users/johanvandenhoogen/ETH/Projects/AMF/data/20210119_VTX_listToModel.csv')['MAARJAM_ID'])

# List of thresholds to use as mask
threshold_list = [0.25, 0.5, 0.75, 0.9, 0.95]

# create list with species + thresholds
mapList = []
for item in threshold_list:
    mapList = mapList + (list(zip(AMF_list, repeat(item))))

# AMF_list, threshold_list = zip(*mapList)

def summary_func(taxon, threshold):
    success = False
    idle = 0

    result = []
    while not success:
        try:

            selectedImage = predictedProbs.filterMetadata('system:index','equals',taxon).toBands().select(taxon+'_mean_probability')

            reducers = ee.Reducer.mean().combine(
              reducer2 = ee.Reducer.median(),
              sharedInputs = True
            ).combine(
              reducer2 = ee.Reducer.min(),
              sharedInputs = True
            ).combine(
              reducer2 = ee.Reducer.max(),
              sharedInputs = True
            ).combine(
              reducer2 = ee.Reducer.stdDev(),
              sharedInputs = True
            );

            unboundedGeo = ee.Geometry.Polygon([[[-180, 88], [180, 88], [180, -88], [-180, -88]]], None, False);

            # Use the combined reducer to get the mean and SD of the image.
            stats = compositeToUse.updateMask(selectedImage.gte(threshold)).reduceRegion(
                reducer = reducers,
                bestEffort = True,
                geometry = unboundedGeo,
                scale = 927.6624232772797
            )

            taxonVarVals = pd.DataFrame(stats.getInfo().items(), columns=['Variable', 'Value'])
            taxonVarVals['Taxon'] = taxon

            fullData = taxonVarVals.pivot(index='Taxon', columns='Variable', values='Value')

            with open('output/threshold_' + str(threshold) + '_AMF_reducedFeatures.csv', 'a') as f:
                fullData.to_csv(f, mode='a', header=f.tell()==0)

            print('processed:',taxon, threshold)
            success = True
        except Exception as e:
                    print(e)
                    success = False
                    idle = (1 if idle > 5 else idle + 1)
                    print("idling for %d" % idle)
                    sleep(idle)


number_of_processes = 40

@contextmanager
def poolcontext(*args, **kwargs):
    pool = multiprocessing.Pool(*args, **kwargs)
    yield pool
    pool.terminate()

if __name__ == '__main__':

    with poolcontext(number_of_processes) as pool:
        results = pool.starmap(summary_func, mapList)
