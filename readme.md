# A global atlas of Earth's mycorrhizal fungi
Which fungi live where? What controls their distribution? Here we globally map the 390 most commonly observed mycorrhizal fungi to identify how different taxa vary in their environmental preferences.

### 1. Data Construction
These scripts transform raw data from GenBank and Global Fungi databases into presence only dataframes, the core input to our mapping pipeline.

### 2. Mapping
These scripts take our formatted and geo-referenced mycorrhizal species presence data, pair them with environmental variables of interest, and then fits machine learning models that are used to map each taxa's distribution globally. We also include code to extract key features of each species' distribution at difference occurrence probability.

### 3. Figure Scripts
These scripts visualize our findings.
