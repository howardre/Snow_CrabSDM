# Modeling the Distribution of Bering Sea Snow Crab

This is part of a NOAA funded project that will create updated species distribution maps for snow crab in the northern Bering Sea. Portions of this work will comprise a chapter in my PhD dissertation.

### Data
Both fishery-independent, fishery-dependent, and oceanographic data were used for this project, along with ROMS output from the Bering10K. Data cleaning and matching processes can be seen [here](code/data_matching.R), though data are not accessible through this respository. Some data are confidential and others must be requested. Specifically, we are using:
- The **Alaska Department of Fish and Game Crab Observer Program** data, which are __confidential__. These data include records of catches of snow crab in the targeted fishery as well as bycatch in other crab fisheries. The fisheries for crab typically occur from December to March.
- The NOAA AFSC **Eastern Bering Sea Bottom Trawl Survey** data, which includes catches of snow crab at stations sampled by this survey. Data are available beginning in 1975 and surveys run during the summer months. May be available upon request.
- The **Bering10K ROMS** output. We are using the latest CMIP6 runs to obtain values for temperature in order to predict next season distributions. These outputs are [publicly available](https://beringnpz.github.io/roms-bering-sea/B10K-dataset-docs/)
- The **NOAA Eastern Bering Sea sediment database**, which provides a comprehensive set of grain sizes in the study region at a 1 km resolution.
- The **ERA5 Reanalysis** sea ice concentration monthly values.

### Methods
Different types of species distribution models (SDM) were compared in order to select the best model. Root mean square error, Spearman's correlation coefficient, and percent deviance explained were used to compare the models. Ultimately, boosted regression trees were selected. This process can be replicated [here](code/model_evaluation.R).
#### Generalized Additive Models (GAMs)
Two types of GAMs were evaluated for use as species distribution models  for snow crab.
1. Delta-type GAMs model presence-absence first using a Bernoulli distribution. Then abundance-only data is log(x+1) transformed and modelled using a Gaussian distribution. Predicted abundance is conditional on the presence-absence from the first model.
2. GAMs using the full set of log transformed data using a Tweedie distribution with a log link. This is a type of Poisson-gamma compound model. 

#### Boosted Regression Trees (BRTs)
Delta-type BRTs were developed in a similar manner to the delta-type GAMs. Predicted abundance for these models was also conditional on the presence-absence from the Bernoulli model. 