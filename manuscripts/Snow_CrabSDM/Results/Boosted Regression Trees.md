**All Models
- Variables: phi, doy, ice, sst, depth, immature individuals, mature individuals
- Gaussian
- Scaled CPUE (due to inability to use large values in gbm package)
- No need to drop variables for either best model
- Want lower deviance, higher correlation

**Females**
_Model 1_ 
- Tree complexity: 5
- Learning rate: 0.05
- Bag fraction: 0.5
- Mean deviance: 0.00043
- SE deviance: 3.81 * 10$^{-5}$
- Mean correlation: 0.49
- SE correlation: 0.025
- Trees: 2050

_Model 2
- Tree complexity: 5
- Learning rate: 0.01
- Bag fraction: 0.5
- Mean deviance: 0.00043
- SE deviance: 2.98 * 10$^{-5}$
- Mean correlation: 0.48
- SE correlation: 0.018
- Trees: 4700

_Model 3 
- Tree complexity: 3
- Learning rate: 0.01
- Bag fraction: 0.5
- Mean deviance: 0.00043
- SE deviance: 3.95 * 10$^{-5}$
- Mean correlation: 0.47
- SE correlation: 0.016
- Trees: 9400

_Model 4
- Tree complexity: 10
- Learning rate: 0.01
- Bag fraction: 0.5
- Mean deviance: 0.00043
- SE deviance: 4.39 * 10$^{-5}$
- Mean correlation: 0.49
- SE correlation: 0.017
- Trees: 2150

_Model 5
- Tree complexity: 5
- Learning rate: 0.1
- Bag fraction: 0.5
- Mean deviance: 0.00044
- SE deviance: 4.18 * 10$^{-5}$
- Mean correlation: 0.48
- SE correlation: 0.022
- Trees: 550

_Model 6
- Tree complexity: 5
- Learning rate: 0.1
- Bag fraction: 0.75
- Mean deviance: 0.00043
- SE deviance: 2.62 * 10$^{-5}$
- Mean correlation: 0.50
- SE correlation: 0.017
- Trees: 900

_Model 7
- Tree complexity: 5
- Learning rate: 0.05
- Bag fraction: 0.75
- Mean deviance: 0.00043
- SE deviance: 4.18 * 10$^{-5}$
- Mean correlation: 0.49
- SE correlation: 0.019
- Trees: 1200

**Males**
_Model 1_
- Tree complexity: 5
- Learning rate: 0.05
- Bag fraction: 0.5
- Mean deviance: 0.0023
- SE deviance: 2.12 * 10$^{-5}$
- Mean correlation: 0.70
- SE correlation: 0.0024

_Model 2
- Tree complexity: 5
- Learning rate: 0.01
- Bag fraction: 0.5
- Mean deviance: 0.0027
- SE deviance: 2.27 * 10$^{-5}$
- Mean correlation: 0.65
- SE correlation: 0.0026

_Model 3 - Chosen due to reasonable number of trees
- Tree complexity: 3
- Learning rate: 0.01
- Bag fraction: 0.5
- Mean deviance: 0.0030
- SE deviance: 2.00 * 10$^{-5}$
- Mean correlation: 0.60
- SE correlation: 0.003

_Model 4
- Tree complexity: 10
- Learning rate: 0.01
- Bag fraction: 0.5
- Mean deviance: 0.0024
- SE deviance: 2.95 * 10$^{-5}$
- Mean correlation: 0.69
- SE correlation: 0.002

_Model 5
- Tree complexity: 10
- Learning rate: 0.01
- Bag fraction: 0.75
- Mean deviance: 0.0024
- SE deviance: 1.85 * 10$^{-5}$
- Mean correlation: 0.69
- SE correlation: 0.0022

_Model 6
- Tree complexity: 5
- Learning rate: 0.1
- Bag fraction: 0.5
- Mean deviance: 0.0023
- SE deviance: 2.26 * 10$^{-5}$
- Mean correlation: 0.71
- SE correlation: 0.0024

_Model 7
- Tree complexity: 5
- Learning rate: 0.1
- Bag fraction: 0.75
- Mean deviance: 0.0023
- SE deviance: 2.17 * 10$^{-5}$
- Mean correlation: 0.71
- SE correlation: 0.0022

**Variable Importance**
_Females_
See figures
Top three are

_Males_
See figures
Top three are phi > sst > immature males
