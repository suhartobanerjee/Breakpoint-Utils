#!/bin/bash

# this can be automated for every env
# exporting each of the env serially
# will parallelize later
# copykat is in infercnv env only.
# since both rely on jags and r-jags

source ~/.bashrc



conda activate digKar
printf "\nExporting currently active env: %s" $(echo $CONDA_DEFAULT_ENV)
conda env export > ./${CONDA_DEFAULT_ENV}_env.yml
printf "\n%s env exported successfully to the environments folder!\n" $(echo $CONDA_DEFAULT_ENV)
