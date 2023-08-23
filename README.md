# TP_predict - Predict TPs and create suspect lists

This collection of scripts allows the user to reproduce the TP prediction and data analyses presented in the following publication:

Trostel, L., Coll, C., Fenner, K., and Hafner, J. (2023). Combining predictive and analytical methods to elucidate pharmaceutical biotransformation in activated sludge. Environ. Sci.: Processes Impacts _25_, 1322–1336. 10.1039/D3EM00161J.
https://doi.org/10.1039/D3EM00161J 

The tools can further be used to perform the same predictions and analyses on a different set of compounds.

## Content

* **TP_prediction**: Script to predict TPs and corresponding biodegradation pathways
* **File_conversion**: Conversion of prediction output to input for suspect screening tools
  * Prediction_output_to_mass_list
  * SMILES_to_mass_and_inclusion_list
* **Additional_analyses**
  * Compare_methods
  * Analyse_cutoff_thresholds

Specific user guidance can be found in the README.md files of the content folders.

## How to
To fetch the code from the git repository, open a terminal and run:
```
$ git clone https://github.com/FennerLabs/TP_predict
```
Go to the newly created directory:
```
$ cd TP_predict
```
To set up TP_predict and install the dependencies, run:
```
$ make
```

## Installation and requirements
The scripts requires rdkit for python, which is easiest installed in a conda environment.
All scripts have been developed and tested in Python version 3.6 on Mac, Linux and Windows operating systems.

### Anaconda step by step guide for non-python users:

1. [Download Anaconda](https://docs.anaconda.com/anaconda/install/index.html) and install it, then run `Anaconda Navigator`
2. create new environment under the `Environment` tab, select python version 3.6.13
3. go to environments, click `play button` on newly created environment, open Terminal
4. run following lines individually (need to confirm: type `y` and press `enter`)(might take a while): `conda install -c rdkit rdkit` and `pip install pubchempy`
5. check if pandas is installed and active according to [this Tutorial](https://docs.anaconda.com/anaconda/navigator/tutorials/pandas/)	
6. open `Anaconda Navigator`, go to `Home` tab, check if `Applications on` is set to the new environment
7. click `gear icon` on `Spyder` > install specific version > 5.0.5 and wait for installation to finish
8. click `launch button` below `Spyder`
