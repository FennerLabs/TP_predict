# Generate Mass and Inclusion Lists from Prediction Output

## Purpose

The purpose of the script is to generate a mass lists, which can be imported into `Compound Discoverer` from multiple `find_best_TP.py` outputs.
Further, the script combines the data of all methods and can select transformation products (TP) based on a scoring system.
It also creates inclusion lists for positive and negative mode with predicted NCE for the use on a QExactivePlus.
The script uses the `PubChem` database and for every compound a request is sent, thus the code can run for some time.

## Input

The scripts needs at least 1 input file, but it supports up to 8 files.
The input file(s) consist of the tsv results from `find_best_TP.py`. The filepaths need to be changed at the start of the script.
The prediction method (enter `EAWAG-PPS` or `envipath`) and package (do not use `.` or `/`) for each file needs to be specified in the script.
Additionally, each optinal file can be considered in the processing (enter "yes") or not.
Moreover, 3 txt files for mapping are required:
* `SMILES_selected_comp.txt`	contains the SMILES of the selected parent compounds
* `name_subst.txt`		contains the full names of the parents
* `code_subst.txt`		contains the short code of the parents (exp. first 3 letters of name)

* Each txt file contains one string per line and the order needs to be matched over the three txt files, so each line of the txt files corresponds to one parent.
These txt files can be easily created by copying from a Excel Worksheet and pasting into a new txt file.

## Scoring System
Can be turned on (`True`) or off (`False`, default). A maximal number of allowed TPs per parent can be set as a variable.
The scoring system removes all TPs, which have a mass below 100 u.
The score is lowered if the TPs has a CAS number, has a low probability or is not predicted by all methods.
If a parent has more than the specified maximal allowed TPs per parent, they get removed starting with lowest score.

## Output

The script generates a multitude of output files. Most of them are pickle files, which contain all the data from a given dictionary.
These can be copied and moved to a different location to run just one part of the code or to save on code running time.
Each of the input files get an individual csv file that can be imported into `Compound Discoverer`.
If more than two input files are processed, 3 csv files are created from the combination of all methods.
2 contain all the data about the parents and TPs before and after the scoring, which is meant to be an overview.
The last one contains only the most important combined data with SMILES, Name, CAS, Formula, monoisotopic mass and InchiKey, which can be imported into `Compound Discoverer` as a mass list.

Additionally, two inclusion lists as well as a txt file is created containing the maximal element count of C, H, O, N, S, P, Cl, Br an I (useful for `Compound Discoverer` workflow).
The output paths set per default to an output folder that is in the same folder as the script.

## Run script in Spyder

1. open the script: `get_mass_list_from_prediction.py`
2. change the pathway of the input files according to the filepath on your harddrive and turn Scoring System on or off.
3. click somewhere in the left window and press `F5` or click `Run File`, if Run settings windows appears, click `Run`
4. wait for the script to finish running (may take a while, `Script finished successfully` marks end of script)


Author: Leo Trostel, 2022