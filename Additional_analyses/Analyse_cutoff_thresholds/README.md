# Analyze cutoff thresholds

## Purpose
This script evaluates the effect of choosing different cutoff thresholds in the TP predictions. 
The evaluated parameters are `generation` and `max TPs`.

## Input
The script uses 5 input files, one for each tested TP prediction model. 
1. Output from `find_best_TPs.py`
2. Output from EAWAG-PPS batch prediction: `resuts_EAWAG-PPS.tsv`
3. List of experimentally confirmed TPs (SMILES): `confirmed_TPs.txt`

## Output
The script outputs the number of true positive predictions for a range of thresholds, for each parameter searately:
* `Cutoff_thresholds_Generation threshold.txt`
* `Cutoff_thresholds_Maximum TP threshold.txt`

The results are plotted as PDF:
* `Plot_cutoff_analysis.pdf`

## Run script in Spyder

1. open the script: `compare_methods_pps.py`
2. change the pathway of the input files according to the filepath on your harddrive
3. click somewhere in the left window and press F5 or click "Run File", if Run settings windows appears, click run
4. wait for the script to finish running (you can see the progress in the bottom right window, "Script finished successfully" marks end of script)

Author: Jasmin Hafner, 2022
