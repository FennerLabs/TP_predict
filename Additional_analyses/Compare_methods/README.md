# Compare performance of TP prediction models

## Purpose
Compares predictions from EAWAG/BBD-PPS and enviPath with found transformation products (TPs) and calculates precision for each method.
The script uses the output of the script "get_mass_list_from_prediction.py"

## Input
The script uses 2 input files:
	o txt file with one SMILES per line of TPs that were found experimentally
	o pickle file of combined dictionary with all predictions made by different methods

## Output
The script outputs four files containing information of number of TPs predicted and found.

## Note
The csv and tsv files can be opened in Microsoft Excel and brought to a more readable format.

## Run script in Spyder

7. open the script: "compare_methods_pps.py"
8. change the pathway of the input files according to the filepath on your harddrive
9. click somewhere in the left window and press F5 or click "Run File", if Run settings windows appears, click run
10. wait for the script to finish running (you can see the progress in the bottom right window, "Script finished successfully" marks end of script)

Author: Leo Trostel, 2022
