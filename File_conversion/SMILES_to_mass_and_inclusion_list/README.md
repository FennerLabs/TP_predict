# Generate Mass and Inclusion Lists

## Purpose

The purpose of the script is to generate a mass lists, which can be imported into `Compound Discoverer`.
It also creates inclusion lists for positive and negative mode with predicted NCE for the use on a QExactivePlus.
The script uses the `PubChem` database and for every compound a request is sent, thus the code can run for some time.

## Input

The scripts needs one txt file as input that contains the SMILES of compounds (one string per line).

## Output

Mass and two inclusion lists as well as an overview with SMILES, Name, CAS, Formula, monoisotopic mass and InchiKey as csv file.
Additionally, a txt file is created containing the maximal element count of C, H, O, N, S, P, Cl, Br an I (useful for `Compound Discoverer` workflow).
The output paths are per default the same folder as the script.

## Run script in Spyder

1. open the script: `smiles_to_mass_and_inclusion_list.py`
2. change the pathway of the input files according to the filepath on your harddrive
3. click somewhere in the left window and press `F5` or click `Run File`, if Run settings windows appears, click `Run`
4. wait for the script to finish running (may take a while, `Script finished successfully` marks end of script)


Author: Leo Trostel, 2022