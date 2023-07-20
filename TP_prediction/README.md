# TP_prediction - Predict Transformation Products (TPs)

TP_prediction predicts TPs and associated biodegradation pathways using the enviPath pathway prediction engine.
To run the script and to save the output of your TP prediction on enviPath, you need a user account on envipath.org.

## Input
Add yor input compounds to ./input/input_structures.tsv. The input format should be the smiles of the compound, followed by its name (or identifier), and separated by  a tab.

## Mandatory settings
* **USERNAME**: Enter here the username of your enviPath account.
* **EP_PACKAGE_ID**: Create a new package for your results on envipath.org and enter its URI here.

## Optional settings
By default, the search will predict the 50 TPs with the highest probability to be observed
according to the relative reasoning model 

You can adapt settings under PATHWAY SEARCH SETTINGS in find_best_TPs.py:

* **EP_MODEL_ID**: URI of enviPath relative reasoning model to be used for prediction
* **MAX_TP**: Maximum number of TPs to predict per input compound
* **PROBABILITY_THRESHOLD**: Lower probability threshold -  any value equal to or lower than the threshold will be excluded
* **INCLUDE_0_PROBABILITIES**: Set probabilities of 0 to 0.01 to continue having a weighting scheme downstream of the pathway
* **MOIETY**: Follow a chemical moiety - only compounds containing this moiety in SMILES will be expanded, e.g., "C(F)(F)F"
* **SORT_TPS_BY_SIZE**: To prioritize small compounds in the node queue
* **FOLLOW_LABELED_ATOM**: Follow labeled atoms - only compounds containing at least one atom labeled with ATOM_LABEL will be expanded 
* **ATOM_LABEL**: Label used to follow atoms, e.g., 14 for radiolabeled carbon

Output settings

* **INPUT_FILE_PATH**: path to input file
* **OUTPUT_DIRECTORY**: path to output directory
* **OUTPUT_FILE_TAG**: Name tag to be added to your output files

## Run prediction
To predict pathways for all compounds specified in the input, run:
```
$ python find_best_TPs.py
```

## Output
The output file containing all predicted pathways and TP information will be stored in ./output.
Each pathway entry starts with '///', followed by the name of the pathway and the link to the pathway entry 
on envipath.org. Each pathway entry is followed by a tab-separated table containing the following information:
* **SMILES**: SMILES of TP
* **name**: Name of TP (automatically generated)
* **combined_probability**: probability of the node (p_node = p_edge * p_node,parent)
* **rules**: List of biotransformation rules used to predict this reaction
* **generation**: number of iteration where TP was generated
* **probability**: probability of the reaction from the parent to this TP (p_edge)
* **parent**: SMILES of parent compound

The output of the TP prediction can be directly used as input for 
* `File_conversion/Prediction_output_to_mass_list/get_mass_list_from_prediction.py`
* `Additional_analyses/Analyse_cutoff_thresholds.py`

Author: Jasmin Hafner, 2022