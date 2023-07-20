# -*- coding: utf-8 -*-
"""
Created on Tue Sep 13 11:51:15 2022

@author: trostele
"""
# Python version 3.6.13
################################################################################################################################################################################################################################################################
# import all necessary packages
import pandas as pd
import rdkit # rdkit is only supported before Python 3.7
from rdkit import Chem
from rdkit.Chem.MolStandardize import rdMolStandardize
from rdkit.Chem import Descriptors
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
import pickle
################################################################################################################################################################################################################################################################

# script to compare predicted TPs with found TPs from different methods

################################################################################################################################################################################################################################################################

# INPUT

found_tp_smiles_input = "./input/found_TP_SMILES.txt"
# contains SMILES of all TPs that were found in samples

pickle_file_data_dict = "./input/data_dict_com_with_CAS.pickle"
# pickle file location of combined dictionary with all predictions made by different methods

package_method_1 = "EAWAG_BBD-PPS_round_2"
package_method_2 = "EAWAG_BBD-PPS_round_2b"
package_method_5 = "enviPath-BBD_1"
package_method_6 = "enviPath-BBD+SOIL_2"
package_method_7 = "enviPath-BBD+SLUDGE_3"
package_method_8 = "enviPath-BBD+SOIL+SLUDGE_4"

# enter used package (source name in combined data dict)
# same input as in "get_mass_list_from_prediction.py"

#!!!
# I used two rounds of predictions with the EAWAG/BBD-PPS, so it combines them
# the script only works with 7 methods, needs to be updated!
#!!!

################################################################################################################################################################################################################################################################

# FUNCTIONS

def canonicalize_smiles(smiles):
    mol = Chem.MolFromSmiles(smiles)  # creates mol object from SMILES
    uncharger = rdMolStandardize.Uncharger() # easier to access
    uncharged = uncharger.uncharge(mol)  # protonates or deprotonates the mol object
    new_smiles = rdkit.Chem.rdmolfiles.MolToSmiles(uncharged)  # converts mol object to canonical SMILES
    can_smiles = Chem.CanonSmiles(new_smiles)
    return can_smiles

def do_pickle(d, pickle_file):
    with open(pickle_file, 'wb') as handle:
        pickle.dump(d, handle, protocol=pickle.HIGHEST_PROTOCOL)

def get_pickle(pickle_file):
    with open(pickle_file, 'rb') as handle:
        d = pickle.load(handle)
    return d

def export(data_dict, smi_list):
    df_name_list_com = []
    df_smi_list_com = []
    df_ID_list_com = []
    df_Formula_list_com = []
    df_MolWeight_list_com = []
    df_name_parent_list_com = []
    df_inchikey_list_com = []
    df_source_list_com = []
    df_alt_parent_list_com = []   
    # add all the data to the lists
    for parent in data_dict:
        df_name_list_com.append(data_dict[parent]["code_parent"][0])
        df_smi_list_com.append(parent)
        df_ID_list_com.append(data_dict[parent]["ID_parent"])
        df_MolWeight_list_com.append(data_dict[parent]["mass_parent"][0])
        df_Formula_list_com.append(data_dict[parent]["Formula_parent"])
        df_name_parent_list_com.append(data_dict[parent]["name"][0])
        df_inchikey_list_com.append(data_dict[parent]["inchi_parent"])
        df_source_list_com.append("")
        df_alt_parent_list_com.append([])
        for tp in data_dict[parent]["TP_dict"]:
            if tp in predicted_and_found:
                df_name_list_com.append(data_dict[parent]["TP_dict"][tp]["code"])
                df_smi_list_com.append(tp)
                df_ID_list_com.append(data_dict[parent]["TP_dict"][tp]["CAS"])
                df_MolWeight_list_com.append(data_dict[parent]["TP_dict"][tp]["mass"])
                df_Formula_list_com.append(data_dict[parent]["TP_dict"][tp]["Formula"])
                df_name_parent_list_com.append(data_dict[parent]["name"][0])
                df_inchikey_list_com.append(data_dict[parent]["TP_dict"][tp]["InchiKey"])
                df_source_list_com.append(data_dict[parent]["TP_dict"][tp]["source_list"])
                df_alt_parent_list_com.append(data_dict[parent]["TP_dict"][tp]["alternative_parent"])        
    df_complete_dict = {"Name of parent":df_name_parent_list_com,"SMILES": df_smi_list_com,
                   "Name": df_name_list_com,
                   "Source": df_source_list_com,
                   "CAS": df_ID_list_com,"Formula": df_Formula_list_com,
                   "MolWeight": df_MolWeight_list_com, "InchiKey":df_inchikey_list_com, "Alternative parent": df_alt_parent_list_com}
    df_complete = pd.DataFrame.from_dict(df_complete_dict)
    df_complete.to_csv("./output/predicted_and_found_TPs.csv", index = False, sep = ",")
 
################################################################################################################################################################################################################################################################

# export overview of predicted TPs that were also found (on txt file) as csv file

found_tp_smiles = []

SMILES_comp_file = open(found_tp_smiles_input)
for line in SMILES_comp_file:
    found_tp_smiles.append(line.rstrip())

found_tp_smiles_canon = []

for tp in found_tp_smiles:
    found_tp_smiles_canon.append(canonicalize_smiles(tp))

data_dict = get_pickle(pickle_file_data_dict)

Bar =  'CCS(=O)(=O)N1CC(CC#N)(n2cc(-c3ncnc4[nH]ccc34)cn2)C1'
Abe =  'CCN1CCN(Cc2ccc(Nc3ncc(F)c(-c4cc(F)c5nc(C)n(C(C)C)c5c4)n3)nc2)CC1'

# remove Abe and Bar from dict because they were not spiked
data_dict.pop(Abe)
data_dict.pop(Bar)

predicted_tp_smiles = []

for parent in data_dict:
    for tp in data_dict[parent]["TP_dict"]:
        predicted_tp_smiles.append(tp)

predicted_and_found = []


for smi in found_tp_smiles_canon:
    if smi in predicted_tp_smiles:
        predicted_and_found.append(smi)

export(data_dict, predicted_and_found)

################################################################################################################################################################################################################################################################

# export precision (= found/predicted) for each method as tsv file

source_list_nest = []

for parent in data_dict:
    for tp in data_dict[parent]["TP_dict"]:
        if tp in predicted_and_found:
            source_list_nest.append(data_dict[parent]["TP_dict"][tp]["source_list"])

source_list = []
for sublist in source_list_nest:
    for item in sublist:
        source_list.append(item)

pps = []
envi_updated = []
envi_soil_updated = []
envi_sludge = []
envi_soil_sludge = []

for source in source_list:
    if "EAWAG_BBD-PPS" in source:
        pps.append(source)
    if package_method_5 in source:
        envi_updated.append(source)    
    if package_method_6 in source:
        envi_soil_updated.append(source)
    if package_method_7 in source:
        envi_sludge.append(source)
    if package_method_8 in source:
        envi_soil_sludge.append(source)

pps_predicted = []
envi_updated_predicted = []
envi_soil_updated_predicted = []
envi_sludge_predicted = []
envi_soil_sludge_predicted = []

for parent in data_dict:
    for tp in data_dict[parent]["TP_dict"]:
        if package_method_1 in data_dict[parent]["TP_dict"][tp]["source_list"] or package_method_2 in data_dict[parent]["TP_dict"][tp]["source_list"]:
            pps_predicted.append(tp)
        if package_method_5 in data_dict[parent]["TP_dict"][tp]["source_list"]:
            envi_updated_predicted.append(tp) 
        if package_method_6 in data_dict[parent]["TP_dict"][tp]["source_list"]:
            envi_soil_updated_predicted.append(tp)
        if package_method_7 in data_dict[parent]["TP_dict"][tp]["source_list"]:
            envi_sludge_predicted.append(tp)
        if package_method_8 in data_dict[parent]["TP_dict"][tp]["source_list"]:
            envi_soil_sludge_predicted.append(tp)

string_precision = "Precision of methods:" + "\n" + "Prediction Method" + "\t" + "Found" + "\t" + "Predicted" + "\t" + "Precision" + "\n" + package_method_1 + "\t" + str(len(pps))+ "\t" + str(len(pps_predicted))+ "\t" + str(round(len(pps)/len(pps_predicted)*100, 2)) + "%\n"+ package_method_5+ "\t" + str(len(envi_updated))+ "\t"  + str(len(envi_updated_predicted))+ "\t" + str(round(len(envi_updated)/len(envi_updated_predicted)*100, 2)) + "%\n"+ package_method_6 + "\t"+ str(len(envi_soil_updated))+ "\t" + str(len(envi_soil_updated_predicted))+ "\t" + str(round(len(envi_soil_updated)/len(envi_soil_updated_predicted)*100, 2)) + "%\n"+ package_method_7+ "\t" + str(len(envi_sludge))+ "\t" + str(len(envi_sludge_predicted))+ "\t" + str(round(len(envi_sludge)/len(envi_sludge_predicted)*100, 2)) + "%\n"+ package_method_8+ "\t" + str(len(envi_soil_sludge))+ "\t" + str(len(envi_soil_sludge_predicted))+ "\t" + str(round(len(envi_soil_sludge)/len(envi_soil_sludge_predicted)*100, 2)) + "%\n"

with open("./output/precision_of_methods.tsv", 'w') as t:
    t.write(string_precision)

################################################################################################################################################################################################################################################################

# export how many times a subset of methods were used to predict found TPs as tsv

source_dict = {}

for source in source_list_nest:
    key = " ".join(source)
    source_dict[key] = 0

for source in source_list_nest:
    key = " ".join(source)
    source_dict[key] += 1

tp_num = 0

for key in source_dict:
    tp_num += source_dict[key]

string_sources = "Sources of TPs\t" + "Times used\n"

for key in source_dict:
    string_sources += key + "\t"
    string_sources += str(source_dict[key]) + "\n"

with open("./output/sources_of_TPs.tsv", 'w') as s:
    s.write(string_sources)

################################################################################################################################################################################################################################################################

# check how many TPs were for a given parent for each method and export as tsv

string_tp_tsv = "Number of predicted TPs per parent (without considering overlap): \n" + "Parent\t" + package_method_1 + "\t"  + package_method_5 + "\t" + package_method_6 + "\t" + package_method_7 + "\t" + package_method_8 + "\t" + "Total predicted TPs" + "\t" + "Found" + "\t" + "Overall Precision" + "\n"

for parent in data_dict:
    list_tp = []
    list_pps = []
    list_envi_updated = []
    list_envi_soil_updated = []
    list_envi_sludge = []
    list_envi_soil_sludge = []
    
    found_tp_per_parent = []
    
    for tp in data_dict[parent]["TP_dict"]:
        list_tp.append(tp)
        if package_method_1 in data_dict[parent]["TP_dict"][tp]["source_list"] or package_method_2 in data_dict[parent]["TP_dict"][tp]["source_list"]:
            list_pps.append(tp)
        if package_method_5 in data_dict[parent]["TP_dict"][tp]["source_list"]:
            list_envi_updated.append(tp)    
        if package_method_6 in data_dict[parent]["TP_dict"][tp]["source_list"]:
            list_envi_soil_updated.append(tp)
        if package_method_7 in data_dict[parent]["TP_dict"][tp]["source_list"]:
            list_envi_sludge.append(tp)
        if package_method_8 in data_dict[parent]["TP_dict"][tp]["source_list"]:
            list_envi_soil_sludge.append(tp)
            
        if tp in found_tp_smiles_canon:
            found_tp_per_parent.append(tp)
    
                
    string_tp_tsv +=  data_dict[parent]["name"][0] + "\t" + str(len(list_pps)) + "\t"+  str(len(list_envi_updated)) + "\t" + str(len(list_envi_soil_updated)) + "\t"+ str(len(list_envi_sludge)) + "\t"+ str(len(list_envi_soil_sludge)) + "\t"  + str(len(list_tp)) + "\t"  + str(len(found_tp_per_parent)) + "\t" + str(round(len(found_tp_per_parent)/len(list_tp)*100, 2)) + "%" + "\n"
   
    list_tp.clear()
    list_pps.clear()
    list_envi_updated.clear()
    list_envi_soil_updated.clear()
    list_envi_sludge.clear()
    list_envi_soil_sludge.clear()
    
with open("./output/number_of_predicted_TPs_per_parent.tsv", 'w') as h:
    h.write(string_tp_tsv)

################################################################################################################################################################################################################################################################

# get combined probabilities of TPs that were predicted and found and only predicted but not found

predicted_and_found_combined_probability = []
only_predicted_combined_probability = []

for parent in data_dict:
    for tp in data_dict[parent]["TP_dict"]:
        if data_dict.get(parent, {}).get("TP_dict", {}).get(tp, {}).get("combined_prob") is not None:
            if tp in predicted_and_found:
                predicted_and_found_combined_probability += ((data_dict[parent]["TP_dict"][tp]["combined_prob"]))
            else:
                only_predicted_combined_probability += ((data_dict[parent]["TP_dict"][tp]["combined_prob"]))

predicted_and_found_combined_probability_float = []
only_predicted_combined_probability_float = []

for prob in predicted_and_found_combined_probability:
    predicted_and_found_combined_probability_float.append(float(prob))

for prob in only_predicted_combined_probability:
    only_predicted_combined_probability_float.append(float(prob))


do_pickle(predicted_and_found_combined_probability_float, "predicted_and_found_combined_probability.pickle")

do_pickle(only_predicted_combined_probability_float, "only_predicted_combined_probability.pickle")


################################################################################################################################################################################################################################################################

list_all = []

for parent in data_dict:
    for tp in data_dict[parent]["TP_dict"]:
        if package_method_1 in data_dict[parent]["TP_dict"][tp]["source_list"] or package_method_2 in data_dict[parent]["TP_dict"][tp]["source_list"]:
            if package_method_5 in data_dict[parent]["TP_dict"][tp]["source_list"]:
                if package_method_6 in data_dict[parent]["TP_dict"][tp]["source_list"]:
                    if package_method_7 in data_dict[parent]["TP_dict"][tp]["source_list"]:
                        if package_method_8 in data_dict[parent]["TP_dict"][tp]["source_list"]:
                            list_all.append(tp)

# print(list_all)
################################################################################################################################################################################################################################################################

# check which rules were used to predict:

rule_list = []    

for parent in data_dict:
    for tp in data_dict[parent]["TP_dict"]:
        rule = data_dict[parent]["TP_dict"][tp]["rule_list"][0]
        rule_split = rule.split(",")
        rule_list += rule_split


# the bt rule classification is not final, those are only the most important rules
add_O_rules = ["bt0063", "bt0023", "bt0003", "bt0242", "bt0243", "bt0193", "bt0014", "bt0259", "bt0374", "bt0005", "bt0332"]
add_H2O_rules = ["bt0067", "bt0350", "bt0430", "bt0024", "bt0021", "bt0020", "bt0389", "bt0373", "bt0391"]
desat_rules = ["bt0002", "bt0001"]

add_O_list = []
add_H2O_list = []
desat_list = []

for rule in rule_list:
    if rule in add_O_rules:
        add_O_list.append(rule)
    if rule in add_H2O_rules:
        add_H2O_list.append(rule)
    if rule in desat_rules:
        desat_list.append(rule)


print("Oxygen addition: ", len(add_O_list), ", water addition: ", len(add_H2O_list), ", desat: ", len(desat_list))


string_rule = "Rules that were used in the predictions:  need to remove all spaces \" \" and '' signs manually \n"

for rule in rule_list:
    string_rule += rule
    string_rule += "\n"

with open("./output/used_rules.txt", 'w') as p:
    p.write(string_rule)



################################################################################################################################################################################################################################################################

t.close()
h.close()
s.close()

print("Script finished successfully")

################################################################################################################################################################################################################################################################

# ⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⣀⣀⣠⣤⣤⣤⣤⣤⣤⣤⣤⣤⣤⣤⣤⣴⣶⣶⣶⣶⣶⣶⣶⣶⣶⣶⣶⣶⣶⣦⣤⣤⣤⣤⣤⣤⣄⣀⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
# ⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢀⣾⡿⠟⠛⠛⠛⠛⠋⠉⠉⠉⠉⠉⠁⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠈⠉⠉⠉⠉⠛⣿⣦⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
# ⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢸⣿⠁⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢠⣄⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢸⣿⡆⠀⠀⠀⠀⠀⠀⠀⠀⠀
# ⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⣿⣿⠀⠀⣴⡄⠀⠀⠀⠀⣠⡄⠀⠀⠀⠀⠀⠀⠶⠀⠀⠀⠀⠀⠈⠁⠀⠀⠀⠀⠀⠀⠀⠀⠀⢰⣿⠀⠀⠀⠀⣿⡇⠀⠀⠀⠀⠀⠀⠀⠀⠀
# ⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⣿⣿⠀⠀⠉⠁⠀⠀⠀⠘⠟⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠈⠁⠀⠀⠀⢸⣿⡇⠀⠀⠀⠀⠀⠀⠀⠀⠀
# ⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⣿⣿⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⣤⠀⢀⣶⣿⣷⣦⣄⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢸⣿⡇⠀⠀⠀⣀⣤⣄⠀⠀⠀
# ⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⣿⣿⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢀⣤⠀⠀⠀⠀⠀⠀⣸⡇⠀⣸⡿⠀⠀⠉⠻⣿⣦⡀⠀⢰⡿⠀⠀⠀⠀⠀⣸⣿⣁⣴⣾⡿⠟⠛⣿⡄⠀⠀
# ⣴⣿⠿⠿⣿⣶⣦⣄⡀⠀⠀⠀⠀⣿⣿⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠘⠛⠀⠀⠀⠀⠀⠀⠉⠁⠀⣿⡇⠀⠀⠀⠀⠈⠻⣿⣆⠀⠀⠀⠀⠀⠀⠀⣿⣿⣿⠟⠁⠀⠀⠀⣿⡇⠀⠀
# ⢿⣧⠀⠀⠀⠀⠉⠛⢿⣶⣄⠀⠀⣿⣿⠀⠀⠀⠀⠀⠙⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⣀⡀⠀⠀⠀⣿⠇⠀⠀⠀⠀⠀⠀⠈⢻⣷⣤⣤⣤⣤⣤⣼⣿⠟⠁⠀⠀⠀⠀⠀⣿⡇⠀⠀
# ⠈⢿⣧⡀⠀⠀⠀⠀⠀⠈⢻⣷⡄⣿⣿⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠛⠃⠀⠀⢰⣿⡆⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⣿⡇⠀⠀
# ⠀⠈⠻⣷⣄⠀⠀⠀⠀⠀⠀⠙⣿⣿⣿⠀⠀⠀⠀⠀⠀⠀⠀⡿⠇⠀⠀⠀⠀⠀⠀⠀⠀⠀⣴⣿⠟⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠹⣿⣦⡀
# ⠀⠀⠀⠘⢿⣷⣄⠀⠀⠀⠀⠀⠘⣿⣿⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠠⣼⡇⠀⠀⣾⡿⠁⠀⠀⠀⠀⢠⣾⠋⠉⢳⡄⠀⠀⠀⠀⠀⠀⠀⠀⢠⣾⠋⠙⢳⡄⠀⠀⠈⢿⣷
# ⠀⠀⠀⠀⠀⠙⢿⣷⣄⡀⠀⠀⠀⢹⣿⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠿⠁⠀⠀⣿⡇⠀⠀⠀⠀⠀⠸⣿⣶⣤⣾⡇⠀⠀⠀⠀⠀⠀⠀⠀⠸⣿⣧⣤⣾⣿⠀⠀⠀⠘⣿
# ⠀⠀⠀⠀⠀⠀⠀⠈⠻⢿⣦⣄⠀⢸⣿⠀⠀⠀⠀⠀⠀⣾⡇⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⣿⡇⠀⠀⠀⠀⠀⠀⠈⠛⠛⠛⠁⠀⠀⠀⢿⣉⣩⠿⠀⠀⠉⠛⠿⠛⠃⠀⠀⠀⠀⣿
# ⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠙⠿⣿⣾⣿⠀⠀⠀⠀⠀⠀⠉⠁⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢿⣧⠀⠀⠀⠀⠀⠀⠀⠀⠀⡄⠀⠀⠀⠀⠀⢸⡇⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢰⣿
# ⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠈⢻⣿⠀⠀⠀⠀⠀⠀⠀⠀⣀⠀⠀⠀⠀⠀⣠⡄⠀⠀⠘⣿⣧⡀⠀⠀⠀⠀⠀⠀⠀⢷⣤⣀⣀⣀⣴⠟⢿⣤⣀⣀⣀⣴⠇⠀⠀⠀⠀⢠⣿⡟
# ⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⣿⣟⠀⠀⠀⠀⠀⠀⠀⠀⠛⠀⠀⠀⠀⠀⠻⠃⠀⠀⠀⠈⢿⣷⣄⡀⠀⠀⠀⠀⠀⠀⠈⠉⠉⠉⠁⠀⠀⠈⠉⠉⠉⠁⠀⠀⠀⢀⣴⣿⠟⠀
# ⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⣽⣿⣶⣦⣤⣤⣤⣤⣤⣤⣤⣤⣤⣤⣤⣤⣤⣤⣤⣤⣤⣤⣤⣽⣿⣿⣷⣶⣦⣤⣤⣤⣤⣤⣤⣤⣤⣤⣤⣤⣤⣤⣤⣴⣶⣶⣾⠿⠛⠁⠀⠀
# ⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⣼⡟⠀⠈⠉⠉⣻⡟⠛⣻⡿⠛⠛⠛⠛⢿⣿⠿⠛⠛⠛⠛⠛⠛⠛⢻⣿⠏⠉⠉⠉⠉⢻⡟⠛⠛⣻⣿⠋⠉⠉⠙⣿⠉⠉⠉⠀⠀⠀⠀⠀⠀⠀
# ⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⣿⣅⠀⠀⣠⣾⠟⠁⠀⣿⠀⠀⠀⢀⣠⣿⠏⠀⠀⠀⠀⠀⠀⠀⠀⢸⣿⡀⠀⠀⠀⣠⣿⠃⠀⠀⣿⡇⠀⠀⠀⢸⣿⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
# ⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠘⠻⠿⠿⠛⠁⠀⠀⠀⠻⢿⣶⣾⠿⠛⠁⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠛⢿⣷⣶⡿⠟⠁⠀⠀⠀⠻⣷⣄⣠⣴⣿⠃⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
# ⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀   ⠀⠈⠛⠛⠉⠀


