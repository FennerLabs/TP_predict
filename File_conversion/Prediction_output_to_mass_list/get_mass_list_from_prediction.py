# -*- coding: utf-8 -*-
"""
Created on Mon Nov 22 16:17:39 2021
Edited on Tue Aug 17 13:55:32 2022

@author: trostele
"""
# start of script
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
import pubchempy as pcp
import re
import copy
################################################################################################################################################################################################################################################################

######################################## see READme file for instructions! #####################################################################################################################################################################################

################################################################################################################################################################################################################################################################
# INPUT FILES:
                         
# first file (mandatory)
file_location_1 = "input/predictions/Eawag_PPS_BBD_results_batch1.tsv"
prediction_method_1 =  "EAWAG-PPS"            # "EAWAG-PPS" or "envipath" allowed
package_method_1 =  "EAWAG_BBD-PPS_batch1"         # enter used package (is used for naming source in combined data dicts and output)
# second file (optional)
file_location_2 = "input/predictions/TP_prediction_BBD+SOIL_top_50.tsv"
prediction_method_2 = "envipath"
package_method_2 = "enviPath-BBD+SOIL_2"
consider_file_2 = "yes"
# third file (optional)
file_location_3 = "input/predictions/TP_prediction_BBD+SLUDGE_top_50.tsv"
prediction_method_3 = "envipath"
package_method_3 = "enviPath-BBD+SLUDGE_3"
consider_file_3 = "yes"
# fourth file (optional)
file_location_4 = "input/predictions/TP_prediction_BBD+SOIL+SLUDGE_top_50.tsv"
prediction_method_4 = "envipath"
package_method_4 = "enviPath-BBD+SOIL+SLUDGE_4"
consider_file_4 = "yes"
# fifth file (optional)
file_location_5 = "input/predictions/Eawag_PPS_BBD_results_batch2.tsv" # PPS predictions may need to be run in batches
prediction_method_5 = "EAWAG-PPS"
package_method_5 = "EAWAG_BBD-PPS_batch2"
consider_file_5 = "yes"
# sixth file (optional)
file_location_6 = "input/predictions/TP_prediction_BBD_top_50.tsv"
prediction_method_6 = "envipath" 
package_method_6 = "enviPath-BBD_1"
consider_file_6 = "yes"
# seventh file (optional)
file_location_7 = ""
prediction_method_7 =  ""           
package_method_7 = ""
consider_file_7 = "no"
# eighth file (optional)
file_location_8 = ""
prediction_method_8 = ""
package_method_8 = ""
consider_file_8 = "no"                         
# mapping files:  
# code_location (mandatory)
code_location = "./input/code_subst.txt"
# SMILES_location (mandatory)
smi_location = "./input/SMILES_selected_comp.txt"
# name_location (mandatory)
name_location = "./input/name_subst.txt"

# turn search for CAS numbers for all compounds on (True) or off (False)
CAS_search = False

# SCORING SYSTEM:
    
scoring_system = False   # True = active, False = inactive (default)
max_TP_per_parent = 50   # add number of maximal allowed TPs per parent (must be an integer!)

# OUTPUT FILES:

# output could be changed
output_location_1 = "./output/CD_masslist_1_" + package_method_1 + ".csv"
output_location_2 = "./output/CD_masslist_2_" + package_method_2 + ".csv"
output_location_3 = "./output/CD_masslist_3_" + package_method_3 + ".csv"
output_location_4 = "./output/CD_masslist_4_" + package_method_4 + ".csv"
output_location_5 = "./output/CD_masslist_5_" + package_method_5 + ".csv"
output_location_6 = "./output/CD_masslist_6_" + package_method_6 + ".csv"
output_location_7 = "./output/CD_masslist_7_" + package_method_7 + ".csv"
output_location_8 = "./output/CD_masslist_8_" + package_method_8 + ".csv"
# output: combined mass list of all methods
output_file_CD_masslist = "./output/CD_masslist_combined.csv"
# output: combined data of all methods
output_file_all_data = "./output/combined_overview.csv"
# output: inclusion list for QExactivePlus for positive mode
output_inclusion_pos = "./output/inclusion_list_pos.csv"
# output: inclusion list for QExactivePlus for negative mode
output_inclusion_neg = "./output/inclusion_list_neg.csv"
# output: max element count file
output_max_element = "./output/max_element_count.txt"
# output: removed TPs above 100 u
output_removed_tps = "./output/removed_TPs_above_100_u.csv"

################################################################################################################################################################################################################################################################
################################################################################################################################################################################################################################################################
# FUNCTIONS

def file_to_csv(input_file, pickle_file, csv_file, file_type):
    """

    :param input_file:
    :param pickle_file:
    :param csv_file:
    :param file_type: 'envipath' or 'EAWAG-PPS'
    """
    print("Converting {} file".format(input_file))
    if file_type == 'envipath':
        D1 = read_enviPath_file_to_dict(input_file)
    elif file_type == 'EAWAG-PPS':
        D1 = read_PPS_file_to_dict(input_file)
    else:
        raise ValueError("Possible values for file_type: 'envipath' or 'EAWAG-PPS'")
    D2 = canonicalize_dict(D1)
    D3 = annotate_dict(D2, file_type)
    do_pickle(D3, pickle_file)
    dict_to_csv(D3, csv_file)
    return D3

def read_enviPath_file_to_dict(input_file):
    envipath_file = open(input_file)
    sep = '\t'
    slash = "///" 
    smiles = "SMILES"
    data_dict_envi = {}
    for line in envipath_file:
        linelist_envi = line.rstrip().split(sep) 
        if line.startswith(slash): #skip pathway line
            substance = ""
            continue
        else:
            if line.startswith(smiles): #skip header line
                continue
            else:
                if len(linelist_envi) == 6:#skip parent line because TP_1 is always first generation (except Atv)
                    substance = linelist_envi[0]
                    continue
                else:
                    if data_dict_envi.get(substance):  #if parent SMILES exists as key then append list
                        data_dict_envi[substance]['TP_list'].append(linelist_envi[0])
                        data_dict_envi[substance]['bt_list'].append(linelist_envi[3])
                        data_dict_envi[substance]['code_TP'].append(linelist_envi[1])
                        data_dict_envi[substance]['combined_prob'].append(linelist_envi[2])
                    else: #otherwise create new entry into dict
                        data_dict_envi[substance] = {'TP_list': [linelist_envi[0]], "TP_list_canon_2":[],"TP_list_canon":[],'combined_prob': [linelist_envi[2]], 'bt_list': [linelist_envi[3]], 'code_TP': [linelist_envi[1]], "code_parent":[], "name" : [], "ID_TP": [], "ID_parent":None, "mass_TP" : [], "mass_parent" : [], "Structure_TP": [], "Structure_parent": None, "Formula_TP":[], "Formula_parent":None, "inchi_TP":[], "inchi_parent": None}
    for parent in code_dict:
        if parent not in data_dict_envi:
            data_dict_envi[parent] = {'TP_list': [], "TP_list_canon": [], "TP_list_canon_2": [],
                                            'bt_list': [], 'code_TP': [], "code_parent": [], "name": [],
                                            "ID_TP": [], "ID_parent": None, "mass_TP": [], "mass_parent": [],
                                            "Structure_TP": [], "Structure_parent": None, "Formula_TP": [],
                                            "Formula_parent": None, "inchi_TP": [], "inchi_parent": None}
    for key in data_dict_envi: 
        if key in code_dict.keys():
            data_dict_envi[key]['code_parent'].append(code_dict.get(key))  # add the code to the data dict
        if key in name_dict.keys(): 
            data_dict_envi[key]["name"].append(name_dict.get(key)) # add name to data_dict
    return data_dict_envi

def canonicalize_smiles(smiles):
    mol = Chem.MolFromSmiles(smiles)  # creates mol object from SMILES
    uncharger = rdMolStandardize.Uncharger() # easier to access
    uncharged = uncharger.uncharge(mol)  # protonates or deprotonates the mol object
    new_smiles = rdkit.Chem.rdmolfiles.MolToSmiles(uncharged)  # converts mol object to canonical SMILES
    can_smiles = Chem.CanonSmiles(new_smiles)
    return can_smiles

def max_element_count(smi_list):
    max_C = 0
    max_N = 0
    max_F = 0
    max_O = 0
    max_S = 0
    max_P = 0
    max_Cl = 0
    max_Br = 0
    max_I = 0
    max_H = 0
    
    for smi in smi_list:
        if smi.count("C") > max_C:
            max_C = smi.count("C")
        if smi.count("N") > max_N:
            max_N = smi.count("N")
        if smi.count("F") > max_F:
            max_F = smi.count("F")
        if smi.count("O") > max_O:
            max_O = smi.count("O")        
        if smi.count("S") > max_S:
            max_S = smi.count("S")
        if smi.count("P") > max_P:
            max_P = smi.count("P")
        if smi.count("Cl") > max_Cl:
            max_Cl = smi.count("Cl")        
        if smi.count("Br") > max_Br:
            max_Br = smi.count("Br")
        if smi.count("I") > max_I:
            max_I = smi.count("I")
        mol = Chem.MolFromSmiles(smi)
        mol = Chem.AddHs(mol)
        mol = Chem.MolToSmiles(mol, allHsExplicit=True)
        hcount = mol.count("H")
        if hcount > max_H:
            max_H = hcount
            
    with open(output_max_element, 'w') as f:
        f.write("max. element count: " + "max C = " + str(max_C) + ", max H = " + str(max_H) + ", max O = " + str(max_O) + ", max N = " + str(max_N)
          + ", max S = " + str(max_S) + ", max P = " + str(max_P) + ", max Cl = " + str(max_Cl) + ", max Br = " + str(max_Br) + ", max I = " + str(max_I))        
    print("max. element count: " + "max C = " + str(max_C) + ", max H = " + str(max_H) + ", max O = " + str(max_O) + ", max N = " + str(max_N)
          + ", max S = " + str(max_S) + ", max P = " + str(max_P) + ", max Cl = " + str(max_Cl) + ", max Br = " + str(max_Br) + ", max I = " + str(max_I))
    return

def suggest_stepped_nce(smi_list):
    MolWeight_list = []
    for compound in smi_list:
            MolWeight_list.append(Chem.Descriptors.ExactMolWt(Chem.MolFromSmiles(compound)))        
    nce_list = []
    for mass in MolWeight_list:
            if mass > 350:
                nce_list.append(15)
            else:
                nce_calc = 5 * round(((mass * -0.41) + 160)/5)
                nce_list.append(nce_calc)
                
    max_nce = max(nce_list)
    if max_nce > 120:
        high_nce = 100
    else:
        high_nce = max_nce - 20
    avg_nce = sum(nce_list)/len(nce_list)
    middle_nce = (5 * round(avg_nce/5)) - 5
    min_nce = min(nce_list)
    if min_nce == 15:
        low_nce = min_nce
    else:
        low_nce = min_nce - 5
    
    if high_nce < 0 or middle_nce - low_nce < 10 or high_nce - middle_nce < 10:
        print("Stepped NCE approach not recommended")
    else:
        print("Suggested Stepped NCE: " + "Low NCE = " + str(low_nce) + ", Middle NCE = " + str(middle_nce) + ", High NCE = " + str(high_nce))
    
    return

def do_pickle(d, pickle_file):
    with open("./output/" + pickle_file, 'wb') as handle:
        pickle.dump(d, handle, protocol=pickle.HIGHEST_PROTOCOL)

def get_pickle(pickle_file):
    with open("./output/" + pickle_file, 'rb') as handle:
        d = pickle.load(handle)
    return d

def canonicalize_dict(D):
    new_D = {}
    for parent in D:
        tp_dict = D[parent]
        tp_dict["TP_list_canon"] = []
        for tp in D[parent]["TP_list"]:
            tp_dict["TP_list_canon"].append(canonicalize_smiles(tp))
        new_D[canonicalize_smiles(parent)] = tp_dict
    return new_D

def annotate_dict(data_dict, data_type):
    for parent in data_dict:
        counter = 1
        for tp in data_dict[parent]["TP_list_canon"]:
            data_dict[parent]['mass_TP'].append(Chem.Descriptors.ExactMolWt(Chem.MolFromSmiles(tp)))
            # add structure of TP to dict
            data_dict[parent]["Structure_TP"].append(Chem.MolFromSmiles(tp))
            # add molecular formula of TP
            data_dict[parent]["Formula_TP"].append(CalcMolFormula(Chem.MolFromSmiles(tp)))
            # add inchikey of TP
            data_dict[parent]["inchi_TP"].append(Chem.inchi.MolToInchiKey(Chem.MolFromSmiles(tp)))
            if data_type == 'EAWAG-PPS': # used for PPS where TPs are not named automatically
                data_dict[parent]['code_TP'].append("TP_" + data_dict[parent]['code_parent'][0] + "_" + str(counter))
                counter += 1
        # add parent mass
        data_dict[parent]['mass_parent'].append(Chem.Descriptors.ExactMolWt(Chem.MolFromSmiles(parent)))
        # add structure of parent to dict
        data_dict[parent]["Structure_parent"] = Chem.MolFromSmiles(parent)
        # add molecular formula of parent
        data_dict[parent]["Formula_parent"] = CalcMolFormula(Chem.MolFromSmiles(parent))
        # add inchikey of parent
        data_dict[parent]["inchi_parent"] = Chem.inchi.MolToInchiKey(Chem.MolFromSmiles(parent))
        
        # !!!
        # add CAS of parent from inchi key, if there is no CAS number for TP then there still needs to be an entry so that later list has same length! 
        if CAS_search == True:
            parent_cas = get_cas_inchi(data_dict[parent]["inchi_parent"])
        if CAS_search == False:
            parent_cas = "CAS search disabled"
        if len(parent_cas) > 0: # only add CAS if it was found
            data_dict[parent]["ID_parent"] = parent_cas          
        else: #otherwise add empty string
            data_dict[parent]["ID_parent"] = ""
        # add CAS of TPs from inchi key
        for tp in data_dict[parent]["inchi_TP"]:
            if CAS_search == True:
                tp_cas = get_cas_inchi(tp)
            if CAS_search == False:
                tp_cas = "CAS search disabled"
            if len(tp_cas) > 0:
                data_dict[parent]["ID_TP"].append(tp_cas)
            else:
                data_dict[parent]["ID_TP"].append("")
        del data_dict[parent]["TP_list"]

    return data_dict

def dict_to_csv(data_dict, output_file):
    # create lists for each column and then create dict with correct layout which can be converted to dataframe using pandas
    df_name_list = []
    df_ID_list = []
    df_Formula_list = []
    df_MolWeight_list = []
    df_Structure_list = []
    # add all the data to the lists
    for parent in data_dict:
        df_name_list.append(data_dict[parent]["code_parent"][0]) #append code of parent first
        for tp in data_dict[parent]["code_TP"]: #then add all the codes of the TPs
            df_name_list.append(tp)
        # add CAS number
        df_ID_list.append(data_dict[parent]["ID_parent"])
        for tp in data_dict[parent]["ID_TP"]:
            df_ID_list.append(tp)
        # add monoisotopic mass
        df_MolWeight_list.append(data_dict[parent]["mass_parent"][0])
        for tp in data_dict[parent]["mass_TP"]:
            df_MolWeight_list.append(tp)
        # add chemical formula
        df_Formula_list.append(data_dict[parent]["Formula_parent"])
        for tp in data_dict[parent]["Formula_TP"]:
            df_Formula_list.append(tp)
        # add mol file
        mol_rep_1 = (Chem.MolToMolBlock(data_dict[parent]["Structure_parent"])).replace("\n", ";") #replace newline character with semicolon
        mol_rep_2 = mol_rep_1[6:] #skip first few spaces
        mol_rep_3 = mol_rep_2[:-1] #remove the last semicolon
        df_Structure_list.append(mol_rep_3) #add string to the list
        for tp in data_dict[parent]["Structure_TP"]:
            tp_mol_rep_1 = (Chem.MolToMolBlock(tp)).replace("\n",";")
            tp_mol_rep_2 = tp_mol_rep_1[6:]
            tp_mol_rep_3 = tp_mol_rep_2[:-1]
            df_Structure_list.append(tp_mol_rep_3)   
    # all lists must be the same length to convert it to dataframe
    assert len(df_name_list) == len(df_Formula_list) == len(df_MolWeight_list) == len(df_Structure_list) == len(df_ID_list), "Error: all lists must be the same length to convert it to dataframe"
    #create dict and convert to dataframe
    df_dict = {"Name": df_name_list, "ID": df_ID_list,"Formula": df_Formula_list,"MolWeight": df_MolWeight_list, "Structure": df_Structure_list}
    df = pd.DataFrame.from_dict(df_dict)
    # export dataframe as csv
    df.to_csv(output_file, index = False, sep = "\t")

def combined_dict_to_csv(data_dict, output_file_CD, output_file_complete):
    # create lists for each column and then create dict with correct which can be converted to dataframe using pandas
    df_name_list_com = []
    df_smi_list_com = []
    df_ID_list_com = []
    df_Formula_list_com = []
    df_MolWeight_list_com = []
    df_Structure_list_com = []
    df_name_parent_list_com = []
    df_inchikey_list_com = []
    df_score_list_com = []
    df_rules_list_com = []
    df_source_list_com = []
    df_alt_parent_list_com = []   
    # change the codes of the TPs, so each TP has its own name
    for parent in data_dict:
        counter = 1
        for tp in data_dict[parent]["TP_dict"]:
            data_dict[parent]["TP_dict"][tp]["code"] = "TP_" + data_dict[parent]['code_parent'][0] + "_" + str(counter)
            counter = counter + 1 
    # add all the data to the lists
    for parent in data_dict:
        df_name_list_com.append(data_dict[parent]["code_parent"][0])
        df_smi_list_com.append(parent)
        df_ID_list_com.append(data_dict[parent]["ID_parent"])
        df_MolWeight_list_com.append(data_dict[parent]["mass_parent"][0])
        df_Formula_list_com.append(data_dict[parent]["Formula_parent"])
        mol_rep_1 = (Chem.MolToMolBlock(data_dict[parent]["Structure_parent"])).replace("\n", ";") #replace newline character with semicolon
        mol_rep_2 = mol_rep_1[6:] #skip first few spaces
        mol_rep_3 = mol_rep_2[:-1] #remove the last semicolon
        df_Structure_list_com.append(mol_rep_3) #add string to the list  
        df_name_parent_list_com.append(data_dict[parent]["name"][0])
        df_inchikey_list_com.append(data_dict[parent]["inchi_parent"])
        df_score_list_com.append("100")
        df_rules_list_com.append("")
        df_source_list_com.append("")
        df_alt_parent_list_com.append([])
        for tp in data_dict[parent]["TP_dict"]: 
            df_name_list_com.append(data_dict[parent]["TP_dict"][tp]["code"])
            df_smi_list_com.append(tp)
            df_ID_list_com.append(data_dict[parent]["TP_dict"][tp]["CAS"])
            df_MolWeight_list_com.append(data_dict[parent]["TP_dict"][tp]["mass"])
            df_Formula_list_com.append(data_dict[parent]["TP_dict"][tp]["Formula"])
            tp_mol_rep_1 = (Chem.MolToMolBlock(data_dict[parent]["TP_dict"][tp]["Structure"])).replace("\n",";")
            tp_mol_rep_2 = tp_mol_rep_1[6:]
            tp_mol_rep_3 = tp_mol_rep_2[:-1]
            df_Structure_list_com.append(tp_mol_rep_3)
            df_name_parent_list_com.append(data_dict[parent]["name"][0])
            df_inchikey_list_com.append(data_dict[parent]["TP_dict"][tp]["InchiKey"])
            df_score_list_com.append(data_dict[parent]["TP_dict"][tp]["score"])
            df_rules_list_com.append(data_dict[parent]["TP_dict"][tp]["rule_list"]) 
            df_source_list_com.append(data_dict[parent]["TP_dict"][tp]["source_list"])
            df_alt_parent_list_com.append(data_dict[parent]["TP_dict"][tp]["alternative_parent"])        
    
    max_element_count(df_smi_list_com)
    
    # all lists must be the same length to convert it to dataframe
    assert len(df_name_list_com) == len(df_Formula_list_com) == len(df_MolWeight_list_com) == len(df_Structure_list_com) == len(df_ID_list_com), "Error: all lists must be the same length to convert it to dataframe"

    #create dict and convert to dataframe
    df_com_dict = {"Name": df_name_list_com, "ID": df_ID_list_com,"Formula": df_Formula_list_com,"MolWeight": df_MolWeight_list_com, "Structure": df_Structure_list_com}
    df_com = pd.DataFrame.from_dict(df_com_dict)
    # export dataframe as csv
    df_com.to_csv(output_file_CD, index = False, sep = "\t")
    df_complete_dict = {"Name of parent":df_name_parent_list_com,"SMILES": df_smi_list_com,
                   "Name": df_name_list_com, "Score": df_score_list_com,
                   "Source": df_source_list_com,
                   "CAS": df_ID_list_com,"Formula": df_Formula_list_com,
                   "MolWeight": df_MolWeight_list_com, "InchiKey":df_inchikey_list_com, "Alternative parent": df_alt_parent_list_com, "bt rules": df_rules_list_com}
    df_complete = pd.DataFrame.from_dict(df_complete_dict)
    # export dataframe as csv
    df_complete.to_csv(output_file_complete, index = False, sep = ",")
    
    # create inclusion list for QExactivePlus
    m_proton = 1.0072756
    df_M_plus_H = []
    df_M_minus_H = []
    df_polarity_pos = []
    df_polarity_neg = []
    df_empty = []
    df_nce_type = []
    df_nce = []
    for mass in df_MolWeight_list_com:
        df_M_plus_H.append(mass + m_proton)
        df_M_minus_H.append(mass - m_proton)
        df_polarity_pos.append("Positive")
        df_polarity_neg.append("Negative")
        df_empty.append(" ")
        df_nce_type.append("NCE")
        if mass > 350:
            df_nce.append(15)
        else:
            nce_calc = 5 * round(((mass * -0.41) + 160)/5)
            df_nce.append(nce_calc)
                   
    d_inclusion_pos = {"Mass [m/z]": df_M_plus_H ,"Formula [M]": df_empty,
                   "Formula type": df_empty, "Species": df_empty,
                   "CS [z]": df_empty, "Polarity": df_polarity_pos, "Start [min]": df_empty, "End [min]": df_empty,
                   "(N)CE": df_nce, "(N)CE type":df_nce_type, "MSX ID": df_empty, "Comment": df_name_list_com}
    df_inclusion_pos = pd.DataFrame.from_dict(d_inclusion_pos)
    # export dataframe as csv
    df_inclusion_pos.to_csv(output_inclusion_pos, index = False, sep = ",")  
    d_inclusion_neg = {"Mass [m/z]": df_M_minus_H ,"Formula [M]": df_empty,
                   "Formula type": df_empty, "Species": df_empty,
                   "CS [z]": df_empty, "Polarity": df_polarity_neg, "Start [min]": df_empty, "End [min]": df_empty,
                   "(N)CE": df_nce, "(N)CE type":df_nce_type, "MSX ID": df_empty, "Comment": df_name_list_com}
    df_inclusion_neg = pd.DataFrame.from_dict(d_inclusion_neg)
    # export dataframe as csv
    df_inclusion_neg.to_csv(output_inclusion_neg, index = False, sep = ",") 
    
    suggest_stepped_nce(df_smi_list_com)
    
    print("Export complete")

def get_cas_inchi(inchi): # add get CAS from inchikey function
    cas_rns = []
    inchi_split = inchi.split("-")[0]
    results = pcp.get_synonyms(inchi_split, 'inchikey')
    for result in results:
        for syn in result.get('Synonym', []):
            match = re.match('(\d{2,7}-\d\d-\d)', syn)
            if match:
                cas_rns.append(match.group(1))
    return cas_rns

def read_PPS_file_to_dict(PPS_file_location):
    PPS_file = open(PPS_file_location)
    sep = '\t'
    line_1 = PPS_file.readline()
    line_list_1 = line_1.rstrip().split(sep)  # rstrip() removes the newline character '\n' at the end of the file
    Settings = {}
    Settings[line_list_1[0]] = line_list_1[1]
    line_2 = PPS_file.readline()
    line_list_2 = line_2.rstrip().split(sep)
    Settings[line_list_2[0]] = line_list_2[1]
    line_3 = PPS_file.readline()
    line_4 = PPS_file.readline()
    PPS_file.readline()
    header_line = PPS_file.readline()
    compound_list = header_line.rstrip().split(sep)
    data = {}
    data_dict = {}
    for line in PPS_file:
        linelist = line.rstrip().split(sep)  # we get the list
        # we know that the first item of the list is the TP, and the following items are biotransformation rules producing the TP from a given compound
        for index, substance in enumerate(compound_list):
            # The first item in the compound list is empty, so let's skip that
            if index == 0:
                continue
            # empty fields at the end of the line are not imported as empty strings, add them manually
            while len(linelist) < len(compound_list):
                linelist.append('')
            set_of_rules = linelist[index]      
            if set_of_rules != '': # Only print if the set of rules is not empty
                if data_dict.get(substance):  # if key exist then append list
                    data_dict[substance]['TP_list'].append(linelist[0])
                    data_dict[substance]['bt_list'].append(set_of_rules)
                else:  # otherwise create new entry into dict
                    data_dict[substance] = {'TP_list': [linelist[0]], "TP_list_canon": [], "TP_list_canon_2": [],
                                            'bt_list': [set_of_rules], 'code_TP': [], "code_parent": [], "name": [],
                                            "ID_TP": [], "ID_parent": None, "mass_TP": [], "mass_parent": [],
                                            "Structure_TP": [], "Structure_parent": None, "Formula_TP": [],
                                            "Formula_parent": None, "inchi_TP": [], "inchi_parent": None}
    for key in data_dict:  
        if key in code_dict.keys():
            data_dict[key]['code_parent'].append(code_dict.get(key))  #add the code to the data dict
        if key in name_dict.keys():
            data_dict[key]["name"].append(name_dict.get(key))   # add name to data_dict like the code before
    for parent in code_dict and name_dict:
        if parent not in data_dict.keys():
            data_dict[parent] = {'TP_list': [], "TP_list_canon": [], "TP_list_canon_2": [],
                                            'bt_list': [], 'code_TP': [], "code_parent": [code_dict.get(parent)], "name": [name_dict.get(parent)],
                                            "ID_TP": [], "ID_parent": None, "mass_TP": [], "mass_parent": [],
                                            "Structure_TP": [], "Structure_parent": None, "Formula_TP": [],
                                            "Formula_parent": None, "inchi_TP": [], "inchi_parent": None}
    return data_dict

def combine_dict (d_1, method_1_package, d_2, method_2_package, d_3, method_3_package, d_4, method_4_package, d_5, method_5_package, d_6, method_6_package, d_7, method_7_package, d_8, method_8_package):
    print("Combining "+ method_1_package + ", " +  method_2_package + ", " +  method_3_package + " , " +  method_4_package + " , " +  method_5_package+ " , " +  method_6_package + " , " +  method_7_package+ " and " +  method_8_package)
    # create new dict from copy of envi dict and then delete TP info
    data_dict_com = copy.deepcopy(d_1) # need deepcopy otherwise it will still change original dict   
    for parent in data_dict_com:
        if data_dict_com.get(parent, {}).get("combined_prob") is None:
            if data_dict_com.get(parent, {}).get("bt_list") is not None:
                del data_dict_com[parent]["bt_list"]
                del data_dict_com[parent]["code_TP"]
                del data_dict_com[parent]["inchi_TP"]
                del data_dict_com[parent]["mass_TP"]
                del data_dict_com[parent]["Formula_TP"]
                del data_dict_com[parent]["Structure_TP"]
                del data_dict_com[parent]["ID_TP"]
                del data_dict_com[parent]["TP_list_canon"]
                del data_dict_com[parent]["TP_list_canon_2"]         
        else:
            if data_dict_com.get(parent, {}).get("bt_list") is not None:
                del data_dict_com[parent]["bt_list"]
                del data_dict_com[parent]["code_TP"]
                del data_dict_com[parent]["inchi_TP"]
                del data_dict_com[parent]["mass_TP"]
                del data_dict_com[parent]["Formula_TP"]
                del data_dict_com[parent]["Structure_TP"]
                del data_dict_com[parent]["ID_TP"]
                del data_dict_com[parent]["TP_list_canon"]
                del data_dict_com[parent]["combined_prob"]
                del data_dict_com[parent]["TP_list_canon_2"]

    # add TP data from first data dict    
    for parent in d_1:
        if d_1.get(parent, {}).get("combined_prob") is not None:
            for index, tp in enumerate(d_1[parent]["TP_list_canon"]): 
                if data_dict_com.get(parent, {}).get("TP_dict") is None:
                    data_dict_com[parent]["TP_dict"] = {}  
                    data_dict_com[parent]["TP_dict"][tp] = {"CAS": d_1[parent]["ID_TP"][index],
                                                                      "rule_list":[d_1[parent]["bt_list"][index]],
                                                                      "mass": d_1[parent]["mass_TP"][index],
                                                                      "Formula": d_1[parent]["Formula_TP"][index],
                                                                      "source_list" :[method_1_package],
                                                                      "code": d_1[parent]["code_TP"][index],
                                                                      "Structure" : d_1[parent]["Structure_TP"][index],
                                                                      "combined_prob": [d_1[parent]["combined_prob"][index]],
                                                                      "score": 100, "InchiKey": d_1[parent]["inchi_TP"][index],
                                                                      "alternative_parent" : []}
                else:
                    if data_dict_com[parent]["TP_dict"].get(tp): # if TP already there then append list
                        if method_1_package not in data_dict_com[parent]["TP_dict"][tp]["source_list"]:
                            data_dict_com[parent]["TP_dict"][tp]["source_list"].append(method_1_package)
                            data_dict_com[parent]["TP_dict"][tp]["combined_prob"].append(d_1[parent]["combined_prob"][index])
                    else:
                        data_dict_com[parent]["TP_dict"][tp] = {"CAS": d_1[parent]["ID_TP"][index],
                                                                      "rule_list":[d_1[parent]["bt_list"][index]],
                                                                      "mass": d_1[parent]["mass_TP"][index],
                                                                      "Formula": d_1[parent]["Formula_TP"][index],
                                                                      "source_list" :[method_1_package],
                                                                      "code": d_1[parent]["code_TP"][index],
                                                                      "Structure" : d_1[parent]["Structure_TP"][index],
                                                                      "combined_prob": [d_1[parent]["combined_prob"][index]],
                                                                      "score": 100, "InchiKey": d_1[parent]["inchi_TP"][index],
                                                                      "alternative_parent" : []}

        else:
            for parent in d_1:
                for index, tp in enumerate(d_1[parent]["TP_list_canon"]):
                    if data_dict_com.get(parent, {}).get("TP_dict") is None:
                        data_dict_com[parent]["TP_dict"] = {}
                        data_dict_com[parent]["TP_dict"][tp] = {"CAS": d_1[parent]["ID_TP"][index],
                                                                      "rule_list":[d_1[parent]["bt_list"][index]],
                                                                      "mass": d_1[parent]["mass_TP"][index],
                                                                      "Formula": d_1[parent]["Formula_TP"][index],
                                                                      "source_list" :[method_1_package],
                                                                      "code": d_1[parent]["code_TP"][index],
                                                                      "Structure" : d_1[parent]["Structure_TP"][index],
                                                                      "combined_prob": [],
                                                                      "score": 100, "InchiKey": d_1[parent]["inchi_TP"][index],
                                                                      "alternative_parent" : []}
                    else:
                        if data_dict_com[parent]["TP_dict"].get(tp): # if TP already there then append list
                            if method_1_package not in data_dict_com[parent]["TP_dict"][tp]["source_list"]:
                                data_dict_com[parent]["TP_dict"][tp]["source_list"].append(method_1_package)
                        else:
                            data_dict_com[parent]["TP_dict"][tp] = {"CAS": d_1[parent]["ID_TP"][index],
                                                                              "rule_list":[d_1[parent]["bt_list"][index]],
                                                                              "mass": d_1[parent]["mass_TP"][index],
                                                                              "Formula": d_1[parent]["Formula_TP"][index],
                                                                              "source_list" :[method_1_package],
                                                                              "code": d_1[parent]["code_TP"][index],
                                                                              "Structure" : d_1[parent]["Structure_TP"][index],
                                                                              "combined_prob": [],
                                                                              "score": 100, "InchiKey": d_1[parent]["inchi_TP"][index],
                                                                              "alternative_parent" : []}
       
    # add TP data from second data dict
    for parent in d_2:
        if d_2.get(parent, {}).get("combined_prob") is not None:
            for index, tp in enumerate(d_2[parent]["TP_list_canon"]): 
                if data_dict_com.get(parent, {}).get("TP_dict") is None:
                    data_dict_com[parent]["TP_dict"] = {}
                    data_dict_com[parent]["TP_dict"][tp] = {"CAS": d_2[parent]["ID_TP"][index],
                                                                      "rule_list":[d_2[parent]["bt_list"][index]],
                                                                      "mass": d_2[parent]["mass_TP"][index],
                                                                      "Formula": d_2[parent]["Formula_TP"][index],
                                                                      "source_list" :[method_2_package],
                                                                      "code": d_2[parent]["code_TP"][index],
                                                                      "Structure" : d_2[parent]["Structure_TP"][index],
                                                                      "combined_prob": [d_2[parent]["combined_prob"][index]],
                                                                      "score": 100, "InchiKey": d_2[parent]["inchi_TP"][index],
                                                                      "alternative_parent" : []}
                else:
                    if data_dict_com[parent]["TP_dict"].get(tp): # if TP already there then append list
                        data_dict_com[parent]["TP_dict"][tp]["combined_prob"].append(d_2[parent]["combined_prob"][index])
                        if d_2[parent]["bt_list"][index] not in data_dict_com[parent]["TP_dict"][tp]["rule_list"]:
                            data_dict_com[parent]["TP_dict"][tp]["rule_list"].append(d_2[parent]["bt_list"][index])
                        if method_2_package not in data_dict_com[parent]["TP_dict"][tp]["source_list"]:
                            data_dict_com[parent]["TP_dict"][tp]["source_list"].append(method_2_package)
                    else:
                        data_dict_com[parent]["TP_dict"][tp] = {"CAS": d_2[parent]["ID_TP"][index],
                                                                      "rule_list":[d_2[parent]["bt_list"][index]],
                                                                      "mass": d_2[parent]["mass_TP"][index],
                                                                      "Formula": d_2[parent]["Formula_TP"][index],
                                                                      "source_list" :[method_2_package],
                                                                      "code": d_2[parent]["code_TP"][index],
                                                                      "Structure" : d_2[parent]["Structure_TP"][index],
                                                                      "combined_prob": [d_2[parent]["combined_prob"][index]],
                                                                      "score": 100, "InchiKey": d_2[parent]["inchi_TP"][index],
                                                                      "alternative_parent" : []}

        else:
            for parent in d_2:
                for index, tp in enumerate(d_2[parent]["TP_list_canon"]):
                    if data_dict_com.get(parent, {}).get("TP_dict") is None:
                        data_dict_com[parent]["TP_dict"] = {}
                        data_dict_com[parent]["TP_dict"][tp] = {"CAS": d_2[parent]["ID_TP"][index],
                                                                      "rule_list":[d_2[parent]["bt_list"][index]],
                                                                      "mass": d_2[parent]["mass_TP"][index],
                                                                      "Formula": d_2[parent]["Formula_TP"][index],
                                                                      "source_list" :[method_2_package],
                                                                      "code": d_2[parent]["code_TP"][index],
                                                                      "Structure" : d_2[parent]["Structure_TP"][index],
                                                                      "combined_prob": [],
                                                                      "score": 100, "InchiKey": d_2[parent]["inchi_TP"][index],
                                                                      "alternative_parent" : []}
                    else:
                        if data_dict_com[parent]["TP_dict"].get(tp): # if TP already there then append list
                            if d_2[parent]["bt_list"][index] not in data_dict_com[parent]["TP_dict"][tp]["rule_list"]:
                                data_dict_com[parent]["TP_dict"][tp]["rule_list"].append(d_2[parent]["bt_list"][index])
                            if method_2_package not in data_dict_com[parent]["TP_dict"][tp]["source_list"]:
                                data_dict_com[parent]["TP_dict"][tp]["source_list"].append(method_2_package)
                        else:
                            data_dict_com[parent]["TP_dict"][tp] = {"CAS": d_2[parent]["ID_TP"][index],
                                                                              "rule_list":[d_2[parent]["bt_list"][index]],
                                                                              "mass": d_2[parent]["mass_TP"][index],
                                                                              "Formula": d_2[parent]["Formula_TP"][index],
                                                                              "source_list" :[method_2_package],
                                                                              "code": d_2[parent]["code_TP"][index],
                                                                              "Structure" : d_2[parent]["Structure_TP"][index],
                                                                              "combined_prob": [],
                                                                              "score": 100, "InchiKey": d_2[parent]["inchi_TP"][index],
                                                                              "alternative_parent" : []}
                        
     # add TP data from third data dict (if specified)                       
    if d_3 != "none":
        for parent in d_3:
            if d_3.get(parent, {}).get("combined_prob") is not None:
                for index, tp in enumerate(d_3[parent]["TP_list_canon"]): 
                    if data_dict_com.get(parent, {}).get("TP_dict") is None:
                        data_dict_com[parent]["TP_dict"] = {}
                        data_dict_com[parent]["TP_dict"][tp] = {"CAS": d_3[parent]["ID_TP"][index],
                                                                      "rule_list":[d_3[parent]["bt_list"][index]],
                                                                      "mass": d_3[parent]["mass_TP"][index],
                                                                      "Formula": d_3[parent]["Formula_TP"][index],
                                                                      "source_list" :[method_3_package],
                                                                      "code": d_3[parent]["code_TP"][index],
                                                                      "Structure" : d_3[parent]["Structure_TP"][index],
                                                                      "combined_prob": [d_3[parent]["combined_prob"][index]],
                                                                      "score": 100, "InchiKey": d_3[parent]["inchi_TP"][index],
                                                                      "alternative_parent" : []}
                    else:
                        if data_dict_com[parent]["TP_dict"].get(tp): # if TP already there then append list
                            data_dict_com[parent]["TP_dict"][tp]["combined_prob"].append(d_3[parent]["combined_prob"][index])
                            if d_3[parent]["bt_list"][index] not in data_dict_com[parent]["TP_dict"][tp]["rule_list"]:
                                data_dict_com[parent]["TP_dict"][tp]["rule_list"].append(d_3[parent]["bt_list"][index])
                            if method_3_package not in data_dict_com[parent]["TP_dict"][tp]["source_list"]:
                                data_dict_com[parent]["TP_dict"][tp]["source_list"].append(method_3_package)
                        else:
                            data_dict_com[parent]["TP_dict"][tp] = {"CAS": d_3[parent]["ID_TP"][index],
                                                                          "rule_list":[d_3[parent]["bt_list"][index]],
                                                                          "mass": d_3[parent]["mass_TP"][index],
                                                                          "Formula": d_3[parent]["Formula_TP"][index],
                                                                          "source_list" :[method_3_package],
                                                                          "code": d_3[parent]["code_TP"][index],
                                                                          "Structure" : d_3[parent]["Structure_TP"][index],
                                                                          "combined_prob": [d_3[parent]["combined_prob"][index]],
                                                                          "score": 100, "InchiKey": d_3[parent]["inchi_TP"][index],
                                                                          "alternative_parent" : []}
    
            else:
                for parent in d_3:
                    for index, tp in enumerate(d_3[parent]["TP_list_canon"]):
                        if data_dict_com.get(parent, {}).get("TP_dict") is None:
                            data_dict_com[parent]["TP_dict"] = {}
                            data_dict_com[parent]["TP_dict"][tp] = {"CAS": d_3[parent]["ID_TP"][index],
                                                                      "rule_list":[d_3[parent]["bt_list"][index]],
                                                                      "mass": d_3[parent]["mass_TP"][index],
                                                                      "Formula": d_3[parent]["Formula_TP"][index],
                                                                      "source_list" :[method_3_package],
                                                                      "code": d_3[parent]["code_TP"][index],
                                                                      "Structure" : d_3[parent]["Structure_TP"][index],
                                                                      "combined_prob": [],
                                                                      "score": 100, "InchiKey": d_3[parent]["inchi_TP"][index],
                                                                      "alternative_parent" : []}
                        else:
                            if data_dict_com[parent]["TP_dict"].get(tp): # if TP already there then append list
                                if d_3[parent]["bt_list"][index] not in data_dict_com[parent]["TP_dict"][tp]["rule_list"]:
                                    data_dict_com[parent]["TP_dict"][tp]["rule_list"].append(d_3[parent]["bt_list"][index])
                                if method_3_package not in data_dict_com[parent]["TP_dict"][tp]["source_list"]:
                                    data_dict_com[parent]["TP_dict"][tp]["source_list"].append(method_3_package)
                            else:
                                data_dict_com[parent]["TP_dict"][tp] = {"CAS": d_3[parent]["ID_TP"][index],
                                                                                  "rule_list":[d_3[parent]["bt_list"][index]],
                                                                                  "mass": d_3[parent]["mass_TP"][index],
                                                                                  "Formula": d_3[parent]["Formula_TP"][index],
                                                                                  "source_list" :[method_3_package],
                                                                                  "code": d_3[parent]["code_TP"][index],
                                                                                  "Structure" : d_3[parent]["Structure_TP"][index],
                                                                                  "combined_prob": [],
                                                                                  "score": 100, "InchiKey": d_3[parent]["inchi_TP"][index],
                                                                                  "alternative_parent" : []}
                            
    # add TP data from fourth data dict (if specified)                             
    if d_4 != "none":
        for parent in d_4:
            if d_4.get(parent, {}).get("combined_prob") is not None:
                for index, tp in enumerate(d_4[parent]["TP_list_canon"]): 
                    if data_dict_com.get(parent, {}).get("TP_dict") is None:
                        data_dict_com[parent]["TP_dict"] = {}
                        data_dict_com[parent]["TP_dict"][tp] = {"CAS": d_4[parent]["ID_TP"][index],
                                                                      "rule_list":[d_4[parent]["bt_list"][index]],
                                                                      "mass": d_4[parent]["mass_TP"][index],
                                                                      "Formula": d_4[parent]["Formula_TP"][index],
                                                                      "source_list" :[method_4_package],
                                                                      "code": d_4[parent]["code_TP"][index],
                                                                      "Structure" : d_4[parent]["Structure_TP"][index],
                                                                      "combined_prob": [d_4[parent]["combined_prob"][index]],
                                                                      "score": 100, "InchiKey": d_4[parent]["inchi_TP"][index],
                                                                      "alternative_parent" : []}
                    else:
                        if data_dict_com[parent]["TP_dict"].get(tp): # if TP already there then append list
                            data_dict_com[parent]["TP_dict"][tp]["combined_prob"].append(d_4[parent]["combined_prob"][index])
                            if d_4[parent]["bt_list"][index] not in data_dict_com[parent]["TP_dict"][tp]["rule_list"]:
                                data_dict_com[parent]["TP_dict"][tp]["rule_list"].append(d_4[parent]["bt_list"][index])
                            if method_4_package not in data_dict_com[parent]["TP_dict"][tp]["source_list"]:
                                data_dict_com[parent]["TP_dict"][tp]["source_list"].append(method_4_package)
                        else:
                            data_dict_com[parent]["TP_dict"][tp] = {"CAS": d_4[parent]["ID_TP"][index],
                                                                          "rule_list":[d_4[parent]["bt_list"][index]],
                                                                          "mass": d_4[parent]["mass_TP"][index],
                                                                          "Formula": d_4[parent]["Formula_TP"][index],
                                                                          "source_list" :[method_4_package],
                                                                          "code": d_4[parent]["code_TP"][index],
                                                                          "Structure" : d_4[parent]["Structure_TP"][index],
                                                                          "combined_prob": [d_4[parent]["combined_prob"][index]],
                                                                          "score": 100, "InchiKey": d_4[parent]["inchi_TP"][index],
                                                                          "alternative_parent" : []}
    
            else:
                for parent in d_4:
                    for index, tp in enumerate(d_4[parent]["TP_list_canon"]):
                        if data_dict_com.get(parent, {}).get("TP_dict") is None:
                            data_dict_com[parent]["TP_dict"] = {}
                            data_dict_com[parent]["TP_dict"][tp] = {"CAS": d_4[parent]["ID_TP"][index],
                                                                      "rule_list":[d_4[parent]["bt_list"][index]],
                                                                      "mass": d_4[parent]["mass_TP"][index],
                                                                      "Formula": d_4[parent]["Formula_TP"][index],
                                                                      "source_list" :[method_4_package],
                                                                      "code": d_4[parent]["code_TP"][index],
                                                                      "Structure" : d_4[parent]["Structure_TP"][index],
                                                                      "combined_prob": [],
                                                                      "score": 100, "InchiKey": d_4[parent]["inchi_TP"][index],
                                                                      "alternative_parent" : []}
                        else:
                            if data_dict_com[parent]["TP_dict"].get(tp): # if TP already there then append list
                                if d_4[parent]["bt_list"][index] not in data_dict_com[parent]["TP_dict"][tp]["rule_list"]:
                                    data_dict_com[parent]["TP_dict"][tp]["rule_list"].append(d_4[parent]["bt_list"][index])
                                if method_4_package not in data_dict_com[parent]["TP_dict"][tp]["source_list"]:
                                    data_dict_com[parent]["TP_dict"][tp]["source_list"].append(method_4_package)
                            else:
                                data_dict_com[parent]["TP_dict"][tp] = {"CAS": d_4[parent]["ID_TP"][index],
                                                                                  "rule_list":[d_4[parent]["bt_list"][index]],
                                                                                  "mass": d_4[parent]["mass_TP"][index],
                                                                                  "Formula": d_4[parent]["Formula_TP"][index],
                                                                                  "source_list" :[method_4_package],
                                                                                  "code": d_4[parent]["code_TP"][index],
                                                                                  "Structure" : d_4[parent]["Structure_TP"][index],
                                                                                  "combined_prob": [],
                                                                                  "score": 100, "InchiKey": d_4[parent]["inchi_TP"][index],
                                                                                  "alternative_parent" : []}               

    # add TP data from fifth data dict (if specified)                             
    if d_5 != "none":
        for parent in d_5:
            if d_5.get(parent, {}).get("combined_prob") is not None:
                for index, tp in enumerate(d_5[parent]["TP_list_canon"]): 
                    if data_dict_com.get(parent, {}).get("TP_dict") is None:
                        data_dict_com[parent]["TP_dict"] = {}
                        data_dict_com[parent]["TP_dict"][tp] = {"CAS": d_5[parent]["ID_TP"][index],
                                                                      "rule_list":[d_5[parent]["bt_list"][index]],
                                                                      "mass": d_5[parent]["mass_TP"][index],
                                                                      "Formula": d_5[parent]["Formula_TP"][index],
                                                                      "source_list" :[method_5_package],
                                                                      "code": d_5[parent]["code_TP"][index],
                                                                      "Structure" : d_5[parent]["Structure_TP"][index],
                                                                      "combined_prob": [d_5[parent]["combined_prob"][index]],
                                                                      "score": 100, "InchiKey": d_5[parent]["inchi_TP"][index],
                                                                      "alternative_parent" : []}
                    else:
                        if data_dict_com[parent]["TP_dict"].get(tp): # if TP already there then append list
                            data_dict_com[parent]["TP_dict"][tp]["combined_prob"].append(d_5[parent]["combined_prob"][index])
                            if d_5[parent]["bt_list"][index] not in data_dict_com[parent]["TP_dict"][tp]["rule_list"]:
                                data_dict_com[parent]["TP_dict"][tp]["rule_list"].append(d_5[parent]["bt_list"][index])
                            if method_5_package not in data_dict_com[parent]["TP_dict"][tp]["source_list"]:
                                data_dict_com[parent]["TP_dict"][tp]["source_list"].append(method_5_package)
                        else:
                            data_dict_com[parent]["TP_dict"][tp] = {"CAS": d_5[parent]["ID_TP"][index],
                                                                          "rule_list":[d_5[parent]["bt_list"][index]],
                                                                          "mass": d_5[parent]["mass_TP"][index],
                                                                          "Formula": d_5[parent]["Formula_TP"][index],
                                                                          "source_list" :[method_5_package],
                                                                          "code": d_5[parent]["code_TP"][index],
                                                                          "Structure" : d_5[parent]["Structure_TP"][index],
                                                                          "combined_prob": [d_5[parent]["combined_prob"][index]],
                                                                          "score": 100, "InchiKey": d_5[parent]["inchi_TP"][index],
                                                                          "alternative_parent" : []}
    
            else:
                for parent in d_5:
                    for index, tp in enumerate(d_5[parent]["TP_list_canon"]):
                        if data_dict_com.get(parent, {}).get("TP_dict") is None:
                            data_dict_com[parent]["TP_dict"] = {}
                            data_dict_com[parent]["TP_dict"][tp] = {"CAS": d_5[parent]["ID_TP"][index],
                                                                      "rule_list":[d_5[parent]["bt_list"][index]],
                                                                      "mass": d_5[parent]["mass_TP"][index],
                                                                      "Formula": d_5[parent]["Formula_TP"][index],
                                                                      "source_list" :[method_5_package],
                                                                      "code": d_5[parent]["code_TP"][index],
                                                                      "Structure" : d_5[parent]["Structure_TP"][index],
                                                                      "combined_prob": [],
                                                                      "score": 100, "InchiKey": d_5[parent]["inchi_TP"][index],
                                                                      "alternative_parent" : []}
                        else:
                            if data_dict_com[parent]["TP_dict"].get(tp): # if TP already there then append list
                                if d_5[parent]["bt_list"][index] not in data_dict_com[parent]["TP_dict"][tp]["rule_list"]:
                                    data_dict_com[parent]["TP_dict"][tp]["rule_list"].append(d_5[parent]["bt_list"][index])
                                if method_5_package not in data_dict_com[parent]["TP_dict"][tp]["source_list"]:
                                    data_dict_com[parent]["TP_dict"][tp]["source_list"].append(method_5_package)
                            else:
                                data_dict_com[parent]["TP_dict"][tp] = {"CAS": d_5[parent]["ID_TP"][index],
                                                                                  "rule_list":[d_5[parent]["bt_list"][index]],
                                                                                  "mass": d_5[parent]["mass_TP"][index],
                                                                                  "Formula": d_5[parent]["Formula_TP"][index],
                                                                                  "source_list" :[method_5_package],
                                                                                  "code": d_5[parent]["code_TP"][index],
                                                                                  "Structure" : d_5[parent]["Structure_TP"][index],
                                                                                  "combined_prob": [],
                                                                                  "score": 100, "InchiKey": d_5[parent]["inchi_TP"][index],
                                                                                  "alternative_parent" : []}               

    # add TP data from sixth data dict (if specified)                             
    if d_6 != "none":
        for parent in d_6:
            if d_6.get(parent, {}).get("combined_prob") is not None:
                for index, tp in enumerate(d_6[parent]["TP_list_canon"]): 
                    if data_dict_com.get(parent, {}).get("TP_dict") is None:
                        data_dict_com[parent]["TP_dict"] = {}
                        data_dict_com[parent]["TP_dict"][tp] = {"CAS": d_6[parent]["ID_TP"][index],
                                                                      "rule_list":[d_6[parent]["bt_list"][index]],
                                                                      "mass": d_6[parent]["mass_TP"][index],
                                                                      "Formula": d_6[parent]["Formula_TP"][index],
                                                                      "source_list" :[method_6_package],
                                                                      "code": d_6[parent]["code_TP"][index],
                                                                      "Structure" : d_6[parent]["Structure_TP"][index],
                                                                      "combined_prob": [d_6[parent]["combined_prob"][index]],
                                                                      "score": 100, "InchiKey": d_6[parent]["inchi_TP"][index],
                                                                      "alternative_parent" : []}
                    else:
                        if data_dict_com[parent]["TP_dict"].get(tp): # if TP already there then append list
                            data_dict_com[parent]["TP_dict"][tp]["combined_prob"].append(d_6[parent]["combined_prob"][index])
                            if d_6[parent]["bt_list"][index] not in data_dict_com[parent]["TP_dict"][tp]["rule_list"]:
                                data_dict_com[parent]["TP_dict"][tp]["rule_list"].append(d_6[parent]["bt_list"][index])
                            if method_6_package not in data_dict_com[parent]["TP_dict"][tp]["source_list"]:
                                data_dict_com[parent]["TP_dict"][tp]["source_list"].append(method_6_package)
                        else:
                            data_dict_com[parent]["TP_dict"][tp] = {"CAS": d_6[parent]["ID_TP"][index],
                                                                          "rule_list":[d_6[parent]["bt_list"][index]],
                                                                          "mass": d_6[parent]["mass_TP"][index],
                                                                          "Formula": d_6[parent]["Formula_TP"][index],
                                                                          "source_list" :[method_6_package],
                                                                          "code": d_6[parent]["code_TP"][index],
                                                                          "Structure" : d_6[parent]["Structure_TP"][index],
                                                                          "combined_prob": [d_6[parent]["combined_prob"][index]],
                                                                          "score": 100, "InchiKey": d_6[parent]["inchi_TP"][index],
                                                                          "alternative_parent" : []}
    
            else:
                for parent in d_6:
                    for index, tp in enumerate(d_6[parent]["TP_list_canon"]):
                        if data_dict_com.get(parent, {}).get("TP_dict") is None:
                            data_dict_com[parent]["TP_dict"] = {}
                            data_dict_com[parent]["TP_dict"][tp] = {"CAS": d_6[parent]["ID_TP"][index],
                                                                      "rule_list":[d_6[parent]["bt_list"][index]],
                                                                      "mass": d_6[parent]["mass_TP"][index],
                                                                      "Formula": d_6[parent]["Formula_TP"][index],
                                                                      "source_list" :[method_6_package],
                                                                      "code": d_6[parent]["code_TP"][index],
                                                                      "Structure" : d_6[parent]["Structure_TP"][index],
                                                                      "combined_prob": [],
                                                                      "score": 100, "InchiKey": d_6[parent]["inchi_TP"][index],
                                                                      "alternative_parent" : []}
                        else:
                            if data_dict_com[parent]["TP_dict"].get(tp): # if TP already there then append list
                                if d_6[parent]["bt_list"][index] not in data_dict_com[parent]["TP_dict"][tp]["rule_list"]:
                                    data_dict_com[parent]["TP_dict"][tp]["rule_list"].append(d_6[parent]["bt_list"][index])
                                if method_6_package not in data_dict_com[parent]["TP_dict"][tp]["source_list"]:
                                    data_dict_com[parent]["TP_dict"][tp]["source_list"].append(method_6_package)
                            else:
                                data_dict_com[parent]["TP_dict"][tp] = {"CAS": d_6[parent]["ID_TP"][index],
                                                                                  "rule_list":[d_6[parent]["bt_list"][index]],
                                                                                  "mass": d_6[parent]["mass_TP"][index],
                                                                                  "Formula": d_6[parent]["Formula_TP"][index],
                                                                                  "source_list" :[method_6_package],
                                                                                  "code": d_6[parent]["code_TP"][index],
                                                                                  "Structure" : d_6[parent]["Structure_TP"][index],
                                                                                  "combined_prob": [],
                                                                                  "score": 100, "InchiKey": d_6[parent]["inchi_TP"][index],
                                                                                  "alternative_parent" : []}               


    # add TP data from seventh data dict (if specified)                             
    if d_7 != "none":
        for parent in d_7:
            if d_7.get(parent, {}).get("combined_prob") is not None:
                for index, tp in enumerate(d_7[parent]["TP_list_canon"]): 
                    if data_dict_com.get(parent, {}).get("TP_dict") is None:
                        data_dict_com[parent]["TP_dict"] = {}
                        data_dict_com[parent]["TP_dict"][tp] = {"CAS": d_7[parent]["ID_TP"][index],
                                                                      "rule_list":[d_7[parent]["bt_list"][index]],
                                                                      "mass": d_7[parent]["mass_TP"][index],
                                                                      "Formula": d_7[parent]["Formula_TP"][index],
                                                                      "source_list" :[method_7_package],
                                                                      "code": d_7[parent]["code_TP"][index],
                                                                      "Structure" : d_7[parent]["Structure_TP"][index],
                                                                      "combined_prob": [d_7[parent]["combined_prob"][index]],
                                                                      "score": 100, "InchiKey": d_7[parent]["inchi_TP"][index],
                                                                      "alternative_parent" : []}
                    else:
                        if data_dict_com[parent]["TP_dict"].get(tp): # if TP already there then append list
                            data_dict_com[parent]["TP_dict"][tp]["combined_prob"].append(d_7[parent]["combined_prob"][index])
                            if d_7[parent]["bt_list"][index] not in data_dict_com[parent]["TP_dict"][tp]["rule_list"]:
                                data_dict_com[parent]["TP_dict"][tp]["rule_list"].append(d_7[parent]["bt_list"][index])
                            if method_7_package not in data_dict_com[parent]["TP_dict"][tp]["source_list"]:
                                data_dict_com[parent]["TP_dict"][tp]["source_list"].append(method_7_package)
                        else:
                            data_dict_com[parent]["TP_dict"][tp] = {"CAS": d_7[parent]["ID_TP"][index],
                                                                          "rule_list":[d_7[parent]["bt_list"][index]],
                                                                          "mass": d_7[parent]["mass_TP"][index],
                                                                          "Formula": d_7[parent]["Formula_TP"][index],
                                                                          "source_list" :[method_7_package],
                                                                          "code": d_7[parent]["code_TP"][index],
                                                                          "Structure" : d_7[parent]["Structure_TP"][index],
                                                                          "combined_prob": [d_7[parent]["combined_prob"][index]],
                                                                          "score": 100, "InchiKey": d_7[parent]["inchi_TP"][index],
                                                                          "alternative_parent" : []}
    
            else:
                for parent in d_7:
                    for index, tp in enumerate(d_7[parent]["TP_list_canon"]):
                        if data_dict_com.get(parent, {}).get("TP_dict") is None:
                            data_dict_com[parent]["TP_dict"] = {}
                            data_dict_com[parent]["TP_dict"][tp] = {"CAS": d_7[parent]["ID_TP"][index],
                                                                      "rule_list":[d_7[parent]["bt_list"][index]],
                                                                      "mass": d_7[parent]["mass_TP"][index],
                                                                      "Formula": d_7[parent]["Formula_TP"][index],
                                                                      "source_list" :[method_7_package],
                                                                      "code": d_7[parent]["code_TP"][index],
                                                                      "Structure" : d_7[parent]["Structure_TP"][index],
                                                                      "combined_prob": [],
                                                                      "score": 100, "InchiKey": d_7[parent]["inchi_TP"][index],
                                                                      "alternative_parent" : []}
                        else:
                            if data_dict_com[parent]["TP_dict"].get(tp): # if TP already there then append list
                                if d_7[parent]["bt_list"][index] not in data_dict_com[parent]["TP_dict"][tp]["rule_list"]:
                                    data_dict_com[parent]["TP_dict"][tp]["rule_list"].append(d_7[parent]["bt_list"][index])
                                if method_7_package not in data_dict_com[parent]["TP_dict"][tp]["source_list"]:
                                    data_dict_com[parent]["TP_dict"][tp]["source_list"].append(method_7_package)
                            else:
                                data_dict_com[parent]["TP_dict"][tp] = {"CAS": d_7[parent]["ID_TP"][index],
                                                                                  "rule_list":[d_7[parent]["bt_list"][index]],
                                                                                  "mass": d_7[parent]["mass_TP"][index],
                                                                                  "Formula": d_7[parent]["Formula_TP"][index],
                                                                                  "source_list" :[method_7_package],
                                                                                  "code": d_7[parent]["code_TP"][index],
                                                                                  "Structure" : d_7[parent]["Structure_TP"][index],
                                                                                  "combined_prob": [],
                                                                                  "score": 100, "InchiKey": d_7[parent]["inchi_TP"][index],
                                                                                  "alternative_parent" : []} 

    # add TP data from eighth data dict (if specified)                             
    if d_8 != "none":
        for parent in d_8:
            if d_8.get(parent, {}).get("combined_prob") is not None:
                for index, tp in enumerate(d_8[parent]["TP_list_canon"]): 
                    if data_dict_com.get(parent, {}).get("TP_dict") is None:
                        data_dict_com[parent]["TP_dict"] = {}
                        data_dict_com[parent]["TP_dict"][tp] = {"CAS": d_8[parent]["ID_TP"][index],
                                                                      "rule_list":[d_8[parent]["bt_list"][index]],
                                                                      "mass": d_8[parent]["mass_TP"][index],
                                                                      "Formula": d_8[parent]["Formula_TP"][index],
                                                                      "source_list" :[method_8_package],
                                                                      "code": d_8[parent]["code_TP"][index],
                                                                      "Structure" : d_8[parent]["Structure_TP"][index],
                                                                      "combined_prob": [d_8[parent]["combined_prob"][index]],
                                                                      "score": 100, "InchiKey": d_8[parent]["inchi_TP"][index],
                                                                      "alternative_parent" : []}
                    else:
                        if data_dict_com[parent]["TP_dict"].get(tp): # if TP already there then append list
                            data_dict_com[parent]["TP_dict"][tp]["combined_prob"].append(d_8[parent]["combined_prob"][index])
                            if d_8[parent]["bt_list"][index] not in data_dict_com[parent]["TP_dict"][tp]["rule_list"]:
                                data_dict_com[parent]["TP_dict"][tp]["rule_list"].append(d_8[parent]["bt_list"][index])
                            if method_8_package not in data_dict_com[parent]["TP_dict"][tp]["source_list"]:
                                data_dict_com[parent]["TP_dict"][tp]["source_list"].append(method_8_package)
                        else:
                            data_dict_com[parent]["TP_dict"][tp] = {"CAS": d_8[parent]["ID_TP"][index],
                                                                          "rule_list":[d_8[parent]["bt_list"][index]],
                                                                          "mass": d_8[parent]["mass_TP"][index],
                                                                          "Formula": d_8[parent]["Formula_TP"][index],
                                                                          "source_list" :[method_8_package],
                                                                          "code": d_8[parent]["code_TP"][index],
                                                                          "Structure" : d_8[parent]["Structure_TP"][index],
                                                                          "combined_prob": [d_8[parent]["combined_prob"][index]],
                                                                          "score": 100, "InchiKey": d_8[parent]["inchi_TP"][index],
                                                                          "alternative_parent" : []}
    
            else:
                for parent in d_8:
                    for index, tp in enumerate(d_8[parent]["TP_list_canon"]):
                        if data_dict_com.get(parent, {}).get("TP_dict") is None:
                            data_dict_com[parent]["TP_dict"] = {}
                            data_dict_com[parent]["TP_dict"][tp] = {"CAS": d_8[parent]["ID_TP"][index],
                                                                      "rule_list":[d_8[parent]["bt_list"][index]],
                                                                      "mass": d_8[parent]["mass_TP"][index],
                                                                      "Formula": d_8[parent]["Formula_TP"][index],
                                                                      "source_list" :[method_8_package],
                                                                      "code": d_8[parent]["code_TP"][index],
                                                                      "Structure" : d_8[parent]["Structure_TP"][index],
                                                                      "combined_prob": [],
                                                                      "score": 100, "InchiKey": d_8[parent]["inchi_TP"][index],
                                                                      "alternative_parent" : []}
                        else:
                            if data_dict_com[parent]["TP_dict"].get(tp): # if TP already there then append list
                                if d_8[parent]["bt_list"][index] not in data_dict_com[parent]["TP_dict"][tp]["rule_list"]:
                                    data_dict_com[parent]["TP_dict"][tp]["rule_list"].append(d_8[parent]["bt_list"][index])
                                if method_8_package not in data_dict_com[parent]["TP_dict"][tp]["source_list"]:
                                    data_dict_com[parent]["TP_dict"][tp]["source_list"].append(method_8_package)
                            else:
                                data_dict_com[parent]["TP_dict"][tp] = {"CAS": d_8[parent]["ID_TP"][index],
                                                                                  "rule_list":[d_8[parent]["bt_list"][index]],
                                                                                  "mass": d_8[parent]["mass_TP"][index],
                                                                                  "Formula": d_8[parent]["Formula_TP"][index],
                                                                                  "source_list" :[method_8_package],
                                                                                  "code": d_8[parent]["code_TP"][index],
                                                                                  "Structure" : d_8[parent]["Structure_TP"][index],
                                                                                  "combined_prob": [],
                                                                                  "score": 100, "InchiKey": d_8[parent]["inchi_TP"][index],
                                                                                  "alternative_parent" : []} 
                                                   
    do_pickle(data_dict_com, "data_dict_com.pickle")                                               
    return data_dict_com

def score_dict(data_dict_com):
    # Scoring system, each TP starts with a score of 100 made up points
    print("Scoring dictionary")
    # check the mass
    for parent in data_dict_com:
        for tp in data_dict_com[parent]["TP_dict"]:
            if data_dict_com[parent]["TP_dict"][tp]["mass"] < 100: # if the mass of the TPs is below 100 u then set the score to 0
                data_dict_com[parent]["TP_dict"][tp]["score"] = data_dict_com[parent]["TP_dict"][tp]["score"] - 100   
    # check if TP has CAS number
    # if TP has CAS number then it was studied before and we want to look for new TPs, however if it has no CAS number then the chance is higher that it is only a wanky prediction and not an actual TP that is observed in the environment
    for parent in data_dict_com:
        for tp in data_dict_com[parent]["TP_dict"]:
            if len(data_dict_com[parent]["TP_dict"][tp]["CAS"]) > 1: # if TP has (at least one) CAS number, then reduce score (empty string "" has len 0)
                data_dict_com[parent]["TP_dict"][tp]["score"] = data_dict_com[parent]["TP_dict"][tp]["score"] - 10   
    # check if TP was predicted by other method
    for parent in data_dict_com:
        for tp in data_dict_com[parent]["TP_dict"]:
            if len(data_dict_com[parent]["TP_dict"][tp]["source_list"]) == 1: # if TP was predicted by only one method
                data_dict_com[parent]["TP_dict"][tp]["score"] = data_dict_com[parent]["TP_dict"][tp]["score"] - 50
            if len(data_dict_com[parent]["TP_dict"][tp]["source_list"]) == 2: # if TP was predicted by two methods
                data_dict_com[parent]["TP_dict"][tp]["score"] = data_dict_com[parent]["TP_dict"][tp]["score"] - 40                
            if len(data_dict_com[parent]["TP_dict"][tp]["source_list"]) == 3: 
                data_dict_com[parent]["TP_dict"][tp]["score"] = data_dict_com[parent]["TP_dict"][tp]["score"] - 30
            if len(data_dict_com[parent]["TP_dict"][tp]["source_list"]) == 4: 
                data_dict_com[parent]["TP_dict"][tp]["score"] = data_dict_com[parent]["TP_dict"][tp]["score"] - 20    
            if len(data_dict_com[parent]["TP_dict"][tp]["source_list"]) == 5: 
                data_dict_com[parent]["TP_dict"][tp]["score"] = data_dict_com[parent]["TP_dict"][tp]["score"] - 15
            if len(data_dict_com[parent]["TP_dict"][tp]["source_list"]) == 6: 
                data_dict_com[parent]["TP_dict"][tp]["score"] = data_dict_com[parent]["TP_dict"][tp]["score"] - 10
            if len(data_dict_com[parent]["TP_dict"][tp]["source_list"]) == 7: 
                data_dict_com[parent]["TP_dict"][tp]["score"] = data_dict_com[parent]["TP_dict"][tp]["score"] - 5
                # if TP is predicted by all 8 methods then it doesn't change the score    
    # check the combined probability (is saved as str so must convert to float)(only for enviPath files possible)
    for parent in data_dict_com:
        for tp in data_dict_com[parent]["TP_dict"]: # check first if it comes from enviPath method because otherwise it doesnt have a probability
            if data_dict_com.get(parent, {}).get(tp, {}).get("combined_prob") is not None:
                # each time a condition is fullfilled, the score is reduced (careful, the penalty is summed up!)
                if float(data_dict_com[parent]["TP_dict"][tp]["combined_prob"]) < 0.3:
                    data_dict_com[parent]["TP_dict"][tp]["score"] = data_dict_com[parent]["TP_dict"][tp]["score"] - 5
                if float(data_dict_com[parent]["TP_dict"][tp]["combined_prob"]) < 0.1:
                    data_dict_com[parent]["TP_dict"][tp]["score"] = data_dict_com[parent]["TP_dict"][tp]["score"] - 5
                if float(data_dict_com[parent]["TP_dict"][tp]["combined_prob"]) < 0.01:
                    data_dict_com[parent]["TP_dict"][tp]["score"] = data_dict_com[parent]["TP_dict"][tp]["score"] - 5
                if float(data_dict_com[parent]["TP_dict"][tp]["combined_prob"]) < 0.005:
                    data_dict_com[parent]["TP_dict"][tp]["score"] = data_dict_com[parent]["TP_dict"][tp]["score"] - 10
                if float(data_dict_com[parent]["TP_dict"][tp]["combined_prob"]) < 0.001:
                    data_dict_com[parent]["TP_dict"][tp]["score"] = data_dict_com[parent]["TP_dict"][tp]["score"] - 10
    # now the unwanted TPs will be deleted (delete the entries with very low score)
    # first copy the data dict com, so we can iterate over the copy, otherwise it gives runtime error: dictionary changed size during iteration
    data_dict_com_copy = copy.deepcopy(data_dict_com)  
    for parent in data_dict_com_copy:
        for tp in data_dict_com_copy[parent]["TP_dict"]:
            if data_dict_com_copy[parent]["TP_dict"][tp]["score"] <= 0:
                del data_dict_com[parent]["TP_dict"][tp]                 
    # check for dublicates, if TP is found for other parent then append alternative_parent list of that parent
    for parent in data_dict_com:
        for tp in data_dict_com[parent]["TP_dict"]:
            for parent2 in data_dict_com:
                if tp in data_dict_com[parent2]["TP_dict"]:
                    if data_dict_com[parent]["code_parent"][0] not in data_dict_com[parent2]["TP_dict"][tp]["code"]:
                        # data_dict_com[parent2]["TP_dict"][tp]["alternative_parent"].append(data_dict_com[parent]["code_parent"])
                        data_dict_com[parent2]["TP_dict"][tp]["alternative_parent"] = data_dict_com[parent]["code_parent"]
                else:
                    continue

    # create new dict to which the removed TPs are added
    removed_tps_dict = {}             
    # if a parent has more than max allowed TPs, they get removed starting with lowest score
    # first copy the data dict com, so we can iterate over the copy, otherwise it gives runtime error: dictionary changed size during iteration
    data_dict_com_copy_3 = copy.deepcopy(data_dict_com)   
    for parent in data_dict_com_copy_3:
        temp_tp_list = [] # create temporary lists containing the TP smiles and corresponding score of one parent
        temp_score_list = []      
        # create entries for every parent
        removed_tps_dict[parent] = {'TP_list': [], "parent_name" : [], "tp_name" : []}  
        for tp in data_dict_com_copy_3[parent]["TP_dict"]:
            temp_tp_list.append(tp)
            temp_score_list.append(data_dict_com_copy_3[parent]["TP_dict"][tp]["score"])
            while len(temp_tp_list) > max_TP_per_parent: # if parent has more than max_TP_per_parent TPs proceed
                index = temp_score_list.index(min(temp_score_list)) # find the first index of the lowest score 
                # add smiles of TP that is deleted
                removed_tps_dict[parent]['TP_list'].append(temp_tp_list[index])
                removed_tps_dict[parent]["parent_name"].append(data_dict_com_copy_3[parent]["code_parent"])
                removed_tps_dict[parent]["tp_name"].append(data_dict_com_copy_3[parent]["TP_dict"][tp]["code"])
                del data_dict_com[parent]["TP_dict"][(temp_tp_list[index])] # index of score is the same as the corresponding TP
                temp_score_list.pop(index) # now remove the score and the TP smiles from temporary list
                temp_tp_list.pop(index)
    do_pickle(data_dict_com, "data_dict_com_scored.pickle")
    parent_removed_list = []
    tp_removed_list = []
    tp_code_removed_list = []
    parent_code_removed_list = []
    for parent in removed_tps_dict:
        for tp in removed_tps_dict[parent]['TP_list']:
            parent_removed_list.append(parent)
            tp_removed_list.append(tp)
        for name in removed_tps_dict[parent]['parent_name']:
            parent_code_removed_list.append(name)
        for code in removed_tps_dict[parent]['tp_name']:
            tp_code_removed_list.append(code)

    df_removed_dict = {"Parent Code": parent_code_removed_list  ,"Parent SMILES": parent_removed_list, "TP Code" : tp_code_removed_list, "TP SMILES": tp_removed_list}
    df_removed = pd.DataFrame.from_dict(df_removed_dict)
    df_removed.to_csv(output_removed_tps , index = False, sep = "\t")    
    
    return data_dict_com

################################################################################################################################################################################################################################################################

# START SCRIPT

assert type(file_location_1) == str, "file_location_1 must be a string!"
assert type(prediction_method_1) == str, "prediction_method_1 must be a string!"
assert type(package_method_1) == str, "package_method_1 must be a string!"
assert type(file_location_2) == str, "file_location_2 must be a string!"
assert type(prediction_method_2) == str, "prediction_method_2 must be a string!"
assert type(package_method_2) == str, "package_method_2 must be a string!"
assert type(file_location_3) == str, "file_location_3 must be a string!"
assert type(prediction_method_3) == str, "prediction_method_3 must be a string!"
assert type(package_method_3) == str, "package_method_3 must be a string!"
assert type(file_location_4) == str, "file_location_4 must be a string!"
assert type(prediction_method_4) == str, "prediction_method_4 must be a string!"
assert type(package_method_4) == str, "package_method_4 must be a string!"
assert type(file_location_5) == str, "file_location_5 must be a string!"
assert type(prediction_method_5) == str, "prediction_method_5 must be a string!"
assert type(package_method_5) == str, "package_method_5 must be a string!"
assert type(file_location_6) == str, "file_location_6 must be a string!"
assert type(prediction_method_6) == str, "prediction_method_6 must be a string!"
assert type(package_method_6) == str, "package_method_6 must be a string!"
assert type(file_location_7) == str, "file_location_7 must be a string!"
assert type(prediction_method_7) == str, "prediction_method_7 must be a string!"
assert type(package_method_7) == str, "package_method_7 must be a string!"
assert type(file_location_8) == str, "file_location_8 must be a string!"
assert type(prediction_method_8) == str, "prediction_method_8 must be a string!"
assert type(package_method_8) == str, "package_method_8 must be a string!"
assert type(code_location) == str, "code_location must be a string!"
assert type(smi_location) == str, "smi_location must be a string!"
assert type(name_location) == str, "name_location must be a string!"
assert type(output_location_1) == str, "output_location_1 must be a string!"
assert type(output_location_2) == str, "output_location_2 must be a string!"
assert type(output_location_3) == str, "output_location_3 must be a string!"
assert type(output_location_4) == str, "output_location_4 must be a string!"
assert type(output_location_5) == str, "output_location_5 must be a string!"
assert type(output_location_6) == str, "output_location_6 must be a string!"
assert type(output_file_CD_masslist) == str, "output_file_CD_masslist must be a string!"
assert type(output_file_all_data) == str, "output_file_all_data must be a string!"
assert type(output_inclusion_pos) == str, "output_inclusion_pos must be a string!"
assert type(output_inclusion_neg) == str, "output_inclusion_neg must be a string!"
assert type(scoring_system) == bool, "scoring_sytem must be either 'True' or 'False'!"
assert type(max_TP_per_parent) == int, "max_TP_per_parent must be an integer (e.g. 15 or 50)!"

def load_mapping_files(code_location, smi_location, name_location):
    # read files with code, name and SMILES of selected parents to get dictionaries, the order is important!
    code_file = open(code_location)
    code_list = []
    for line2 in code_file:
        code_list.append(line2.rstrip())
    SMILES_comp_file = open(smi_location)
    SMILES_list = []
    for line3 in SMILES_comp_file:
        SMILES_list.append(line3.rstrip())
    name_file = open(name_location)
    name_list = []
    for line4 in name_file:
        name_list.append(line4.rstrip())
    code_dict = dict(zip(SMILES_list, code_list))  # dictionary with code of selected compounds with corresponding SMILES
    name_dict = dict(zip(SMILES_list, name_list))  # dictionary with name of selected compounds with corresponding SMILES
    return code_dict, name_dict

################################################################################################################################################################################################################################################################
#    MAIN
################################################################################################################################################################################################################################################################
code_dict, name_dict = load_mapping_files(code_location, smi_location, name_location)

# read first file
data_dict_1 = file_to_csv(file_location_1, "data_dict_" + package_method_1 + ".pickle", output_location_1, prediction_method_1) 
# read second file (if available)
if consider_file_2 == "yes":
    data_dict_2 = file_to_csv(file_location_2, "data_dict_" + package_method_2 + ".pickle", output_location_2, prediction_method_2)
# read third file (if available)
if consider_file_3 == "yes":
    data_dict_3 = file_to_csv(file_location_3, "data_dict_" + package_method_3 + ".pickle", output_location_3, prediction_method_3)
# read fourth file (if available)
if consider_file_4 == "yes":
    data_dict_4 = file_to_csv(file_location_4, "data_dict_" + package_method_4 + ".pickle", output_location_4, prediction_method_4)  
# read fifth file (if available)
if consider_file_5 == "yes":
    data_dict_5 = file_to_csv(file_location_5, "data_dict_" + package_method_5 + ".pickle", output_location_5, prediction_method_5)  
# read sixth file (if available)
if consider_file_6 == "yes":
    data_dict_6 = file_to_csv(file_location_6, "data_dict_" + package_method_6 + ".pickle", output_location_6, prediction_method_6)     
# read seventh file (if available)
if consider_file_7 == "yes":
    data_dict_7 = file_to_csv(file_location_7, "data_dict_" + package_method_7 + ".pickle", output_location_7, prediction_method_7)    
# read eighth file (if available)
if consider_file_8 == "yes":
    data_dict_8 = file_to_csv(file_location_8, "data_dict_" + package_method_8 + ".pickle", output_location_8, prediction_method_8)
        
# import saved dictionaries  
# data_dict_1 = get_pickle("data_dict_1.pickle")
# data_dict_2 = get_pickle("data_dict_2.pickle")
# data_dict_3 = get_pickle("data_dict_3.pickle")
# data_dict_4 = get_pickle("data_dict_4.pickle")
# data_dict_5 = get_pickle("data_dict_5.pickle")
# data_dict_6 = get_pickle("data_dict_6.pickle")
# data_dict_7 = get_pickle("data_dict_7.pickle")
# data_dict_8 = get_pickle("data_dict_8.pickle")

# combine all data dicts into one and apply scoring system
# also, generate mass list for CD, csv file with all the data and iclusion lists for QExactivePlus
if scoring_system == True:
    if consider_file_2 == "yes":
        if consider_file_3 == "yes":
            if consider_file_4 == "yes":
                if consider_file_5 == "yes":
                    if consider_file_6 == "yes":
                        if consider_file_7 == "yes":
                            if consider_file_8 == "yes":
                                data_dict_com = combine_dict(data_dict_1, package_method_1, data_dict_2, package_method_2, data_dict_3, package_method_3, data_dict_4, package_method_4, data_dict_5, package_method_5, data_dict_6, package_method_6, data_dict_7, package_method_7, data_dict_8, package_method_8)
                                data_dict_com_scored = score_dict(data_dict_com)
                                combined_dict_to_csv(data_dict_com_scored, output_file_CD_masslist, output_file_all_data)
                            else:
                                data_dict_com = combine_dict(data_dict_1, package_method_1, data_dict_2, package_method_2, data_dict_3, package_method_3, data_dict_4, package_method_4, data_dict_5, package_method_5, data_dict_6, package_method_6, data_dict_7, package_method_7, "none", "none")
                                data_dict_com_scored = score_dict(data_dict_com)
                                combined_dict_to_csv(data_dict_com_scored, output_file_CD_masslist, output_file_all_data)                                
                        else:
                            data_dict_com = combine_dict(data_dict_1, package_method_1, data_dict_2, package_method_2, data_dict_3, package_method_3, data_dict_4, package_method_4, data_dict_5, package_method_5, data_dict_6, package_method_6, "none", "none", "none", "none")
                            data_dict_com_scored = score_dict(data_dict_com)
                            combined_dict_to_csv(data_dict_com_scored, output_file_CD_masslist, output_file_all_data)
                    else:
                        data_dict_com = combine_dict(data_dict_1, package_method_1, data_dict_2, package_method_2, data_dict_3, package_method_3, data_dict_4, package_method_4, data_dict_5, package_method_5, "none", "none", "none", "none", "none", "none")
                        data_dict_com_scored = score_dict(data_dict_com)
                        combined_dict_to_csv(data_dict_com_scored, output_file_CD_masslist, output_file_all_data)
                else:
                    data_dict_com = combine_dict(data_dict_1, package_method_1, data_dict_2, package_method_2, data_dict_3, package_method_3, data_dict_4, package_method_4, "none", "none", "none", "none", "none", "none", "none", "none")
                    data_dict_com_scored = score_dict(data_dict_com)
                    combined_dict_to_csv(data_dict_com_scored, output_file_CD_masslist, output_file_all_data)
            else:
                data_dict_com = combine_dict(data_dict_1, package_method_1, data_dict_2, package_method_2, data_dict_3, package_method_3, "none", "none", "none", "none", "none", "none", "none", "none", "none", "none")
                data_dict_com_scored = score_dict(data_dict_com)
                combined_dict_to_csv(data_dict_com_scored, output_file_CD_masslist, output_file_all_data)
        else:
            data_dict_com = combine_dict(data_dict_1, package_method_1, data_dict_2, package_method_2, "none", "none", "none", "none", "none", "none", "none", "none", "none", "none", "none", "none")
            data_dict_com_scored = score_dict(data_dict_com)
            combined_dict_to_csv(data_dict_com_scored, output_file_CD_masslist, output_file_all_data)

if scoring_system == False:
    if consider_file_2 == "yes":
        if consider_file_3 == "yes":
            if consider_file_4 == "yes":
                if consider_file_5 == "yes":
                    if consider_file_6 == "yes":
                        if consider_file_7 == "yes":
                            if consider_file_8 == "yes":
                                data_dict_com = combine_dict(data_dict_1, package_method_1, data_dict_2, package_method_2, data_dict_3, package_method_3, data_dict_4, package_method_4, data_dict_5, package_method_5, data_dict_6, package_method_6, data_dict_7, package_method_7, data_dict_8, package_method_8)
                                combined_dict_to_csv(data_dict_com, output_file_CD_masslist, output_file_all_data)
                            else:
                                data_dict_com = combine_dict(data_dict_1, package_method_1, data_dict_2, package_method_2, data_dict_3, package_method_3, data_dict_4, package_method_4, data_dict_5, package_method_5, data_dict_6, package_method_6, data_dict_7, package_method_7, "none", "none")
                                combined_dict_to_csv(data_dict_com, output_file_CD_masslist, output_file_all_data)                                
                        else:
                            data_dict_com = combine_dict(data_dict_1, package_method_1, data_dict_2, package_method_2, data_dict_3, package_method_3, data_dict_4, package_method_4, data_dict_5, package_method_5, data_dict_6, package_method_6, "none", "none", "none", "none")
                            combined_dict_to_csv(data_dict_com, output_file_CD_masslist, output_file_all_data)
                    else:
                        data_dict_com = combine_dict(data_dict_1, package_method_1, data_dict_2, package_method_2, data_dict_3, package_method_3, data_dict_4, package_method_4, data_dict_5, package_method_5, "none", "none", "none", "none", "none", "none")
                        combined_dict_to_csv(data_dict_com, output_file_CD_masslist, output_file_all_data)
                else:
                    data_dict_com = combine_dict(data_dict_1, package_method_1, data_dict_2, package_method_2, data_dict_3, package_method_3, data_dict_4, package_method_4, "none", "none", "none", "none", "none", "none", "none", "none")
                    combined_dict_to_csv(data_dict_com, output_file_CD_masslist, output_file_all_data)
            else:
                data_dict_com = combine_dict(data_dict_1, package_method_1, data_dict_2, package_method_2, data_dict_3, package_method_3, "none", "none", "none", "none", "none", "none", "none", "none", "none", "none")
                combined_dict_to_csv(data_dict_com, output_file_CD_masslist, output_file_all_data)
        else:
            data_dict_com = combine_dict(data_dict_1, package_method_1, data_dict_2, package_method_2, "none", "none", "none", "none", "none", "none", "none", "none", "none", "none", "none", "none")
            combined_dict_to_csv(data_dict_com, output_file_CD_masslist, output_file_all_data)

###############################################################################################################################################################################################################################################################
print("Script finished successfully")
# end of script
###############################################################################################################################################################################################################################################################
# ░░░░░░░░███████████████░░░░░░░░
# ░░░░░█████████████████████░░░░░
# ░░░░████████████████████████░░░
# ░░░██████████████████████████░░
# ░░█████████████████████████████
# ░░███████████▀░░░░░░░░░████████
# ░░███████████░░░░░░░░░░░░░░░███
# ░████████████░░░░░░░░░░░░░░░░██
# ░█░░███████░░░░░░░░░░░▄▄░░░░░██
# █░░░░█████░░░░░░▄███████░░██░░█
# █░░█░░░███░░░░░██▀▀░░░░░░░░██░█
# █░░░█░░░░░░░░░░░░▄██▄░░░░░░░███
# █░░▄█░░░░░░░░░░░░░░░░░░█▀▀█▄░██
# █░░░░░░░░░░░░░░░░░░░░░░█░░░░██░
# ░███░░░░░░░░░░░░░░░░░░░█░░░░█░░
# ░░█░█░░░░░░░█░░░░░██▀▄░▄██░░░█░
# ░░█░█░░░░░░█░░░░░░░░░░░░░░░░░█░
# ░░░██░░░░░░█░░░░▄▄▄▄▄▄░░░░░░█░░
# ░░░██░░░░░░░█░░█▄▄▄▄░▀▀██░░█░░░
# ░░░██░░░░░░░█░░▀████████░░█░░░░
# ░░█░░█░░░░░░░█░░▀▄▄▄▄██░░█░░░░░
# ░░█░░░█░░░░░░░█░░░░░░░░░█░░░░░░
# ░█░░░░░█░░░░░░░░░░░░░░░░█░░░░░░
# ░░░░░░░░█░░░░░░█░░░░░░░░█░░░░░░
# ░░░░░░░░░░░░░░░░████████░░░░░░░
