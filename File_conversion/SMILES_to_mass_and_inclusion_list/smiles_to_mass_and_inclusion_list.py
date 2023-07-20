# -*- coding: utf-8 -*-
"""
Created on Wed Aug 17 09:25:44 2022

@author: trostele
"""

# start of script
# Python version 3.6.13
################################################################################################################################
# import all necessary packages
import pandas as pd
import rdkit # rdkit is only supported before Python 3.7
from rdkit import Chem
from rdkit.Chem.MolStandardize import rdMolStandardize
from rdkit.Chem import Descriptors
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
import pubchempy as pcp
import re
################################################################################################################################
# INPUT

smiles_file_location = "./smiles.txt" # add the path to the txt file that contains the SMILES

################################################################################################################################

# FUNCTIONS

def canonicalize_smiles(smiles):
    mol = Chem.MolFromSmiles(smiles)  # creates mol object from SMILES
    uncharger = rdMolStandardize.Uncharger() # easier to access
    uncharged = uncharger.uncharge(mol)  # protonates or deprotonates the mol object
    new_smiles = rdkit.Chem.rdmolfiles.MolToSmiles(uncharged)  # converts mol object to canonical SMILES
    can_smiles = Chem.CanonSmiles(new_smiles)
    return can_smiles

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

def smi_list_to_csv (smi_location):
    print("Converting {} file".format(smi_location))
    # create lists for each column and then create dict with correct which can be converted to dataframe using pandas
    SMILES_comp_file = open(smi_location)
    df_smi_list_com = []
    for line in SMILES_comp_file:
        df_smi_list_com.append(canonicalize_smiles(line.rstrip()))
 
    df_name_list_com = []
    df_ID_list_com = []
    df_Formula_list_com = []
    df_MolWeight_list_com = []
    df_Structure_list_com = []
    df_inchikey_list_com = []   
    # change the codes of the TPs, so each TP has its own name
    counter = 1
    for compound in df_smi_list_com:
        df_name_list_com.append("Compound_" + str(counter))
        counter = counter + 1
        df_MolWeight_list_com.append(Chem.Descriptors.ExactMolWt(Chem.MolFromSmiles(compound)))
        df_Formula_list_com.append(CalcMolFormula(Chem.MolFromSmiles(compound)))
        df_inchikey_list_com.append(Chem.inchi.MolToInchiKey(Chem.MolFromSmiles(compound)))      
        df_ID_list_com.append(get_cas_inchi(Chem.inchi.MolToInchiKey(Chem.MolFromSmiles(compound))))
        mol_rep_1 = (Chem.MolToMolBlock(Chem.MolFromSmiles(compound)).replace("\n",";"))
        mol_rep_2 = mol_rep_1[6:]
        mol_rep_3 = mol_rep_2[:-1] 
        df_Structure_list_com.append(mol_rep_3)


    # all lists must be the same length to convert it to dataframe
    assert len(df_name_list_com) == len(df_Formula_list_com) == len(df_MolWeight_list_com) == len(df_Structure_list_com) == len(df_ID_list_com), "Error: all lists must be the same length to convert it to dataframe"

    #create dict and convert to dataframe
    df_com_dict = {"Name": df_name_list_com, "ID": df_ID_list_com,"Formula": df_Formula_list_com,"MolWeight": df_MolWeight_list_com, "Structure": df_Structure_list_com}
    df_com = pd.DataFrame.from_dict(df_com_dict)
    # export dataframe as csv
    df_com.to_csv("./CD_masslist.csv", index = False, sep = "\t")
    df_complete_dict = {"SMILES": df_smi_list_com, "Name": df_name_list_com,"CAS": df_ID_list_com,"Formula": df_Formula_list_com,"MolWeight": df_MolWeight_list_com, "InchiKey":df_inchikey_list_com}
    df_complete = pd.DataFrame.from_dict(df_complete_dict)
    # export dataframe as csv
    df_complete.to_csv("./overview.csv", index = False, sep = ",")
    
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
    df_inclusion_pos.to_csv("./inclusion_list_pos.csv", index = False, sep = ",")  
    d_inclusion_neg = {"Mass [m/z]": df_M_minus_H ,"Formula [M]": df_empty,
                   "Formula type": df_empty, "Species": df_empty,
                   "CS [z]": df_empty, "Polarity": df_polarity_neg, "Start [min]": df_empty, "End [min]": df_empty,
                   "(N)CE": df_nce, "(N)CE type":df_nce_type, "MSX ID": df_empty, "Comment": df_name_list_com}
    df_inclusion_neg = pd.DataFrame.from_dict(d_inclusion_neg)
    # export dataframe as csv
    df_inclusion_neg.to_csv("./inclusion_list_neg.csv", index = False, sep = ",") 
    print("Export complete")
    smi_list = df_smi_list_com
    return smi_list
    
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
            
    with open('./max_element_count.txt', 'w') as f:
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

################################################################################################################################

# MAIN

# enter the pathway to the smiles text file below
smi_list = smi_list_to_csv(smiles_file_location)
# prints out max element count and creates txt file
max_element_count(smi_list)
# prints out suggested values for stepped NCE approach if possible
suggest_stepped_nce(smi_list)

################################################################################################################################
print("Script finished successfully")
# end of script
################################################################################################################################