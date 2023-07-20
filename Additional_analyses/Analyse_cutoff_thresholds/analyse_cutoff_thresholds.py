# Script by Jasmin Hafner, Eawag, October 2022
# Analyze TP prediction output to find maxTP / generation threshold that gives best precision
# Input:    1) Output from find_best_TPs.py
#           2) List of experimentally confirmed TPs (SMILES)

from rdkit import Chem
from rdkit.Chem.MolStandardize import rdMolStandardize
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np

# input
input_folder = 'input/'
output_folder = 'output/'
filename_list = ['TP_prediction_BBD_top_50.tsv',
                 'TP_prediction_BBD+SOIL_top_50.tsv',
                 'TP_prediction_BBD+SLUDGE_top_50.tsv',
                 'TP_prediction_BBD+SOIL+SLUDGE_top_50.tsv']
filename_PPS = 'results_EAWAG-PPS.tsv'
experimentally_confirmed = 'confirmed_TPs.txt'

color_palette = ['#B9252E', '#056586', '#8EC1B8', '#423E2D', '#C7AF5E']
# pps: #db4344, bbd: '#0054b6', 'bbd+soil': 'bbd+sludge':, 'bbd+soil+sludge:
# https://davidmathlogic.com/colorblind/#%23B9252E-%23056586-%238EC1B8-%23423E2D-%23C7AF5E

def __main__():
    confirmed_TPs = file_to_list(input_folder+experimentally_confirmed)
    TP_df = collect_data(confirmed_TPs)
    plot_data(TP_df)

def file_to_list(input_dir):
    open_file = open(input_dir)
    TP_list = []
    for line in open_file:
        TP_list.append(line.rstrip())
    return TP_list

def canonicalize_smiles(smiles):
    mol = Chem.MolFromSmiles(smiles)  # creates mol object from SMILES
    uncharger = rdMolStandardize.Uncharger() # easier to access
    mol = uncharger.uncharge(mol)  # protonates or deprotonates the mol object
    new_smiles = Chem.rdmolfiles.MolToSmiles(mol)  # converts mol object to canonical SMILES
    can_smiles = Chem.CanonSmiles(new_smiles)
    return can_smiles

def get_dictionary(TP_data, filename, confirmed_TPs):
    TP_file = open(input_folder + filename)
    for line in TP_file:
        if '/// Pathway name: ' in line:
            rank = 0
            continue
        elif 'SMILES' in line:
            continue

        # get data
        split = line.rstrip().split('\t')
        smiles = split[0]
        name = split[1]
        combined_probability = float(split[2])
        generation = int(split[4])
        probability = float(split[5])
        data_origin = filename.split('_')[2]
        rank+=1
        if not 'TP' in name:
            continue

        # check if TP has been experimentally confirmed
        can_smiles = canonicalize_smiles(smiles)
        confirmed = False
        if can_smiles in confirmed_TPs:
            confirmed = True

        # generate dictionary
        compound_data = {'smiles': can_smiles, 'name': name, 'combined_p': combined_probability, 'p': probability,
                         'generation': generation, 'confirmed': confirmed, 'data_origin': data_origin, 'rank': rank}

        TP_data.append(compound_data)
    return TP_data

def get_stats(df, x, column_name, origin):
    this_df = df.loc[df['data_origin'] == origin]
    precision = []
    count = []
    origin_list = []
    total_count_list = []
    for x_val in x:
        if column_name == 'generation':
            cutoff_data = this_df.loc[this_df[column_name]<=x_val]
        elif column_name in ['combined_p', 'p']:
            cutoff_data = this_df.loc[this_df[column_name] >= x_val]
        elif column_name == 'rank':
            cutoff_data = this_df.loc[this_df[column_name] <= x_val]

        if True not in cutoff_data['confirmed'].value_counts():
            true_count = 0
        else:
            true_count = cutoff_data['confirmed'].value_counts()[True]
        total_count = len(cutoff_data)

        if total_count == 0:
            prec = 0
        else:
            prec = true_count/total_count*100
        if origin == 'EAWAG-PPS':
            new_origin = origin
        else:
            new_origin = 'eP-'+origin


        precision.append(prec)
        count.append(true_count)
        origin_list.append(new_origin)
        total_count_list.append(total_count)
    return precision, count, origin_list, total_count_list

def plot_data_generation(df):
    max_val = df['generation'].max()
    x = list(range(0, max_val))
    gen_data = pd.DataFrame()
    data_origin_list = ['EAWAG-PPS', 'BBD', 'BBD+SOIL', 'BBD+SLUDGE', 'BBD+SOIL+SLUDGE']
    for origin in data_origin_list:
        y_prec, y_count, origin_list, total_count = get_stats(df, x,  'generation', origin)
        new_gen_data = pd.DataFrame({'Generation threshold': x, 'Precision (in %)': y_prec,
                                     'Correctly predicted TPs': y_count, 'Prediction method': origin_list,
                                     'Predicted compounds': total_count})
        all_true = df.loc[df['data_origin'] == origin]['confirmed'].value_counts()[True]
        all = len(df.loc[df['data_origin'] == origin])
        overall_precision = all_true/all
        gen_data = gen_data.append(new_gen_data)
        print('origin: {}, precision: {}, TPs_found: {}, total_tps: {}'.format(origin,overall_precision, all_true, all))
    # plot data
    column_name = 'Generation threshold'
    gen_data.to_csv('{}Cutoff_thresholds_{}.txt'.format(output_folder, column_name), sep='\t')
    fig, axs = plt.subplots(2, 2, figsize=(6, 6), sharey='row')
    sns.set_context('paper')
    gen_data = gen_data.loc[gen_data['Generation threshold']<=3]
    sns.lineplot(data=gen_data, x=column_name, y='Precision (in %)',
                 hue='Prediction method', ax=axs[1][0], palette=color_palette, legend=False)
    this = sns.lineplot(data=gen_data, x=column_name, y='Correctly predicted TPs', hue='Prediction method', ax=axs[0][0],
                 palette=color_palette, legend=True)
    sns.move_legend(this, "center left", bbox_to_anchor=(-0.03, 1.18), ncol=3, borderaxespad=0)
    return plt, axs


def plot_data_max_TPs(df, plt, axs):
    x = list(np.arange(0, 51, 1))
    gen_data = pd.DataFrame()
    data_origin_list = ['BBD', 'BBD+SOIL', 'BBD+SLUDGE', 'BBD+SOIL+SLUDGE']
    for origin in data_origin_list:
        y_prec, y_count, origin_list, total_count = get_stats(df, x, 'rank', origin)
        new_gen_data = pd.DataFrame({'Maximum TP threshold': x, 'Precision (in %)': y_prec,
                                     'Correctly predicted TPs': y_count, 'Prediction method': origin_list,
                                     'Predicted compounds': total_count})
        gen_data = gen_data.append(new_gen_data)
    column_name = 'Maximum TP threshold'

    sns.lineplot(data=gen_data, x=column_name, y='Precision (in %)', hue='Prediction method', ax=axs[1][1],
                 palette=color_palette[1:], legend=False)
    sns.lineplot(data=gen_data, x=column_name, y='Correctly predicted TPs', hue='Prediction method', ax=axs[0][1],
                 palette=color_palette[1:], legend=False)
    gen_data.to_csv('{}Cutoff_thresholds_{}.txt'.format(output_folder, column_name), sep='\t')
    # plt.tight_layout()
    plt.savefig('{}Plot_cutoff_analysis.pdf'.format(output_folder, column_name))
    plt.close()

def add_PPS_data(data, confirmed_TPs):
    pps_file = open(input_folder + filename_PPS)
    pps_file.readline() # header
    for line in pps_file:
        split = line.rstrip().split('\t')
        can_smiles = canonicalize_smiles(split[0])
        generation = int(split[1])
        confirmed = False
        if can_smiles in confirmed_TPs:
            confirmed = True
        compound = {'smiles': can_smiles, 'name': '', 'combined_p': np.NaN, 'p': np.NaN,
         'generation': generation, 'confirmed': confirmed, 'data_origin': 'EAWAG-PPS'}
        data.append(compound)
    return data

def collect_data(confirmed_TPs):
    TP_data = []  # {'smiles': {'name': '', 'p': '', 'combined_p': '' }}
    for filename in filename_list:
        TP_data = get_dictionary(TP_data, filename, confirmed_TPs)
    TP_data = add_PPS_data(TP_data, confirmed_TPs)
    TP_df = pd.DataFrame(TP_data)
    TP_df.to_csv(output_folder + 'TP_list_cutoff_analysis.txt', sep='\t')
    return TP_df

def plot_data(TP_df):
    plt, axs = plot_data_generation(TP_df)
    plot_data_max_TPs(TP_df, plt, axs)


__main__()
