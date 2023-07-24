from enviPath_python import enviPath
from enviPath_python.objects import *

def load_input(input_path):
    """
    Load input smiles and names
    :param input_path: path to input file
    :return: list of dictionaries containing smiles and name of input compounds
    """
    f = open(input_path)
    input_list = []
    for line in f:
        if line not in ['', '\n']:
            line_split = line.rstrip().split('\t')
            input_list.append({'smiles': line_split[0], 'name': line_split[1]})
    return input_list

def upload_envipath_pathway(eP, result, pkg):
    """
    Upload resulting pathway dictionary to enviPath
    :param eP: enviPath object
    :param result: result list of dictionaries from pathway prediction
    :param pkg: package object where results should be uploaded
    :return: dictionary {'name': pathway name, 'id': URI of pathway}
    """
    assert 'anonymous' not in str(eP.who_am_i()), 'Upload not possible when not logged in'
    source = result[0]
    pkg.add_compound(smiles=source['smiles'],name=source['name'])
    pathway = Pathway.create(pkg, smiles=source['smiles'], name=source['name'], root_node_only=True)
    # Add the observed degradation product as a second node
    for TP in result[1:]:
        # print('adding to pw:', TP['name'], TP['smiles'], TP['generation'], TP['parent_smiles'])
        pathway.add_node(smiles=TP['smiles'], name=TP['name'], depth=TP['generation'])
        # check for the case of multiple parents:
        for parent in TP['parent_smiles'].split(','):
            pathway.add_node(smiles=parent)
            pathway.add_edge(smirks='{}>>{}'.format(parent, TP['smiles']))

    print('New pathway created for {}: {}'.format(source['name'], pathway.id))
    return {'name': source['name'], 'id': pathway.id}

def expand_smiles(smiles, rr):
    """
    Get all potential TPs by applying enviPath biotransformation and relative reasoning rules
    :param smiles: input smiles
    :param rr: relative reasoning object
    :return: list of dictionaries for each predicted TP: {'smiles': smiles,
                                                          'name': rule name,
                                                          'probability': relative reasoning probability}
    """
    res = rr.classify_smiles(smiles)
    # sort by probability
    res.sort(reverse=True, key=lambda x: x['probability'])
    return res

def clean_result(result_dict):
    """
    Sorts TP list for output
    :param result_dict: result dictionary
    :return: sorted and named list of TPs
    """
    result_list = list(result_dict.values())
    result_list.sort(reverse=True, key=lambda x: x['combined_probability'])
    result_list.sort(key=lambda x: x['generation']) # make sure that source compound is first
    # get name of source compound
    source_name = result_list[0]['name']
    source_smiles = result_list[0]['smiles']
    TP_count = 0
    D = {source_smiles: source_name}
    new_result_list = []
    new_result_list.append(result_list[0])
    for res in result_list[1:]:
        new_res = res
        TP_count += 1
        new_name = 'TP_{}_{}'.format(source_name, TP_count)
        new_res['name'] = new_name
        multiple_parents = res['parent_smiles'].split(',')
        for p in multiple_parents:
            new_res['parent_name'] = D[p]
            D[res['smiles']] = new_name
        new_result_list.append(new_res)
    return new_result_list

def result_to_compound_dict(result):
    """
    Translates result from enviPath node expansion into a compound dictionary
    :param result: list of dictionaries with predicted TP information
    :return: dictionary of TPs
    """
    compound_dict = {}
    for r in result:
        probability = float(r['probability'])
        for product_smiles in r['products']:
            if product_smiles not in compound_dict.keys():
                compound_dict[product_smiles] = {'rules' : r['name'], 'rule_IDs': r['id'], 'probability': probability, 'smiles': product_smiles}
            else:
                # check if there's a rule with better probability
                if probability > compound_dict[product_smiles]['probability']:
                    # update probability and rules associated to this probability
                    compound_dict[product_smiles]['probability'] = probability
                    compound_dict[product_smiles]['rules'] = r['name']
                    compound_dict[product_smiles]['rule_IDs'] = r['id']
    return compound_dict

def update_compound_entry(compound_entry, this_combined_probability, rules, rule_IDs, this_generation, parent_smiles,
                          size_metric, size_value):
    """
    Update the compound entry with new information
    :param compound_entry: dictionary of compound information
    :param this_combined_probability: new combined probability
    :param rules: new rules
    :param rule_IDs: new rule IDs
    :param this_generation: new generation
    :param parent_smiles: new parent compound
    :param parent_compound: new parent compound
    :param size_metric: size metric
    :param size_value: new size value
    :return: updated compound entry
    """
    if compound_entry['combined_probability'] < this_combined_probability:
        compound_entry['combined_probability'] = this_combined_probability
        compound_entry['rules'] = rules
        compound_entry['rule_IDs'] = rule_IDs
        compound_entry['generation'] = this_generation
        compound_entry['parent_smiles'] = parent_smiles
        compound_entry[size_metric] = size_value
    elif compound_entry['combined_probability'] == this_combined_probability:
        compound_entry['rules'] += ',{}'.format(rules)
        compound_entry['rule_IDs'] += ',{}'.format(rule_IDs)
        compound_entry['parent_smiles'] += ',{}'.format(parent_smiles)
    return compound_entry
