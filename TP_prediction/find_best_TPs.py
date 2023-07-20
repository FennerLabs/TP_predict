import sys
sys.path.insert(0, '../src/envipath-python/enviPath_python/')
sys.path.insert(0, '../src/envipath-python/')
from enviPath_python.enviPath import *
from enviPath_python.enviPath import *
from util import *
import getpass

#---------------------------#
# enviPath SETTINGS         #
#---------------------------#

# Define the instance to use
INSTANCE_HOST = 'https://envipath.org'
USERNAME = '' # TODO USER: enter your username

# MODEL # TODO USER: Select model for relative reasoning. Default: Standard model BBD - ML - ECC - 2022
# New enviPath - default relative reasoning model on envipath: BBD - ML - ECC - 2022
# EP_MODEL_ID = 'https://envipath.org/package/de0cdca1-c3ff-44ed-8ffd-f29c269bfa55/relative-reasoning/646afb6c-6cfc-4d4b-8d22-e196d849ec73'
# New enviPath - BBD+SOIL relative reasoning model on envipath: BBD+SOIL - ML - ECC - 2022
# EP_MODEL_ID = 'https://envipath.org/package/de0cdca1-c3ff-44ed-8ffd-f29c269bfa55/relative-reasoning/0aec0115-941d-4e33-a844-4431d3ec598d'
# New enviPath - BBD+SOIL relative reasoning model on envipath: BBD+SLUDGE - ML - ECC - 2022
# EP_MODEL_ID = 'https://envipath.org/package/de0cdca1-c3ff-44ed-8ffd-f29c269bfa55/relative-reasoning/2a41e599-f962-4268-b7d1-5f1cb252b937'
# New enviPath - BBD+SOIL relative reasoning model on envipath: BBD+SOIL+SLUDGE - ML - ECC - 2022
# EP_MODEL_ID = 'https://envipath.org/package/de0cdca1-c3ff-44ed-8ffd-f29c269bfa55/relative-reasoning/76bd8654-e02f-4fbd-98df-bf61411f9b92'
# Old model - old default relative reasoning model on envipath from 2021: BBD - ML - ECC
EP_MODEL_ID = 'https://envipath.org/package/32de3cf4-e3e6-4168-956e-32fa5ddb0ce1/relative-reasoning/edaf8d8c-430a-4277-848b-3e163a86febf'

# PACKAGE - # TODO USER: prepare a new package (manually) and add it's URI here - make sure it is empty when running script!
EP_PACKAGE_ID = 'https://envipath.org/package/de0cdca1-c3ff-44ed-8ffd-f29c269bfa55' # Test package
# List of output packages used for Sludge TP paper
# EP_PACKAGE_ID = 'https://envipath.org/package/0915fad3-b889-4aa8-ac98-0707b717be57' # Package for results using BBD - ML - ECC - 2022 model
# EP_PACKAGE_ID = 'https://envipath.org/package/80cf58b1-21e2-4c28-9cc6-dc69c6445bdf' # Package for results using BBD+SOIL - ML - ECC - 2022 model
# EP_PACKAGE_ID = 'https://envipath.org/package/7d64aa85-2e3c-413f-a538-4d5f2bfd4662' # Package for results using BBD+SLUDGE - ML - ECC - 2022 model
# EP_PACKAGE_ID = 'https://envipath.org/package/11f2acd5-5209-4d49-ad77-93f6f6965886' # Package for results using BBD+SOIL+SLUDGE - ML - ECC - 2022 model


#---------------------------#
# PATHWAY SEARCH SETTINGS   #
#---------------------------#
# These are the default settings used for the Sludge TP paper.
# They can be modified to direct the pathway search towards a specific objective.

# Maximum number of TPs to predict
MAX_TP = 20

# Lower probability threshold
PROBABILITY_THRESHOLD = 0 # any value equal to or lower than the threshold will be excluded

# Set probabilities of 0 to 0.01 to continue having a weighting scheme downstream of the pathway
INCLUDE_0_PROBABILITIES = False

# Follow moiety - only compounds containing moiety in SMILES will be expanded
MOIETY = "" # e.g., "C(F)(F)F"

# To prioritize small compounds in the queue
SORT_TPS_BY_SIZE = False

# Follow labeled atoms
FOLLOW_LABELED_ATOM = False
ATOM_LABEL = '14'

# Print as a reaction file additionally to the list of TPs
PRINT_REACTION_FILE = True

#---------------------------#
# FILE PATH SETTINGS        #
#---------------------------#
# Input/output files
INPUT_FILE_PATH = 'input/input_structures.tsv'
OUTPUT_DIRECTORY = 'output/'
OUTPUT_FILE_TAG = 'TEST'


#---------------------------#
# CONNECT TO ENVIPATH       #
#---------------------------#
eP = enviPath(INSTANCE_HOST)
password = getpass.getpass()
eP.login(USERNAME, password)


#---------------------------#
# FUNCTIONS                 #
#---------------------------#
def __main__(rr_id, pkg_id, tag):
    """
    Main function, predicts pathways for a list of input smiles
    Output: pathways are saved to specified enviPath package, TP list as .tsv file to output folder
    :param rr_id: URI of enviPath relative reasoning mode to be used
    :param pkg_id: URI of enviPath package to store resulting pathways
    :param tag: string tag to attach to output files for identification
    """
    rr = RelativeReasoning(eP.requester, id=rr_id)
    pkg = Package(eP.requester, id=pkg_id)
    input_list = load_input(INPUT_FILE_PATH)
    output_file_path_TP = '{}TP_prediction_{}_top_{}.tsv'.format(OUTPUT_DIRECTORY, tag, MAX_TP)
    output_file_path_reaction = '{}reaction_prediction_{}_top_{}.tsv'.format(OUTPUT_DIRECTORY, tag, MAX_TP)
    outfile_TP = open(output_file_path_TP, 'w')
    outfile_reaction = open(output_file_path_reaction, 'w')
    rule_dict = {}
    for compound_input in input_list:
        result = predict_TPs(compound_input['smiles'], compound_input['name'], rr)
        result_list = clean_result(result) # sort and name TPs
        pathway_info = upload_envipath_pathway(eP, result_list, pkg)
        print_TP_file(outfile_TP, result_list, pathway_info) # continuous writing of result file
        if PRINT_REACTION_FILE:
            rule_dict = print_reaction_file(outfile_reaction, result_list, rule_dict, pathway_info)



def print_TP_file(open_TP_file, result, pathway = None):
    """
    Prints output to open file
    :param open_file: writable file object
    :param result: clean result dictionary from predict_TPs() function
    :param pathway: pathway URI from upload_envipath_pathway() function
    """
    open_TP_file.write('///') # signifies new pathway entry
    if pathway:
        open_TP_file.write(' Pathway name: {}, Pathway id: {}'.format(pathway['name'], pathway['id']))
    open_TP_file.write('\n')
    header= ('SMILES\tname\tcombined_probability\trules\tgeneration\tprobability\tparent\n')
    open_TP_file.write(header)
    for TP in result:
        open_TP_file.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(TP['smiles'],
                                              TP['name'],
                                              TP['combined_probability'],
                                              TP['rules'],
                                              TP['rule_IDs'],
                                              TP['generation'],
                                              TP['probability'],
                                              TP['parent_smiles'])
                        )

def print_reaction_file(open_file, result, rule_dict, pathway=None):
    open_file.write('///') # signifies new pathway entry
    if pathway:
        open_file.write(' Pathway name: {}, Pathway id: {}'.format(pathway['name'], pathway['id']))
    open_file.write('\n')
    header= ('substrate_name\tproduct_name\tsubstrate_SMILES\tproduct_SMILES\trules\tEC_numbers\tgeneration\tcombined_probability\tprobability\n')
    open_file.write(header)
    for TP in result[1:]: # ignore first line which is parent
        EC_list, rule_dict = get_EC_list(TP, rule_dict)
        open_file.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(
                                              TP['parent_name'],
                                              TP['name'],
                                              TP['parent_smiles'],
                                              TP['smiles'],
                                              TP['rules'],
                                              ','.join(EC_list),
                                              TP['generation'],
                                              TP['combined_probability'],
                                              TP['probability']))
    return rule_dict

def get_EC_list(TP, rule_dict):
    rule_id_list = TP['rule_IDs'].split(',')
    EC_list = []
    for rule_id in rule_id_list:
        if not rule_dict.get(rule_id):
            ec_list_for_rule = []
            rule = ParallelCompositeRule(eP.requester, id = rule_id)
            try:
                rule_ECs = rule._get('ecNumbers')
            except:
                print('Warning: JsonDecodeError for ', rule_id)
            else:
                for ec in rule_ECs:
                    if 'KEGG' in ec['linkingMethod']:
                        ec_list_for_rule.append(ec['ecNumber'])
                rule_dict[rule_id] = ec_list_for_rule
        else:
            ec_list_for_rule = rule_dict[rule_id]

        EC_list.extend(ec_list_for_rule)
    return EC_list, rule_dict


def predict_TPs(input_smiles, input_name, rr):
    """
    Pathway prediction for single compound
    :param input_smiles: input smiles of parent compound
    :param input_name: name of parent compound
    :param rr: relative reasoning object
    :return: dictionary of resulting TPs
    """
    print('\n### PREDICT TPs FOR COMPOUND {} ###\n'.format(input_name))
    num_TP = -1 # counter starts at -1, because source compound is also in the TP list
    validated_TPs = {}  # container for resulting predictions
    queued_items = [{'probability': 1, 'combined_probability': 1, 'smiles': input_smiles, 'generation': 0, 'parent_smiles': '',
                     'rules': '', 'rule_IDs': '', 'name': input_name, 'size': len(input_smiles)}]
    queue = [input_smiles]  # queue is updated after each cycle to have top TP first, list of smiles
    while num_TP < MAX_TP:
        if len(queue) == 0:
            print('\nEmpty queue - The exploration of has converged at {} predicted TPs'.format(num_TP))
            return validated_TPs # stop TP prediction
        smiles = queue.pop(0) # get top item in queue
        data = queued_items.pop(0) # remove data from queued items
        result_list = expand_smiles(smiles, rr)  # create children
        TP_dict = result_to_compound_dict(result_list)
        queue, queued_items, validated_TPs = update_queue(queue, queued_items, validated_TPs, TP_dict, data)
        validated_TPs[smiles] = data
        num_TP += 1
    return validated_TPs


def update_queue(_queue,_queued_items, _validated_TPs, _TPs, _parent_data):
    """
    Update queue with TPs predicted in current iteration
    :param _queue: ordered list of smiles to explore
    :param _queued_items: ordered list of compound dictionaries, same order as _queue
    :param _validated_TPs: list of already validated TPs for resulting pathway
    :param _TPs: predicted TPs from current iteration, to be evaluated and added to queue
    :param _parent_data: compound dictionary of the parent compound of _TPs
    :return: new_queue: new ordered list of smiles to explore
    :return _queued_items: new ordered list of compound dictionaries
    :return: _validated_TPs: updated list of already validated TPs
    """
    parent_probability = _parent_data['combined_probability']
    parent_generation = _parent_data['generation']
    parent_smiles = _parent_data['smiles']
    queue_before = len(_queue)
    for smiles in _TPs:
        data = _TPs[smiles]
        # If the probability is 0 , we don't consider the TP further
        this_probability = data['probability']
        if this_probability <= PROBABILITY_THRESHOLD:
            continue
        # If a moiety is given and it is not in SMILES, we don't follow the TP further
        if MOIETY not in smiles:
            continue
        if FOLLOW_LABELED_ATOM and ATOM_LABEL not in smiles:
            continue
        if INCLUDE_0_PROBABILITIES and this_probability == 0:
            this_probability = 0.01
        # add combined probability
        this_combined_probability = parent_probability * this_probability
        this_generation = parent_generation + 1
        rules = data['rules']
        rule_IDs = data['rule_IDs']
        # first, check if compound already in validated. if yes, update
        if smiles in _validated_TPs.keys():
            _validated_TPs[smiles] = update_compound_entry(_validated_TPs[smiles],
                                                           this_combined_probability, rules, rule_IDs,
                                                           this_generation, parent_smiles, size_metric='size',
                                                           size_value=len(smiles))
        # next, check if compound is already in queue. if yes, update
        elif smiles in _queue:
            index = _queue.index(smiles)
            assert smiles == _queued_items[index]['smiles'], \
                'smiles {} does not match smiles in {}'.format(smiles, _queued_items[index])
            _queued_items[index] = update_compound_entry(_queued_items[index],
                                                           this_combined_probability, rules, rule_IDs,
                                                           this_generation, parent_smiles, size_metric='size',
                                                           size_value=len(smiles))
        # else, add new item to queue
        else:
            data['combined_probability'] = this_combined_probability
            data['generation'] = this_generation
            data['parent_smiles'] = parent_smiles
            data['carbon_count'] = smiles.upper().count('C')
            _queued_items.append(data)
            _queue.append(smiles)
            assert len(_queued_items) == len(_queue)
    # First sort by size
    if SORT_TPS_BY_SIZE:
        _queued_items.sort(reverse=False, key=lambda x: x['size'])
    # order dict by combined probability
    _queued_items.sort(reverse=True, key=lambda x: x['combined_probability'])
    queue_after = len(_queue)
    print('Added {} smiles to queue'.format(queue_after - queue_before))
    new_queue = [] # resetting queue
    [new_queue.append(x['smiles']) for x in _queued_items]
    print ('New queue for compound', parent_smiles)
    for q in new_queue:
        print(q, _queued_items[new_queue.index(q)]['combined_probability'])

    return new_queue, _queued_items, _validated_TPs


#---------------------------#
# MAIN                      #
#---------------------------#
__main__(rr_id = EP_MODEL_ID, pkg_id = EP_PACKAGE_ID, tag = OUTPUT_FILE_TAG)


