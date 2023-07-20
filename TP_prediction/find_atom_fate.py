from enviPath_python.enviPath import *
from enviPath_python.enviPath import *
from util import *
import getpass
from rdkit import Chem
from rdkit.Chem import Descriptors

#---------------------------#
# enviPath SETTINGS         #
#---------------------------#
# enviPath instance - Define the instance to use
INSTANCE_HOST = 'https://envipath.org'
USERNAME = 'anonymous' # TODO USER: enter your username

# MODEL - Define model for relative reasoning and output package
# New enviPath - default relative reasoning model on envipath: BBD - ML - ECC - 2022
EP_MODEL_ID = 'https://envipath.org/package/de0cdca1-c3ff-44ed-8ffd-f29c269bfa55/relative-reasoning/646afb6c-6cfc-4d4b-8d22-e196d849ec73'

# PACKAGE - prepare a new package (manually) and add it's URI here - make sure it is empty / new  when running script!
EP_PACKAGE_ID = '' # TODO USER: enter the URI of your results package here to save results to enviPath


#---------------------------#
# PATHWAY SEARCH SETTINGS   #
#---------------------------#
# Maximum number of TPs to predict
MAX_TP = 50

# Lower probability threshold
PROBABILITY_THRESHOLD = -1 # any value equal to or lower than the threshold will be excluded

# Set probabilities of 0 to 0.01 to continue having a weighting scheme downstream of the pathway
INCLUDE_0_PROBABILITIES = False

# Follow moiety - only compounds containing moiety in SMILES will be expanded
MOIETY = "" # e.g., "C(F)(F)F"

# To prioritize small compounds in the
SORT_TPS_BY_SIZE = True
SIZE_METRIC = 'carbon_count' # other options: 'size', 'mw'

# Follow labeled atoms
FOLLOW_LABELED_ATOM = True
ATOM_LABEL = '14'
EXPLORE_NON_ZERO_PATHWAYS = False

#---------------------------#
# FILE PATH SETTINGS        #
#---------------------------#
# Input/output files
INPUT_FILE_PATH = '../input/input_structures.tsv'
OUTPUT_DIRECTORY = '../output/'
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
    rr = RelativeReasoning(eP.requester, id=rr_id)
    pkg = Package(eP.requester, id=pkg_id)
    input_list_raw = load_input(INPUT_FILE_PATH)
    input_list = label_input_list(input_list_raw)
    output_file_path = '{}TP_prediction_{}_top_{}.tsv'.format(OUTPUT_DIRECTORY, tag, MAX_TP)
    outfile = open(output_file_path, 'w')
    for compound_input in input_list:
        result = find_path_to_stable_moiety(compound_input['smiles'], compound_input['name'], rr)
        result_list = clean_result(result) # sort and name TPs
        pathway_info = upload_envipath_pathway(eP, result_list, pkg)
        print_result(outfile, result_list, pathway_info) # continuous writing of result file

def label_input_list(input_list):
    isotopomer_list = []
    for index, entry in enumerate(input_list):
        mol = Chem.MolFromSmiles(entry['smiles'])
        # Iterate over the atoms
        c_count = 1
        for atom in mol.GetAtoms():
            if atom.GetSymbol() != 'C':
                continue
            # For each atom, set the property "atomLabel" to a custom value, let's say a1, a2, a3,...
            atom.SetIsotope(14)
            smiles = Chem.MolToSmiles(mol, kekuleSmiles=True)
            name = entry['name']+'_C{}'.format(c_count)
            isotopomer_list.append({'smiles': smiles, 'name': name})
            # reset structure
            atom.SetIsotope(False)
            # mol = Chem.Cleanup(mol)
            c_count +=1
    return isotopomer_list

def remove_isomer_information(smiles):
    mol = Chem.MolFromSmiles(smiles)
    for atom in mol.GetAtoms():
        atom.SetIsotope(False)
    # mol = Chem.Cleanup(mol)
    new_smiles = Chem.MolToSmiles(mol, kekuleSmiles=True)
    return new_smiles

def get_molecular_weight(smiles):
    mol = Chem.MolFromSmiles(smiles)
    mw = Descriptors.MolWt(mol)
    return mw

def get_carbon_count(smiles):
    mol = Chem.MolFromSmiles(smiles)
    # Iterate over the atoms to count carbons
    c_count = 0
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'C':
            c_count += 1
    return c_count

def get_size_value(smiles):
    if SIZE_METRIC == 'size':
        val = len(smiles)
    elif SIZE_METRIC == 'mw':
        val = get_molecular_weight()
    elif SIZE_METRIC == 'carbon_count':
        val = get_carbon_count(smiles)
    else:
        raise ValueError("SIZE_METRIC needs to be set to size, mw, or carbon_count")
    return val

def print_result(open_file, result, pathway = None):
    header = ('Name\tSMILES\tID\tPersistent_moiety\tPersistent_moiety_unlabeled\tcombined_probability\tgeneration\n')
    open_file.write(header)
    smallest_TP = result[0]
    for TP in result:
        if TP[SIZE_METRIC] < smallest_TP[SIZE_METRIC]:
            smallest_TP = TP
    smallest_TP_unlabeled = remove_isomer_information(smallest_TP['smiles'])
    open_file.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(pathway['name'],
                                                          result[0]['smiles'],
                                                          pathway['id'],
                                                          smallest_TP['smiles'],
                                                          smallest_TP_unlabeled,
                                                          smallest_TP['combined_probability'],
                                                          smallest_TP['generation'])
                    )

def find_path_to_stable_moiety(input_smiles, input_name, rr):
    print('\n### FIND PATH FOR COMPOUND {} ###\n'.format(input_name))
    num_TP = -1 # counter starts at -1, because source compound is also in the TP list
    validated_TPs = {}  # container for resulting predictions
    queued_items = [{'probability': 1, 'combined_probability': 1, 'smiles': input_smiles, 'generation': 0, 'parent': '',
                     'rules': '', 'name': input_name, SIZE_METRIC: get_size_value(input_smiles)}]
    queue = [input_smiles]  # queue is updated after each cycle to have top TP first, list of smiles
    current_size = 10000
    current_stable_moiety = ''
    secondary_exploration = False
    while num_TP < MAX_TP:
        if len(queue) == 0:
            print('\nEmpty queue - The exploration of has converged at {} predicted TPs'.format(num_TP))
            return validated_TPs # stop TP prediction
        smiles = queue.pop(0) # get top item in queue
        data = queued_items.pop(0) # get data from queued items for new TP to explore
        # if new smiles to expand is bigger than previous smiles -> return
        if data[SIZE_METRIC] > current_size and data['combined_probability'] == 0:
            # if any compounds in the queue still have combined p > 0, then explore those:
            if EXPLORE_NON_ZERO_PATHWAYS:
                queue, queued_items = remove_zero_probabilities(queue, queued_items)
                secondary_exploration = True
                smiles = queue.pop(0)  # get top item in queue
                data = queued_items.pop(0)  # get data from queued items for new TP to explore
            # if new smiles to expand is bigger than previous smiles -> return
            if data[SIZE_METRIC] > current_size and data['combined_probability'] == 0:
                print('\nNo more smaller molecules in queue '
                      '- The exploration of has converged at {} predicted TPs'.format(num_TP))
                print('\nSMILES: {}, Stable moiety: {}, Number of TPs: {}'.format(smiles,
                                                                                  current_stable_moiety,
                                                                                  num_TP))
                return validated_TPs # stop TP prediction
        result_list = expand_smiles(smiles, rr)  # create children
        TP_dict = result_to_compound_dict(result_list)
        queue, queued_items, validated_TPs = update_queue(queue, queued_items, validated_TPs,
                                                          TP_dict, data, secondary_exploration)
        validated_TPs[smiles] = data
        current_size = data[SIZE_METRIC]
        current_stable_moiety = smiles
        num_TP += 1
    return validated_TPs

def remove_zero_probabilities(_queue, _queued_items):
    print("clean queue...")
    queue = []
    queued_items = []
    for entry in _queued_items:
        if entry['combined_probability'] > 0:
            queue.append(entry['smiles'])
            queued_items.append(entry)
    return queue, queued_items

def update_queue(_queue,_queued_items, _validated_TPs, _TPs, _parent_data, _secondary_exploration):
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
        if _secondary_exploration and this_probability == 0:
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
        # first, check if compound already in validated. if yes, update
        if smiles in _validated_TPs.keys():
            _validated_TPs[smiles] = update_compound_entry(_validated_TPs[smiles],
                                                           this_combined_probability, rules,
                                                           this_generation, parent_smiles, size_metric=SIZE_METRIC,
                                                           size_value=get_size_value(smiles))

        # next, check if compound is already in queue. if yes, update
        elif smiles in _queue:
            index = _queue.index(smiles)
            assert smiles == _queued_items[index]['smiles'], \
                'smiles {} does not match smiles in {}'.format(smiles, _queued_items[index])
            _queued_items[index] = update_compound_entry(_queued_items[index],
                                                           this_combined_probability, rules,
                                                           this_generation, parent_smiles, size_metric=SIZE_METRIC,
                                                           size_value=get_size_value(smiles))

        # else, add new item to queue
        else:
            data['combined_probability'] = this_combined_probability
            data['generation'] = this_generation
            data['parent'] = parent_smiles
            data[SIZE_METRIC] = get_size_value(smiles)
            _queued_items.append(data)
            _queue.append(smiles)
            assert len(_queued_items) == len(_queue)

    # first order dict by combined probability
    _queued_items.sort(reverse=True, key=lambda x: x['combined_probability'])
    # then sort by size
    if SORT_TPS_BY_SIZE:
        _queued_items.sort(reverse=False, key=lambda x: x[SIZE_METRIC])

    queue_after = len(_queue)
    print('Added {} smiles to queue'.format(queue_after - queue_before))
    new_queue = [] # resetting queue
    [new_queue.append(x['smiles']) for x in _queued_items]
    print ('New queue for compound', parent_smiles)
    for q in new_queue:
        print(q, _queued_items[new_queue.index(q)]['combined_probability'],
              _queued_items[new_queue.index(q)][SIZE_METRIC])

    return new_queue, _queued_items, _validated_TPs


#---------------------------#
# MAIN                      #
#---------------------------#
__main__(rr_id = EP_MODEL_ID, pkg_id = EP_PACKAGE_ID, tag = OUTPUT_FILE_TAG)

