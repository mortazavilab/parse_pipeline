import pandas as pd
import os

# from the parse pipeline
def load_bc_dict(fname, verb=False):
	""" Load barcode edit dict
	"""
	with open(fname, 'r') as INFILE:
	    bc_dict = json.load(INFILE)
	    if verb:
	        print(f"Loaded {fname}")

	# Top level has int keys and holds default dicts
	new_dict = {}
	for k, v in bc_dict.items():
	    new_dict[int(k)] = defaultdict(list, bc_dict[k])

	return new_dict

def get_bc_round_set(kit, chemistry):
	KIT_INT_DICT = {'custom_1': 1, 'WT': 48, 'WT_mini': 12, 'WT_mega': 96}
	kit_n = KIT_INT_DICT[kit]
	if kit_n == 12:
		bc_round_set = [['bc1','n24_v4'], ['bc2','v1'], ['bc3','v1']]
	if kit_n == 96:
		bc_round_set = [['bc1','n192_v4'], ['bc2','v1'], ['bc3','v1']]
	if kit_n == 48:
		bc_round_set = [['bc1','v2'], ['bc2','v1'], ['bc3','v1']]

	if kit == 'WT' and chemistry == 'v2':
		bc_round_set = [['bc1', 'n96_v4'], ['bc2', 'v1'], ['bc3', 'v1']]

	return bc_round_set

def load_barcodes(kit, chemistry):
	"""
	Load the barcodes. Adapted from the Parse biosciences pipeline.

	Returns:
		edit_dict_set (dict): Dict for barcode<1,2,3> with
			key: query bc
			item: corrected bc
	"""
	pkg_path = os.path.dirname(__file__)
	bc_path = '/'.join(pkg_path.split('/')[:-1])+'/barcodes/'

	bc_round_set = get_bc_round_set(kit, chemistry)

	edit_dict_set = {}
	for entry in bc_round_set:
		bc = entry[0]
		ver = entry[1]
		fname = bc_path + 'bc_dict_{}.json'.format(ver)
		edit_dict = load_bc_dict(fname)
		edit_dict_set[bc] = edit_dict

	return edit_dict_set

# From the Parse biosciences pipeline
def load_barcodes_set(kit, chemistry):
	"""
	Load the barcodes. Adapted from the Parse biosciences pipeline.
	"""
	pkg_path = os.path.dirname(__file__)
	bc_path = '/'.join(pkg_path.split('/')[:-1])+'/barcodes/'

	bc_round_set = get_bc_round_set(kit, chemistry)

	bc_set = {}
	for entry in bc_round_set:
		bc = entry[0]
		ver = entry[1]
		fname = bc_path + '/bc_data_{}.csv'.format(ver)
		bc_df = pd.read_csv(fname)
		bcs = set(bc_df.sequence.tolist())
		bc_set[bc] = bcs

	return bc_set

def get_bc1_matches(kit, chemistry):
    pkg_path = os.path.dirname(os.path.dirname(__file__))
    # pkg_path = '/'.join(pkg_path.split('/')[:-1])
    bc_round_set = get_bc_round_set(kit, chemistry)

    # determine file to use
    for entry in bc_round_set:
        if entry[0] == 'bc1':
            ver = entry[1]

    # read in and restructure such that each dt bc is
    # matched with its randhex partner from the same well
    fname = pkg_path+'/barcodes/bc_data_{}.csv'.format(ver)
    df = pd.read_csv(fname)
    df.loc[df.well.duplicated(keep=False)].sort_values(by='well')
    drop_cols = ['bci', 'uid', 'type']
    bc1_dt = df.loc[df['type'] == 'T'].drop(drop_cols, axis=1)
    bc1_dt.rename({'sequence': 'bc1_dt'}, axis=1, inplace=True)
    bc1_randhex = df.loc[df['type'] == 'R'].drop(drop_cols, axis=1)
    bc1_randhex.rename({'sequence': 'bc1_randhex'}, axis=1, inplace=True)
    bc_df = bc1_dt.merge(bc1_randhex, on='well')

    return bc_df

def get_bcs(bc, kit, chemistry):
    """
    Parameters:
        bc (int): {1,2,3}
    """
    pkg_path = os.path.dirname(os.path.dirname(__file__))
    bc_round_set = get_bc_round_set(kit, chemistry)

    bc_name = f'bc{bc}'

    # determine file to use
    for entry in bc_round_set:
        if entry[0] == bc_name:
            ver = entry[1]

    fname = pkg_path+'/barcodes/bc_data_{}.csv'.format(ver)
    df = pd.read_csv(fname)
    df.loc[df.well.duplicated(keep=False)].sort_values(by='well')

    if bc == 2 or bc == 3:
        assert len(df['type'].unique().tolist()) == 1

    drop_cols = ['bci', 'uid', 'type']
    df.drop(drop_cols, axis=1, inplace=True)
    df.rename({'sequence': bc_name}, axis=1, inplace=True)

    return df
