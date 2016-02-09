

# soo, I accidentally scrambled the orders of the entries of the tica.npz and dtrajs.npz archives
# here I'm remedying that!

# what happened was something like this:
#
#X = [np.ones(50)*i for i in range(100)]
#np.savez_compressed('test.npz',*X)
#X_ = [x[1] for x in np.load('test.npz').items()]
#np.savez_compressed('test2.npz',*X_)
#Y = [x[1] for x in np.load('test2.npz').items()]
#
# and afterwards I (incorrectly) assumed X and Y would be in the same order

import numpy as np

def find_one_step_mapping(filename):
    ''' if I read the arrays in a compressed archive in the wrong order, using
    
    scrambled_list = [x[1] for x in np.load('test.npz').items()]
    
    what mapping from indices to indices
    
    will make sorted_list[i] = scrambled_list[mapping[i]] ?
    '''
    
    keys = np.load(filename).keys()
    sorted_keys = sorted(keys,key=lambda i:int(i[i.index('_')+1:]))
    
    mapping = dict()
    
    for i in range(len(keys)):
        index = keys.index(sorted_keys[i])
        mapping[i] = index
    
    return mapping



def load_dtrajs_from_npz(filename):
    scrambled_items = np.load(filename).items()
    #dtrajs = [y[1] for y in np.load(filename).items()]
    
    # since in both collections the number and names of the arrays were the same, and the dictionary
    # hash applied to the string names of the arrays, 
    # I only need to know one mapping and apply it twice
    mapping = find_one_step_mapping(filename)
    dtrajs = [scrambled_items[mapping[mapping[i]]][1] for i in range(len(scrambled_items))]

    return dtrajs
