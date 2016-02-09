
import cPickle
kmeds = cPickle.load('kmeds.pickle')

import mdtraj as md

Y = [y.atom_slice(heavy_atoms) for y in md.load(filenames)]
dtrajs = kmeds.transform(Y)


