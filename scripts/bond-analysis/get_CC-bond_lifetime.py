#!/usr/bin/python3
"""
python3 gen_bond_length.py --fxyz *.xyz
--rcut $rcut
"""

import argparse
import matplotlib.pyplot as plt
import sys
import numpy as np
from ase.io import read, write
from ase import Atoms
from ase.neighborlist import NeighborList
from scipy.spatial.distance import cdist

from asaplib.io import str2bool

def read_frame(filedesc):
    # 1-3
    comment = filedesc.readline()
    comment = filedesc.readline()
    comment = filedesc.readline()
    # 4
    natoms = int(filedesc.readline())
    #print(natoms)
    # 5
    comment = filedesc.readline()
    # 6,7,8
    cell = np.zeros((3,3),float)
    # line 1
    size = filedesc.readline().split()
    cell[0,0] = abs(float(size[1])-float(size[0]))
    # line 2
    size = filedesc.readline().split()
    cell[1,1] = abs(float(size[1])-float(size[0]))
    # line 1
    size = filedesc.readline().split()
    cell[2,2] = abs(float(size[1])-float(size[0]))

    #print(cell)
    # 9
    comment = filedesc.readline()
    # ITEM: ATOMS type x y z
    names = np.zeros(natoms, int)
    q = np.zeros((natoms,3),float)
    for i in range(natoms):
        line = filedesc.readline().split()
        if line[0] == 'C':
            names[i] = 6
        else:
            names[i] = 1
        q[i,0] = float(line[1])
        q[i,1] = float(line[2])
        q[i,2] = float(line[3])
    #print(len(names))
    return [natoms, cell, names, q]

def main(fxyz, prefix, cutoff):

    #frames = read(fxyz, index=':')
    ixyz = open(fxyz,"r")
    frames = []
    while True:
        try:
            [ na, cell, names, pos] = read_frame(ixyz)
            if np.sum(cell) > 0:
                pbc = [True, True, True]
            else:
                pbc = [False,False,False]
            #print(names, cell, pos, pbc)
            frtemp = Atoms(numbers=names,cell=cell,positions=pos,pbc=pbc)
            frames.append(frtemp)
            #print(frtemp)
        except:
            break

    # number of the carbon atoms has to be the same
    # and the order needs to be the same
    all_carbon_bonds = {}
    c_index = np.where(frames[0].get_atomic_numbers()==6)[0]
    print(c_index)
    for central_atom in c_index:
        all_carbon_bonds[central_atom] = {}

    # intialize the neighborlist pbject
    natoms = len(frames[0].get_positions())    
    r_cut_list =  np.zeros(natoms)
    r_cut_list[c_index] = cutoff
    #print(r_cut_list)
    nl = NeighborList(r_cut_list, skin=0., sorted=False, self_interaction=False,
                 bothways=True)

    # store the lifetime of each carbon bond
    lifetime_counter = []
    # drop 1/3 of the total frames before recording the bond life time
    frame_drop = len(frames)/3.

    for num_frame,frame in enumerate(frames):
        # build neighborlist
        nl.update(frame)

        for central_atom in c_index:
            nb_indices, _ = nl.get_neighbors(central_atom)

            remove_bond = []
            # check if the old bonds are still there
            for old_nb_indices in all_carbon_bonds[central_atom]:
                if old_nb_indices not in nb_indices:
                    remove_bond.append(old_nb_indices)
                    if num_frame > frame_drop:
                        lifetime_counter.append(all_carbon_bonds[central_atom][old_nb_indices])
            for old_nb_indices in remove_bond:
                del all_carbon_bonds[central_atom][old_nb_indices]

            # exiting carbon bond
            for nb_now in nb_indices:
                if nb_now in all_carbon_bonds[central_atom]:
                    all_carbon_bonds[central_atom][nb_now] += 1
                else: # new bond 
                    all_carbon_bonds[central_atom][nb_now] = 1

        print(int(100*num_frame/len(frames)),"%") 

    np.savetxt(str(prefix)+'-CC-bond-life-time.dat', lifetime_counter, fmt='%i')

    remaining_lifetime = []
    # also dump out the bond lifetime of the remaining bonds
    for central_atom in c_index:
        for old_nb_indices in all_carbon_bonds[central_atom]:
            remaining_lifetime.append(all_carbon_bonds[central_atom][old_nb_indices])
    np.savetxt(str(prefix)+'-CC-bond-remaining-life-time.dat', remaining_lifetime, fmt='%i')

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-fxyz', type=str, required=True, help='Location of lammpstrj file')
    parser.add_argument('--prefix', type=str, default='ASAP', help='Filename prefix')
    parser.add_argument('--rcut', type=float, default=0.8, help='Cutoff radius/2')
    args = parser.parse_args()

    main(args.fxyz, args.prefix, args.rcut)
