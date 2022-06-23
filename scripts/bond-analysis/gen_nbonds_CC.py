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
    ncarbon = 0
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
    ncarbon = len(np.where(frames[0].get_atomic_numbers()==6)[0])
    #print(ncarbon)
    nframe = len(frames)
    print(nframe)
    carbon_bonds_all = np.zeros((ncarbon,nframe),int)
    for num_frame,frame in enumerate(frames):
        natoms = len(frame.get_positions())
        # the index of the carbon atoms    
        c_array = np.where(frame.get_atomic_numbers()==6)
        c_index = c_array[0]
        # build neighborlist
        r_cut_list =  np.zeros(natoms)
        r_cut_list[c_array] = cutoff
        #print(r_cut_list)
        nl = NeighborList(r_cut_list, skin=0., sorted=False, self_interaction=False,
                 bothways=True)    
        nl.update(frame)

        # compute displacements r_ij
        avg_rij = np.zeros(natoms)
        for c_i,central_atom in enumerate(c_index):
            indices, offsets = nl.get_neighbors(central_atom)
            carbon_bonds_all[c_i,num_frame] = len(indices)
            #displacements = []
            #for i, offset in zip(indices, offsets):
            #    displacements.append(np.linalg.norm(frame.positions[i] + np.dot(offset, frame.get_cell()) - frame.positions[central_atom] ))
            #print(displacements)
            #try:
            #    avg_rij[central_atom] = np.mean(displacements)
            #except:
            #    avg_rij[central_atom] = 0.0
        #frame.new_array('avg_bond', avg_rij)
        print(int(100*num_frame/nframe),"%") 

    #write(str(prefix) + ".xyz", frames)
    np.savetxt(str(prefix)+'-n-carbon-bonds.dat', carbon_bonds_all, fmt='%i')

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-fxyz', type=str, required=True, help='Location of lammpstrj file')
    parser.add_argument('--prefix', type=str, default='ASAP', help='Filename prefix')
    parser.add_argument('--rcut', type=float, default=0.8, help='Cutoff radius/2')
    args = parser.parse_args()

    main(args.fxyz, args.prefix, args.rcut)
