from rdkit import Chem
from rdkit.Chem.Draw import IPythonConsole
import os
import numpy as np
import matplotlib.pyplot as plt
import json
import pickle
import collections
from collections import OrderedDict
from matplotlib.pyplot import cm


CHARPROTSET = { "A": 1, "C": 2, "B": 3, "E": 4, "D": 5, "G": 6, 
                "F": 7, "I": 8, "H": 9, "K": 10, "M": 11, "L": 12, 
                "O": 13, "N": 14, "Q": 15, "P": 16, "S": 17, "R": 18, 
                "U": 19, "T": 20, "W": 21, 
                "V": 22, "Y": 23, "X": 24, 
                "Z": 25 }

CHARPROTLEN = 25

CHARCANSMISET = { "#": 1, "%": 2, ")": 3, "(": 4, "+": 5, "-": 6, 
             ".": 7, "1": 8, "0": 9, "3": 10, "2": 11, "5": 12, 
             "4": 13, "7": 14, "6": 15, "9": 16, "8": 17, "=": 18, 
             "A": 19, "C": 20, "B": 21, "E": 22, "D": 23, "G": 24,
             "F": 25, "I": 26, "H": 27, "K": 28, "M": 29, "L": 30, 
             "O": 31, "N": 32, "P": 33, "S": 34, "R": 35, "U": 36, 
             "T": 37, "W": 38, "V": 39, "Y": 40, "[": 41, "Z": 42, 
             "]": 43, "_": 44, "a": 45, "c": 46, "b": 47, "e": 48, 
             "d": 49, "g": 50, "f": 51, "i": 52, "h": 53, "m": 54, 
             "l": 55, "o": 56, "n": 57, "s": 58, "r": 59, "u": 60,
             "t": 61, "y": 62}

CHARCANSMILEN = 62

CHARISOSMISET = {"#": 29, "%": 30, ")": 31, "(": 1, "+": 32, "-": 33, "/": 34, ".": 2, 
                "1": 35, "0": 3, "3": 36, "2": 4, "5": 37, "4": 5, "7": 38, "6": 6, 
                "9": 39, "8": 7, "=": 40, "A": 41, "@": 8, "C": 42, "B": 9, "E": 43, 
                "D": 10, "G": 44, "F": 11, "I": 45, "H": 12, "K": 46, "M": 47, "L": 13, 
                "O": 48, "N": 14, "P": 15, "S": 49, "R": 16, "U": 50, "T": 17, "W": 51, 
                "V": 18, "Y": 52, "[": 53, "Z": 19, "]": 54, "\\": 20, "a": 55, "c": 56, 
                "b": 21, "e": 57, "d": 22, "g": 58, "f": 23, "i": 59, "h": 24, "m": 60, 
                "l": 25, "o": 61, "n": 26, "s": 62, "r": 27, "u": 63, "t": 28, "y": 64}

CHARISOSMILEN = 64


## ######################## ##
#
#  Encoding Helpers
#
## ######################## ## 

#  Y = -(np.log10(Y/(math.pow(math.e,9))))

def one_hot_smiles(line, MAX_SMI_LEN, smi_ch_ind):
    X = np.zeros((MAX_SMI_LEN, len(smi_ch_ind))) #+1

    for i, ch in enumerate(line[:MAX_SMI_LEN]):
        X[i, (smi_ch_ind[ch]-1)] = 1 

    return X #.tolist()

def one_hot_sequence(line, MAX_SEQ_LEN, smi_ch_ind):
    X = np.zeros((MAX_SEQ_LEN, len(smi_ch_ind))) 
    for i, ch in enumerate(line[:MAX_SEQ_LEN]):
        X[i, (smi_ch_ind[ch])-1] = 1

    return X #.tolist()


def label_smiles(line, MAX_SMI_LEN, smi_ch_ind):
    X = np.zeros(MAX_SMI_LEN)
    for i, ch in enumerate(line[:MAX_SMI_LEN]): #    x, smi_ch_ind, y
        X[i] = smi_ch_ind[ch]

    return X #.tolist()

def label_sequence(line, MAX_SEQ_LEN, smi_ch_ind):
    X = np.zeros(MAX_SEQ_LEN)

    for i, ch in enumerate(line[:MAX_SEQ_LEN]):
        X[i] = smi_ch_ind[ch]

    return X #.tolist()




# Set the following ones properly
f = open('d:\\Y_PDB.txt', 'r')    # file which contain label of PDB data  (uploaded in git repository)
f1 = open('d:\\name_PDB.txt', 'r')  # file which contain the name of samples in pdb dataset  (uploaded in git repository)
path_pdb = 'C:\\Users\\abbasi\\Desktop\\pdbbind_v2016.tar\\pdbbind_v2016\\v2016\\'  # the path of dataset


XD = []  # matrix which stores drugs
XT = []  # matrix which stores proteins

y_value = [] 
y_unit = []
y_real = []
flag = []



for i in range(0, len(f)):
    try:
        tmp = f.readline()
        value = tmp.split('>=')
        if len(value)<2:
            value = tmp.split('>')
            if len(value)<2:
               value = tmp.split('<=')
               if len(value)<2:
                  value = tmp.split('<')
                  if len(value)<2:
                     value = tmp.split('=')
                     if len(value)<2:
                        value = tmp.split('~')
                        if len(value)<2:
                           print(value)
        
        y_value.append(value[1][:-3])
        y_unit.append(value[0])
        if value[0] == 'IC50':
            #   flag.append(0)
            #if value[0] == 'Kd':
            #   flag.append(1)
            #if value[0] == 'Ki':
            #   flag.append(2)       
            print(i)
            
               
            folder = f1.readline()
            folder = folder.split('\n')[0]
            m = Chem.MolFromMolFile(path_pdb+folder+'\\'+folder+'_ligand.sdf')
            #print(Chem.MolToSmiles(m))
            #print('here1')
            
            #print('here2')

            from Bio.PDB import PDBParser
              # this is the module described below
            d3to1 = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
            'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N', 
            'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W', 
            'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M', 'HOH': 'B'}
            # read
            sloppyparser = PDBParser(
                PERMISSIVE=True#, structure_builder=xpdb.SloppyStructureBuilder()
            )
            structure = sloppyparser.get_structure("MD_system", path_pdb + folder+'/'+folder+'_protein.pdb')
            #print(structure.get_chains(structure.get_models()))
            struct = PDBParser(QUIET=True).get_structure('tmp', path_pdb + folder+'/'+folder+'_protein.pdb')
            #print(struct)
            count=0
            for model in structure:
                for chain in model:
                    if count==0:
                        seq = []
                        for residue in chain:
                            seq.append(d3to1[residue.resname])
                        #print(seq)
                        count=count+1
            
            XD.append(label_smiles(Chem.MolToSmiles(m), 100, CHARISOSMISET))
            XT.append(label_sequence(seq, 1000, CHARPROTSET))
            if value[1][-3:-1] == 'pM':
               y_real.append(float(value[1][:-3])* (0.001))
            if value[1][-3:-1] == 'uM':
               y_real.append(float(value[1][:-3])* (10^(3)))
            if value[1][-3:-1] == 'mM':
               y_real.append(float(value[1][:-3])* (10^6))
            if value[1][-3:-1] == 'fM':
               y_real.append(float(value[1][:-3])* (0.000001))
            if value[1][-3:-1] == 'nM':
               y_real.append(float(value[1][:-3]))
            
    except:
        continue
import numpy as np
from scipy.io import savemat, loadmat


print(type(y_value))
print(type(y_unit))
print(set(y_unit))
print(type(y_real))
print(np.sum(np.array(flag)==0))
print(np.sum(np.array(flag)==1))
print(np.sum(np.array(flag)==2))

print(len(y_real))
print(len(XD))
print(len(XT))

pdb={'XT': XT, 'XD': XD, 'Y': y_real}
savemat('pdb_xt.mat',pdb)
