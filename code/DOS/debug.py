import json
import os
import numpy as np
# import scipy

# load dos.json
with open('dos.json', 'r') as f:
    dos = json.load(f)

cellinfo = dos['AtomInfo']['Atoms']

temp1 = 0
temp2 = 0
spin1 = [0]*9
spin2 = []
for i in range(len(cellinfo)):
    temp1 = 0
    for i1 in range(len(dos['DosInfo']['Spin1']['ProjectDos'])):
        if i+1 == dos['DosInfo']['Spin1']['ProjectDos'][i1]['AtomIndex']:
            spin1[temp1] = dos['DosInfo']['Spin1']['ProjectDos'][i1]['Contribution']
            temp1 += 1
            if np.sum(dos['DosInfo']['Spin1']['ProjectDos'][i1]['Contribution']) == 0:
                print(cellinfo[i]['Element'], i+1)
    if "Spin2" in dos['DosInfo']:
        for i1 in range(len(dos['DosInfo']['Spin2']['ProjectDos'])):
            if i+1 == dos['DosInfo']['Spin2']['ProjectDos'][i1]['AtomIndex']:
                spin2.append(dos['DosInfo']['Spin2']
                             ['ProjectDos'][i1]['Contribution'])

for i in range(9):
    print(i+1, np.sum(spin1[i]))