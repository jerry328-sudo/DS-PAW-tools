'''## functions
class Read_As_File:
def dos_read
def coordinate_transform
'''

import json
import os
import numpy as np

# read the as file class

class Read_As_File:
    '''- a function to deal with the as file
    - Read_As_File(file = 'structure.as', path = os.getcwd())
    - target: the line number of the Cartesian coordinates
    - n_atoms: the number of atoms
    - atom_type: the type of atoms
    - I: the lines of the as file
    - Ic: the lines of the as file after split
    '''
    def __init__(self, file = 'structure.as', path = os.getcwd()):
        self.file = path+'\\'+file
        self.Ic = []
        self.I = []
        self.n_atoms = 0
        self.target = 6
        self.atom_type = {}
        self.__read_as_file()

    def __read_as_file(self):
        '''- read the as file'''
        print(self.file)
        with open(self.file, 'r') as f:
            self.I = f.readlines()
            for i in range(len(self.I)):
                self.Ic.append(self.I[i].split())

        # count the number of atoms
        self.n_atoms = int(self.Ic[1][0])
        self.target = len(self.I)
        for i in range(len(self.I)):
            if 'Cartesian' in self.Ic[i]:
                self.target = i
            if i > self.target:
                if self.Ic[i][0] not in self.atom_type.keys():
                    self.atom_type[self.Ic[i][0]] = 1
                else:
                    self.atom_type[self.Ic[i][0]] += 1
    
    def write_as_file(self, path = os.getcwd(), file = 'structure.as'):
        '''- write the as file'''
        with open(path + '\\' + file, 'w') as f:
            for i in range(len(self.Ic)):
                f.write(' '.join(self.Ic[i])+'\n')


def dos_read( path, file = 'dos.json'):
    '''数据结构：
    dosinfo{
        'DosEnergy': dos['DosInfo']['DosEnergy'],
        'EFermi': dos['DosInfo']['EFermi'],
        'NumberOfDos': dos['DosInfo']['NumberOfDos'],
        'Lattice': dos['AtomInfo']['Lattice'],
        'Spin1': dos['DosInfo']['Spin1']['Dos'],
        'Spin2': dos['DosInfo']['Spin2']['Dos'],
        'cellinfo':{
            {
            'Element': 'C',
            'Position': [0.1946, 0.6121, 0.4334],
            'spin1': {'py': ......, 'pz': ......, 'px': ......,......},
            'spin1': {'py': ......, 'pz': ......, 'px': ......,......},
            }
        }
    }
    '''

    filepath = path + '\\' + file
    # load dos.json
    with open( filepath , 'r') as f:
        dos = json.load(f)

    cellinfo = dos['AtomInfo']['Atoms']


    spin1 = [0]*9
    spin2 = [0]*9
    for i in range(len(cellinfo)):
        temp1 = 0
        temp2 = 0
        for i1 in range(len(dos['DosInfo']['Spin1']['ProjectDos'])):
            if i+1 == dos['DosInfo']['Spin1']['ProjectDos'][i1]['AtomIndex']:
                spin1[temp1]=dos['DosInfo']['Spin1']['ProjectDos'][i1]['Contribution']
                temp1 += 1
        if "Spin2" in dos['DosInfo']:
            for i1 in range(len(dos['DosInfo']['Spin2']['ProjectDos'])):
                if i+1 == dos['DosInfo']['Spin2']['ProjectDos'][i1]['AtomIndex']:
                    spin2[temp2]=dos['DosInfo']['Spin2']['ProjectDos'][i1]['Contribution']
                    temp2 += 1

        # 写入自旋轨道态密度
        cellinfo[i]['spin1']={dos['DosInfo']['Orbit'][0]:spin1[0]}
        if "Spin2" in dos['DosInfo']:
                cellinfo[i]['spin2']={dos['DosInfo']['Orbit'][0]:spin2[0]}
        for j in range(9):
            cellinfo[i]['spin1'][dos['DosInfo']['Orbit'][j]]=spin1[j]
            if "Spin2" in dos['DosInfo']:
                cellinfo[i]['spin2'][dos['DosInfo']['Orbit'][j]]=spin2[j]
        # p轨道数据
        p = np.array(spin1[1]) + np.array(spin1[2]) + np.array(spin1[3])
        # d轨道数据
        d = np.array(spin1[4]) + np.array(spin1[5]) + np.array(spin1[6]) + np.array(spin1[7]) + np.array(spin1[8])
        cellinfo[i]['spin1']['p']=p.tolist()
        cellinfo[i]['spin1']['d']=d.tolist()
        if "Spin2" in dos['DosInfo']:
            p = np.array(spin2[1]) + np.array(spin2[2]) + np.array(spin2[3])
            d = np.array(spin2[4]) + np.array(spin2[5]) + np.array(spin2[6]) + np.array(spin2[7]) + np.array(spin2[8])
            cellinfo[i]['spin2']['p']=p.tolist()
            cellinfo[i]['spin2']['d']=d.tolist()

    dosinfo = {
        'DosEnergy': dos['DosInfo']['DosEnergy'],
        'EFermi': dos['DosInfo']['EFermi'],
        'NumberOfDos': dos['DosInfo']['NumberOfDos'],
        'Lattice': dos['AtomInfo']['Lattice'],
        'Spin1': dos['DosInfo']['Spin1']['Dos'],
    }
    if "Spin2" in dos['DosInfo']:
        dosinfo['Spin2'] = dos['DosInfo']['Spin2']['Dos']
    dosinfo['cellinfo'] = cellinfo
    return dosinfo

def coordinate_transform(dos,coor):
    '''将坐标从晶胞坐标系转换为直角坐标系'''
    za = np.array(dos['Lattice'][0:3])
    zb = np.array(dos['Lattice'][3:6])
    zc = np.array(dos['Lattice'][6:9])
    ax = za[0]*coor[0] + zb[0]*coor[1] + zc[0]*coor[2]
    ay = za[1]*coor[0] + zb[1]*coor[1] + zc[1]*coor[2]
    az = za[2]*coor[0] + zb[2]*coor[1] + zc[2]*coor[2]
    return [ax,ay,az]