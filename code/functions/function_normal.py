'''## functions
class Read_As_File:
def dos_read
def coordinate_transform
'''

import json
import os
import numpy as np
import scipy.interpolate as si
import matplotlib.pyplot as plt
import h5py

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

    def __init__(self, file='structure.as', path=os.getcwd()):
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

    def write_as_file(self, path=os.getcwd(), file='structure.as'):
        '''- write the as file'''
        with open(path + '\\' + file, 'w') as f:
            for i in range(len(self.Ic)):
                f.write(' '.join(self.Ic[i])+'\n')

    def set_fix_atom(self, methode=0, direction=0, thick=0, fixed_element=[], freedom='T T T', filename=0, path=os.getcwd()):
        '''- set the fix atom
        - methode: 0: error
                   1: fix the atom by element
                   2: fix the atom by coordinate
        - direction: 0: x;
                     1: y;
                     2: z;
        - thick: the thickness of the fixed atom
        - fixed_element: the element of the fixed atom
        - freedom: the freedom of the fixed atom
        - filename: the name of the fixed atom file
        - path: the path of the fixed atom file
        '''
        count_free_atom = 0
        count_fix_atom = 0
        count_atom = 0
        if methode == 0:
            raise ValueError('The methode is not defined')
        elif methode == 1:
            if len(fixed_element) == 0:
                raise ValueError('The fixed element is not defined')
            else:
                for i in range(len(self.Ic)):
                    if 'Cartesian' in self.Ic[i]:
                        if 'Fix' not in self.I[i]:
                            self.Ic[i].append('Fix_x Fix_y Fix_z')
                            temp = 1
                        else:
                            temp = 0
                            pass
                    if i > self.target:
                        count_atom += 1
                        if self.Ic[i][0] in fixed_element:
                            if temp == 1:
                                self.Ic[i].append(freedom)
                                count_fix_atom += 1
                            elif temp == 0:
                                del self.Ic[i][-3:]
                                self.Ic[i].append(freedom)
                                count_fix_atom += 1
                        else:
                            if temp == 1:
                                self.Ic[i].append('F F F')
                                count_free_atom += 1
                            elif temp == 0:
                                count_free_atom += 1
                                pass
                if count_atom != self.n_atoms:
                    raise ValueError('The number of atoms is not correct!')
                if filename == '0':
                    self.write_as_file(path=path)
                else:
                    self.write_as_file(path=path, file=filename)
                with open(path + '\\' + 'atom_count.txt', 'w') as f:
                    f.write('Total number of atoms: {}\n'.format(count_atom))
                    f.write('Number of fixed atoms: {}\n'.format(count_fix_atom))
                    f.write('Number of free atoms: {}\n'.format(count_free_atom))
        elif methode == 2:
            if direction not in [0, 1, 2]:
                raise ValueError('The direction is not defined')
            for i in range(len(self.Ic)):
                if 'Cartesian' in self.Ic[i]:
                    if 'Fix' not in self.I[i]:
                        self.Ic[i].append('Fix_x Fix_y Fix_z')
                        temp = 1
                    else:
                        temp = 0
                        pass
                if i > self.target:
                    count_atom += 1
                    if float(self.Ic[i][direction+1]) < thick:
                        if temp == 1:
                            self.Ic[i].append(freedom)
                            count_fix_atom += 1
                        elif temp == 0:
                            del self.Ic[i][-3:]
                            self.Ic[i].append(freedom)
                            count_fix_atom += 1
                    else:
                        if temp == 1:
                            self.Ic[i].append('F F F')
                            count_free_atom += 1
                        elif temp == 0:
                            count_free_atom += 1
                            pass
            if count_atom != self.n_atoms:
                raise ValueError('The number of atoms is not correct!')
            if filename == '0':
                self.write_as_file(path=path)
            else:
                self.write_as_file(path=path, file=filename)
            with open(path + '\\' + 'atom_count.txt', 'w') as f:
                f.write('Total number of atoms: {}\n'.format(count_atom))
                f.write('Number of fixed atoms: {}\n'.format(count_fix_atom))
                f.write('Number of free atoms: {}\n'.format(count_free_atom))
        else:
            raise ValueError('The methode is not defined')




def dos_read(path, file='dos.json'):
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
    with open(filepath, 'r') as f:
        dos = json.load(f)

    cellinfo = dos['AtomInfo']['Atoms']

    spin1 = [0]*9
    spin2 = [0]*9
    for i in range(len(cellinfo)):
        temp1 = 0
        temp2 = 0
        for i1 in range(len(dos['DosInfo']['Spin1']['ProjectDos'])):
            if i+1 == dos['DosInfo']['Spin1']['ProjectDos'][i1]['AtomIndex']:
                spin1[temp1] = dos['DosInfo']['Spin1']['ProjectDos'][i1]['Contribution']
                temp1 += 1
        if "Spin2" in dos['DosInfo']:
            for i1 in range(len(dos['DosInfo']['Spin2']['ProjectDos'])):
                if i+1 == dos['DosInfo']['Spin2']['ProjectDos'][i1]['AtomIndex']:
                    spin2[temp2] = dos['DosInfo']['Spin2']['ProjectDos'][i1]['Contribution']
                    temp2 += 1

        # 写入自旋轨道态密度
        cellinfo[i]['spin1'] = {dos['DosInfo']['Orbit'][0]: spin1[0]}
        if "Spin2" in dos['DosInfo']:
            cellinfo[i]['spin2'] = {dos['DosInfo']['Orbit'][0]: spin2[0]}
        for j in range(9):
            cellinfo[i]['spin1'][dos['DosInfo']['Orbit'][j]] = spin1[j]
            if "Spin2" in dos['DosInfo']:
                cellinfo[i]['spin2'][dos['DosInfo']['Orbit'][j]] = spin2[j]
        # p轨道数据
        p = np.array(spin1[1]) + np.array(spin1[2]) + np.array(spin1[3])
        # d轨道数据
        d = np.array(spin1[4]) + np.array(spin1[5]) + \
            np.array(spin1[6]) + np.array(spin1[7]) + np.array(spin1[8])
        cellinfo[i]['spin1']['p'] = p.tolist()
        cellinfo[i]['spin1']['d'] = d.tolist()
        if "Spin2" in dos['DosInfo']:
            p = np.array(spin2[1]) + np.array(spin2[2]) + np.array(spin2[3])
            d = np.array(spin2[4]) + np.array(spin2[5]) + \
                np.array(spin2[6]) + np.array(spin2[7]) + np.array(spin2[8])
            cellinfo[i]['spin2']['p'] = p.tolist()
            cellinfo[i]['spin2']['d'] = d.tolist()

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


def coordinate_transform(dos, coor):
    '''将坐标从晶胞坐标系转换为直角坐标系'''
    za = np.array(dos['Lattice'][0:3])
    zb = np.array(dos['Lattice'][3:6])
    zc = np.array(dos['Lattice'][6:9])
    ax = za[0]*coor[0] + zb[0]*coor[1] + zc[0]*coor[2]
    ay = za[1]*coor[0] + zb[1]*coor[1] + zc[1]*coor[2]
    az = za[2]*coor[0] + zb[2]*coor[1] + zc[2]*coor[2]
    return [ax, ay, az]

def plot_neb_barrier(neb_h5, interp_kinds="cubic", data_clean=False, data_save=False):
    '''给出插值的样条曲线的阶数
    'zero' 、'nearest'零阶
    'slinear' 、'linear'线性
    'quadratic' 、'cubic'二阶和三阶样条曲线,更高阶的曲线可以直接使用整数值指定'''
    neb = h5py.File(neb_h5, 'r')
    data_form = ' %19.6f'
    reaction_coordinate = neb['BarrierInfo']['ReactionCoordinate'][:]
    energy = neb['BarrierInfo']['TotalEnergy'][:]

    x = []
    for c in reaction_coordinate:
        if len(x) > 0:
            x.append(x[-1] + c)
        else:
            x.append(c)

    y = [x-energy[0] for x in energy]

    # 清理数据，挑选最大值
    maxE = 0
    for i in range(len(y)):
        if y[i] > y[maxE]:
            maxE = i
    if data_clean:
        x = [0, 1, 2]
        # x = [x[0], x[maxE], x[-1]]
        y = [y[0], y[maxE], y[-1]]

    inter_f = si.interp1d(x,y,kind=interp_kinds)
    xnew = np.linspace(x[0],x[-1],100)
    ynew = inter_f(xnew)
    if data_save:
        with open("neb_data.txt","w") as file:
            file.write('raw data'+'\n')
            for i in range(len(x)):
                file.write(data_form % x[i] + '\t'+ data_form % y[i] + '\n')
            file.write('interpolation data'+'\n')
            for i in range(len(xnew)):
                file.write(data_form % xnew[i] + " "+ data_form % ynew[i]+"\n")

    # 样条插值
    # tck = si.splrep(x, y, s=0)
    # xnew = np.arange(x[0], x[-1], (x[-1]-x[0])/100)
    # ynew = si.splev(xnew, tck, der=0)

    plt.plot(xnew,ynew,c="b")
    plt.scatter(x,y,c="r")
    plt.xlabel("Reaction Coordinate")
    plt.ylabel("Energy")