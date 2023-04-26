from ase import Atoms, Atom
from ase.io.vasp import write_vasp
import h5py
import numpy as np
import os
import sys
current_path = os.path.dirname(os.path.abspath(__file__))
parent_path = os.path.dirname(current_path)
sys.path.append(parent_path)
import functions.aseFunction as af

def tran(fix):
    if fix == 0:
        return 'F F F'
    elif fix == 1:
        return 'T T T'

# get the directory where the py file is located
# dir = os.path.dirname(os.path.realpath(__file__))
# os.chdir(dir)
dir = os.getcwd()

# load DS-PAW file
data = h5py.File('frequency.h5', 'r')
CoordinateType = [x.decode('UTF-8') for x in data["AtomInfo"]["CoordinateType"]]
CoordinateType = "".join(CoordinateType)
if CoordinateType == "Cartesian":
    pass
elif CoordinateType == "Direct":
    exit("Please use the Cartesian coordinate system")

Elements = [x.decode('UTF-8') for x in data["AtomInfo"]["Elements"]]
temp1 = [] # 将Elements中的符号转换为正常的格式
temp2 = []
for letter in Elements:
    if letter != ";":
        temp2.append(letter)
    elif letter == ";":
        temp1.append("".join(temp2))
        temp2 = []
Elements = temp1
# 读取结构信息
Fix = [x.astype(np.float64) for x in data["AtomInfo"]["Fix"]]
Mag = [x.astype(np.float64) for x in data["AtomInfo"]["Mag"]]
Position = [x.astype(np.float64) for x in data["AtomInfo"]["Position"]]
Position = np.array(Position)
Position = Position.reshape(-1, 3)
Lattice = [x.astype(np.float64) for x in data["AtomInfo"]["Lattice"]]
Lattice = np.array(Lattice)
Lattice = Lattice.reshape(-1, 3)
# 读取频率信息
fre_length = len(data["FrequencyInfo"]["EigenValues"]["Frequency"])
Frequency = []
FrequencyType = []
Unit = []
EigenVectors = []
for i in range(1, fre_length+1):
    # 频率值
    value = data["FrequencyInfo"]["EigenValues"]["Frequency"][str(i)][0]
    Frequency.append(value)
    # 频率类型
    value = [x.decode('UTF-8') for x in data["FrequencyInfo"]["EigenValues"]["FrequencyTpye"][str(i)]]
    value = "".join(value)
    FrequencyType.append(value)
    # 频率单位
    value = [x.decode('UTF-8') for x in data["FrequencyInfo"]["EigenValues"]["Unit"][str(i)]]
    value = "".join(value)
    Unit.append(value)
    # 频率对应的本征矢量
    value = [x.astype(np.float64) for x in data["FrequencyInfo"]["EigenVectors"][str(i)]]
    EigenVectors.append(value)
# 构造原子类
atoms = Atoms(cell = Lattice, pbc = True)
for i in range(len(Elements)):
    atom = Atom(symbol = Elements[i], 
                position = Position[i], 
                magmom = Mag[i], tag = Fix[i])
    atoms.append(atom)



# Determine whether the file folder exists
if not os.path.exists('vaspfiles'):

    # make a new directory in the current directory
    os.mkdir(dir + '/vaspfiles')

    # set the work path as the new directory
    os.chdir(dir + '/vaspfiles')
else:
    os.chdir(dir + '/vaspfiles')

frelen = len(Frequency)
Frequency = np.array(Frequency)
# make a new file
with open('OUTCAR', 'w') as f:
    f.write(''' Eigenvectors and eigenvalues of the dynamical matrix
 ----------------------------------------------------
 

''')
    for i in range(frelen):
        # 零点能单位换算：1(THz) = 6.28(2PiTHz) = 33.356 cm-1 = 4.135meV
        a = Frequency[i]
        p = np.around(a * 2 * np.pi, decimals=6)
        cm = np.around(a * 33.356, decimals=6)
        mev = np.around(a * 4.135, decimals=6)
        a = np.around(a, decimals=6)

        if FrequencyType[i] == 'f':
            f.write('   '+str(i+1)+' '+ FrequencyType[i] +'  '+'='+'%+12s'
                %str(a.tolist())+' '+Unit[i]+'%+13s'%str(p.tolist())
                +' 2PiTHz'+'%+13s'%str(cm.tolist())+' cm-1'+'%+13s'%str(mev.tolist())+' meV'+'\n'+'\n')
        else:
            f.write('   '+str(i+1)+' '+FrequencyType[i]+'='+'%+12s'
                %str(a.tolist())+' '+Unit[i]+'%+13s'%str(p.tolist())
                +' 2PiTHz'+'%+13s'%str(cm.tolist())+' cm-1'+'%+13s'%str(mev.tolist())+' meV'+'\n'+'\n')
        f.write('             X         Y         Z           dx          dy          dz'+'\n')
        for j in range(len(Elements)):
            f.write("      " + ('%.09f'%Position[j][0])[:11] + "  " + ('%.09f'%Position[j][1])[:11] + "  "
                    + ('%.09f'%Position[j][2])[:11] + "    " + ('%.09f'%EigenVectors[i][j][0])[:11] + "  "
                    + ('%.09f'%EigenVectors[i][j][1])[:11] + "  " + ('%.09f'%EigenVectors[i][j][2])[:11] + "\n")

write_vasp("CONTCAR", atoms)