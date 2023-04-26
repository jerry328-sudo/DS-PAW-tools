# %%
from ase import Atoms, Atom
from ase.visualize import view
import re
import ase.units as units
import h5py
import numpy as np

def molar_mass(element):
    """
    Calculation of molar mass
    - param element: Element
    - return: Molar mass
    """
    mass_list = {'H': 1.00794, 'He': 4.002602, 'Li': 6.941, 'Be': 9.012182, 'B': 10.811, 
                'C': 12.0107, 'N': 14.0067, 'O': 15.9994, 'F': 18.9984032, 'Ne': 20.1797, 
                'Na': 22.98976928, 'Mg': 24.3050, 'Al': 26.9815386, 'Si': 28.0855, 
                'P': 30.973762, 'S': 32.065, 'Cl': 35.453, 'Ar': 39.948, 'K': 39.0983, 
                'Ca': 40.078, 'Sc': 44.955912, 'Ti': 47.867, 'V': 50.9415, 'Cr': 51.9961, 
                'Mn': 54.938045, 'Fe': 55.845, 'Co': 58.933195, 'Ni': 58.6934, 'Cu': 63.546, 
                'Zn': 65.38, 'Ga': 69.723, 'Ge': 72.64, 'As': 74.92160, 'Se': 78.96, 
                'Br': 79.904, 'Kr': 83.798, 'Rb': 85.4678, 'Sr': 87.62, 'Y': 88.90585, 
                'Zr': 91.224, 'Nb': 92.90638, 'Mo': 95.96, 'Tc': 98, 'Ru': 101.07, 
                'Rh': 102.90550, 'Pd': 106.42, 'Ag': 107.8682, 'Cd': 112.411, 'In': 114.818, 
                'Sn': 118.710, 'Sb': 121.760, 'Te': 127.60, 'I': 126.90447, 'Xe': 131.293, 
                'Cs': 132.9054519, 'Ba': 137.327, 'La': 138.90547, 'Ce': 140.116, 
                'Pr': 140.90765, 'Nd': 144.242, 'Pm': 145, 'Sm': 150.36, 'Eu': 151.964,
                'Gd': 157.25, 'Tb': 158.92535, 'Dy': 162.500, 'Ho': 164.93032, 'Er': 167.259,
                'Tm': 168.93421, 'Yb': 173.054, 'Lu': 174.9668, 'Hf': 178.49, 'Ta': 180.94788,
                'W': 183.84, 'Re': 186.207, 'Os': 190.23, 'Ir': 192.217, 'Pt': 195.084,
                'Au': 196.966569, 'Hg': 200.59, 'Tl': 204.3833, 'Pb': 207.2, 'Bi': 208.98040,
                'Po': 209, 'At': 210, 'Rn': 222, 'Fr': 223, 'Ra': 226, 'Ac': 227, 'Th': 232.03806,
                'Pa': 231.03588, 'U': 238.02891, 'Np': 237, 'Pu': 244, 'Am': 243, 'Cm': 247,
                'Bk': 247, 'Cf': 251, 'Es': 252, 'Fm': 257, 'Md': 258, 'No': 259, 'Lr': 262,
                'Rf': 261, 'Db': 262, 'Sg': 266, 'Bh': 264, 'Hs': 277, 'Mt': 268, 'Ds': 281,
                'Rg': 272, 'Cn': 285, 'Nh': 284, 'Fl': 289, 'Mc': 288, 'Lv': 293, 'Ts': 294,
                'Og': 294}
    return mass_list[element]

def write_structure(structure, file="structure.as"):
    '''将结构写入structure.as文件, structure为ase.Atoms类型, 默认所有原子全部放开'''
    with open(file, "w") as f:
        f.write("Total number of atoms\n")
        f.write(str(structure.get_global_number_of_atoms()) + "\n")
        f.write("Lattice\n")
        for i in range(3):
            f.write(" "+" ".join([str(num)
                    for num in structure.get_cell()[i]]) + "\n")
        f.write("Cartesian Mag Fix_x Fix_y Fix_z\n")
        for i in range(structure.get_global_number_of_atoms()):
            f.write(" ".join([structure[i].symbol, " ".join(
                [str(num) for num in structure[i].position.tolist()]), str(structure[i].magmom), "F F F"]) + "\n")


def read_hzw(file = None, magmoms = None):
    '''读取hzw文件, 返回ase.Atoms类型, magmoms为字典类型, key为元素符号, value为对应元素的玻尔磁子数'''
    # %%
    with open(file, "r") as f:
        lines = f.readlines()
        lines_split = [line.split() for line in lines]

    # %%
    # 设置参数

    # %%
    structure = Atoms(pbc=True)
    atom_line_index = 0
    uni_cell_line_index = 0
    atom_number_line_index = 0
    cell = []

    for i in range(len(lines)):
        # 确定各个参数的行数
        if re.match("% Uni-cell vector", lines[i]):
            uni_cell_line_index = i
        if re.match("%Total number of device_structure", lines[i]):
            atom_number_line_index = i
        if re.match("%Atom site", lines[i]):
            atom_line_index = i
            break

    atom_index = 0
    for i in range(len(lines)):
        # 读取晶胞信息
        if i > uni_cell_line_index and i < uni_cell_line_index + 4:
            cell.append([float(line_split) for line_split in lines_split[i] if line_split != ""])
        if i == atom_number_line_index: # 将cell写入structure
            structure.set_cell(cell)

        # 读取原子信息
        if i > atom_line_index:
            if len(lines_split[i]) == 0:
                break
            position = lines_split[i][1:4]
            position = [float(num) for num in position]
            symbol = lines_split[i][0]
            if symbol in magmoms.keys():
                atom_temp = Atom(symbol=lines_split[i][0], position=position, magmom=magmoms[symbol], index=atom_index)
            else:
                atom_temp = Atom(symbol=lines_split[i][0], position=position, magmom=0, index=atom_index)
            structure.append(atom_temp)
            atom_index += 1
    return structure


class freLoad:
    """def __init__(self, filename = "frequency.h5"):
    filename;
    atoms;
    fre_length;
    Frequency;
    FrequencyType;
    Unit;
    EigenVectors"""
    filename = None
    atoms = None
    fre_length = None
    Frequency = []
    FrequencyType = []
    Unit = []
    EigenVectors = []
    def __init__(self, filename = "frequency.h5"):
        self.filename = filename
        self.atoms = self.__structureLoad()
        self.__frequencyLoad()
    def __structureLoad(self):
        data = h5py.File('frequency.h5', 'r')
        CoordinateType = [x.decode('UTF-8') for x in data["AtomInfo"]["CoordinateType"]]
        CoordinateType = "".join(CoordinateType)
        if CoordinateType == "Cartesian":
            pass
        else:
            raise ValueError("CoordinateType is not Cartesian")
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
        # 构造原子类
        atoms = Atoms(cell = Lattice, pbc = True)
        for i in range(len(Elements)):
            atom = Atom(symbol = Elements[i], 
                        position = Position[i], 
                        mass = molar_mass(Elements[i]),
                        magmom = Mag[i], tag = Fix[i])
            atoms.append(atom)
        return atoms
    def __frequencyLoad(self):
        data = h5py.File('frequency.h5', 'r')
        self.fre_length = len(data["FrequencyInfo"]["EigenValues"]["Frequency"])
        for i in range(1, self.fre_length+1):
            # 频率值
            value = data["FrequencyInfo"]["EigenValues"]["Frequency"][str(i)][0]
            self.Frequency.append(value)
            # 频率类型
            value = [x.decode('UTF-8') for x in data["FrequencyInfo"]["EigenValues"]["FrequencyTpye"][str(i)]]
            value = "".join(value)
            self.FrequencyType.append(value)
            # 频率单位
            value = [x.decode('UTF-8') for x in data["FrequencyInfo"]["EigenValues"]["Unit"][str(i)]]
            value = "".join(value)
            self.Unit.append(value)
            # 频率对应的本征矢量
            value = [x.astype(np.float64) for x in data["FrequencyInfo"]["EigenVectors"][str(i)]]
            self.EigenVectors.append(value)