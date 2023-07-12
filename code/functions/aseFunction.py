#!/usr/bin/env python
"""
This module provides some useful functions for ASE.
Here is a list of functions:
    - molar_mass: Calculation of molar mass
    - write_structure(structure, file="structure.as")(strongly not recommended): Write structure to file
    - write_CHGCAR(structure, grid, chgInfo, filename = "CHGCAR")
    - write_as_file(path=os.getcwd(), filename='structure.as', structure = None)
    - read_hzw(file = None, magmoms = None)
    - read_as(file = "structure.as")
    - class freLoad
        - __init__(self, filename = "frequency.h5")
        - write_as_file(self, path=os.getcwd(), filename='structure.as')
    - class relaxLoad
        - __init__(self, filename = "relax.h5")
        - write_as_file(self, path=os.getcwd(), filename='structure_new.as', save_kind = 'final')
    - class pchargeLoad
        - __init__(self, filename = "pcharge.h5")
    - class scfLoad
        - __init__(self, filename = "scf.h5")
        - write_bands_occupation(self)
"""


from ase import Atoms, Atom
from ase.visualize import view
import re
import ase.units as units
import h5py
import numpy as np
import os
# np.set_printoptions(suppress=True)

def tranf2num(fix):
    s1 = 0
    s2 = 0
    s3 = 0
    if fix[0] == "T":
        s1 = 100
    if fix[1] == "T":
        s2 = 10
    if fix[2] == "T":
        s3 = 1
    return s1+s2+s3

def tran(fix):
    if fix == 0:
        return "F F F"
    elif fix == 1:
        return "F F T"
    elif fix == 10:
        return "F T F"
    elif fix == 11:
        return "F T T"
    elif fix == 100:
        return "T F F"
    elif fix == 101:
        return "T F T"
    elif fix == 110:
        return "T T F"
    elif fix == 111:
        return "T T T"

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

def write_CHGCAR(structure, grid, chgInfo, filename = "CHGCAR"):
    '''structure: Atoms类'''
    latt_form = ' %19.6f'
    cform = ' %19.16f'
    element_form = ' %5s'
    count = {}
    with open(filename, 'w') as f:
        f.write("".join(structure.get_chemical_formula(mode='hill')))
        f.write("\n")
        f.write("1.0\n")
        for i in range(3):
            f.write(latt_form % structure.cell[i][0])
            f.write(latt_form % structure.cell[i][1])
            f.write(latt_form % structure.cell[i][2])
            f.write("\n")
        for i in structure.get_chemical_symbols():
            count[i] = count.get(i, 0) + 1
        for i in count:
            f.write(element_form % i)
        f.write("\n")
        for i in count:
            f.write(element_form % count[i])
        f.write("\n")
        f.write("Cartesian\n")
        for i in range(len(structure)):
            f.write(cform % structure.positions[i][0])
            f.write(cform % structure.positions[i][1])
            f.write(cform % structure.positions[i][2])
            f.write("\n")
        f.write("\n")
        f.write('  ' + " ".join([str(x) + " " for x in grid]))
        f.write("\n")
        for i, num in enumerate(chgInfo):
            f.write(' %19.16e' % num)
            if (i + 1) % 5 == 0:
                f.write("\n")
            else:
                f.write(" ")

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

def read_as(file = "structure.as"):
    with open(file, "r") as f:
        lines = f.readlines()
        lines_split = [line.split() for line in lines]
    # read the lattice
    lattice = []
    for i in range(3,6):
        lattice.append([float(num) for num in lines_split[i][:3]])
    atoms = Atoms(cell = lattice, pbc = True)
    # read the atoms
    sign = 0
    if "Mag" in lines_split[6]:
        signm = 1
    else:
        signm = 0
    if "Fix_x" in lines_split[6]:
        signf = 10
    else:
        signf = 0
    sign = signm + signf
    for i in range(7, len(lines_split)):
        element = lines_split[i][0]
        position = [float(num) for num in lines_split[i][1:4]]
        if sign == 0:
            atom = Atom(symbol = element, position = position, mass = molar_mass(element))
        elif sign == 1:
            atom = Atom(symbol = element, position = position, mass = molar_mass(element), magmom = float(lines_split[i][4]))
        elif sign == 10:
            fix = tranf2num(lines_split[i][4:7])
            atom = Atom(symbol = element, position = position, mass = molar_mass(element), tag = fix)
        elif sign == 11:
            fix = tranf2num(lines_split[i][5:8])
            atom = Atom(symbol = element, position = position, mass = molar_mass(element), magmom = float(lines_split[i][4]), tag = fix)
        atoms.append(atom)
    return atoms

def write_as_file(path=os.getcwd(), filename='structure.as', structure = None):
    '''- write the as file
    - path: the path to save the file
    - file: the name of the file
    - structure: ASE Atoms object
    '''
    latt_form = ' %19.6f'
    cform = ' %19.16f'
    element_form = ' %2s'
    atoms = structure
    with open(path + '\\' + filename, 'w') as f:
        f.write("Total number of atoms" + "\n")
        f.write(str(len(atoms)) + "\n")
        f.write("Lattice" + "\n")
        for i in range(3):
            f.write(latt_form % atoms.cell[i][0] + " " + latt_form % atoms.cell[i][1] + " " + latt_form % atoms.cell[i][2] + "\n")
        f.write("Cartesian Mag")
        if atoms[0].tag != -1:
            f.write(" Fix_x Fix_y Fix_z" + "\n")
        else:
            f.write("\n")
        for i in range(len(atoms)):
            f.write(element_form % atoms[i].symbol + " ")
            f.write(cform % atoms[i].position[0] + " " + cform % atoms[i].position[1] + " " + cform % atoms[i].position[2] + " ")
            f.write(cform % atoms[i].magmom)
            if atoms[0].tag != -1:
                f.write(" " + tran(atoms[i].tag) + "\n")
            else:
                f.write("\n")

class freLoad:
    """def __init__(self, filename = "frequency.h5"):
    filename;
    atoms;
    fre_length;
    Frequency; np.array
    FrequencyType;
    Unit;
    EigenVectors; np.array"""
    filename = None
    atoms = None
    fre_length = None
    imfreLocate = None
    name='vib'
    Frequency = []
    FrequencyType = []
    Unit = []
    EigenVectors = []
    # position = []
    # elements = []
    def __init__(self, filename = "frequency.h5"):
        self.filename = filename
        self.atoms = self.__structureLoad()
        self.__frequencyLoad()
        for i in range(self.fre_length):
            if self.FrequencyType[i] == 'f/i':
                self.imFreLocate = i
                break
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
        temp1.append("".join(temp2))
        Elements = temp1
        # self.elements = Elements
        # 读取结构信息
        if "Fix" in data["AtomInfo"]:
            Fix = np.array([x.astype(np.float64) for x in data["AtomInfo"]["Fix"]])
            Fix = Fix.reshape(-1, 3)
            Fix = Fix[:, 0]*100 + Fix[:, 1] * 10 + Fix[:, 2]
        else:
            Fix = np.zeros(len(data["AtomInfo"]["Position"]), dtype=int) - 1 # 如果文件中没有Fix信息，则所有固定值设为-1
        if "Mag" in data["AtomInfo"]:
            Mag = [x.astype(np.float64) for x in data["AtomInfo"]["Mag"]]
        else:
            Mag = np.zeros(len(Elements), dtype=np.float64).tolist()
        Position = [x.astype(np.float64) for x in data["AtomInfo"]["Position"]]
        Position = np.array(Position)
        Position = Position.reshape(-1, 3)
        # self.position = Position
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
        self.Frequency = np.array(self.Frequency)
        self.EigenVectors = np.array(self.EigenVectors)

    def write_as_file(self, path=os.getcwd(), filename='structure.as'):
        '''- write the as file'''
        latt_form = ' %19.6f'
        cform = ' %19.16f'
        element_form = ' %2s'
        with open(path + '\\' + filename, 'w') as f:
            f.write("Total number of atoms" + "\n")
            f.write(str(len(self.atoms)) + "\n")
            f.write("Lattice" + "\n")
            for i in range(3):
                f.write(latt_form % self.atoms.cell[i][0] + " " + latt_form % self.atoms.cell[i][1] + " " + latt_form % self.atoms.cell[i][2] + "\n")
            f.write("Cartesian Mag")
            if self.atoms[0].tag != -1:
                f.write(" Fix_x Fix_y Fix_z" + "\n")
            else:
                f.write("\n")
            for i in range(len(self.atoms)):
                f.write(element_form % self.atoms[i].symbol + " ")
                f.write(cform % self.atoms[i].position[0] + " " + cform % self.atoms[i].position[1] + " " + cform % self.atoms[i].position[2] + " ")
                f.write(cform % self.atoms[i].magmom)
                if self.atoms[0].tag != -1:
                    f.write(" " + tran(self.atoms[i].tag) + "\n")
                else:
                    f.write("\n")
    
    def write_jmol(self, name = "vibration"):
        """Writes file for viewing of the modes with jmol."""
        self.name = name
        with open(self.name + '.xyz', 'w') as fd:
            self._write_jmol(fd)

    def _write_jmol(self, fd):
        symbols = self.atoms.get_chemical_symbols()
        freq = self.Frequency
        # Convert to cm^-1 if necessary
        # 暂时默认单位为THz吧，后面有需要再改
        # if self.Unit[0] == "THz":
        #     freq = np.around(freq * 33.35641, decimals=6) # 1 THz = 33.35641 cm^-1
        # else:
        #     raise ValueError("The unit of frequency is not THz.")
        for n in range(len(freq)):
            # 写入原子数量
            fd.write('%6d\n' % len(self.atoms))
            # 确定是否为虚频
            if self.FrequencyType[n] == "f/i":
                c = 'i'
            else:
                c = ' '
            # 写入频率信息
            fd.write('Mode #%d, f = %.4f%s ' % (n, float(freq[n]), c))
            fd.write('%s\n' % self.Unit[n])
            # # 写入强度信息
            # if self.ir:
            #     fd.write(', I = %.4f (D/Å)^2 amu^-1.\n' % self.intensities[n])
            # else:
            #     fd.write('.\n')
            # 写入原子信息
            mode = self.EigenVectors[n]
            for i, pos in enumerate(self.atoms.positions):
                fd.write('%2s %12.5f %12.5f %12.5f %12.5f %12.5f %12.5f\n' %
                         (symbols[i], pos[0], pos[1], pos[2],
                          mode[i, 0], mode[i, 1], mode[i, 2]))

class relaxLoad:
    RelaxTraj = [] # 弛豫轨迹
    Force = [] # 轨迹中每一帧的受力信息
    atoms_initial = None # 初始结构信息
    atoms_final = None # 最终结构信息
    electron_number = None # 电子数
    energy = {}
    def __init__(self, filename = "relax.h5"):
        self.filename = filename
        self.__structureLoad()
    def __structureLoad(self):
        # 读取初始原子信息
        data = h5py.File(self.filename, 'r')
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
        temp1.append("".join(temp2))
        Elements = temp1
        # self.elements = Elements
        # 结构固定信息
        if "Fix" in data["AtomInfo"]:
            Fix = np.array([x.astype(np.float64) for x in data["AtomInfo"]["Fix"]])
            Fix = Fix.reshape(-1, 3)
            Fix = Fix[:, 0]*100 + Fix[:, 1] * 10 + Fix[:, 2]
        else:
            Fix = np.zeros(len(data["AtomInfo"]["Position"]), dtype=int) - 1 # 如果文件中没有Fix信息，则所有固定值设为-1
        # 读取初始磁矩信息
        if "Mag" in data["AtomInfo"]:
            InitialMag = [x.astype(np.float64) for x in data["AtomInfo"]["Mag"]]
        else:
            InitialMag = np.zeros(len(data["AtomInfo"]["Position"]), dtype=int).tolist()
        Position = [x.astype(np.float64) for x in data["AtomInfo"]["Position"]]
        Position = np.array(Position)
        Position = Position.reshape(-1, 3)
        # self.position = Position
        Lattice = [x.astype(np.float64) for x in data["AtomInfo"]["Lattice"]]
        Lattice = np.array(Lattice)
        Lattice = Lattice.reshape(-1, 3)
        # 构造原子类
        self.atoms_initial = Atoms(cell = Lattice, pbc = True)
        for i in range(len(Elements)):
            atom = Atom(symbol = Elements[i], 
                        position = Position[i], 
                        mass = molar_mass(Elements[i]),
                        magmom = InitialMag[i], tag = Fix[i])
            self.atoms_initial.append(atom)
        check = 1
        if "MagInfo" in data:
            if "TotalMagOnAtom" in data["MagInfo"]:
                FinalMag = [x.astype(np.float64) for x in data["MagInfo"]["TotalMagOnAtom"]]
            else:
                check = 0
        else:
            FinalMag = np.zeros(len(Elements), dtype=np.float64).tolist()
        # self.Force = []
        # self.RelaxTraj = []

        # 读取结构弛豫信息
        for i in range(1,data["Structures"]["FinalStep"][0]+1):
            key = "Step-" + str(i)
            # 首先读取受力信息
            force = [x.astype(np.float64) for x in data["Structures"][key]["Force"]]
            force = np.array(force)
            force = force.reshape(-1, 3)
            self.Force.append(force)
            # 读取原子磁矩信息 mag
            if "Mag" in data["Structures"][key]:
                mag = [x.astype(np.float64) for x in data["Structures"][key]["Mag"]]
                if len(mag) == len(Elements):
                    pass
                else:
                    raise ValueError("The spin is not collinear!")
                mag = np.array(mag)
            else:
                mag = np.zeros(len(Elements), dtype=int).tolist()
            # 读取晶胞信息 lattice_relax
            lattice_relax = [x.astype(np.float64) for x in data["Structures"][key]["Lattice"]]
            lattice_relax = np.array(lattice_relax)
            lattice_relax = lattice_relax.reshape(-1, 3)
            # 读取原子位置信息 position_relax
            position_relax = [x.astype(np.float64) for x in data["Structures"][key]["Position"]]
            position_relax = np.array(position_relax)
            position_relax = position_relax.reshape(-1, 3)
            # 虽然它写的是笛卡尔坐标系，但是实际上是分数坐标系
            position_relax = np.dot(position_relax, lattice_relax)
            # 构造当前步骤的原子类 atoms_relax
            atoms_relax = Atoms(cell = lattice_relax, pbc = True)
            for j in range(len(Elements)):
                atom = Atom(symbol = Elements[j], 
                            position = position_relax[j], 
                            mass = molar_mass(Elements[j]),
                            magmom = mag[j], tag = Fix[j])
                atoms_relax.append(atom)
            self.RelaxTraj.append(atoms_relax) # 将当前步骤的原子类加入弛豫轨迹

        # 最终结构信息
        self.atoms_final = self.RelaxTraj[-1]
        # 校验磁矩信息
        if check == 1:
            if sum(self.atoms_final.get_initial_magnetic_moments() - FinalMag) != 0:
                raise ValueError("The magnetic moments are not consistent!(code bug, fix it!)")
        if "Electron" in data:
            self.electron_number = data["Electron"][0].astype(np.float64)
            # 读取能量信息
            self.energy = {
                "EFermi": data["Energy"]["EFermi"][0].astype(np.float64),
                "TotalEnergy": data["Energy"]["TotalEnergy"][0].astype(np.float64),
                "TotalEnergy0": data["Energy"]["TotalEnergy0"][0].astype(np.float64),
                "VdwCorrection": data["VdwCorrection"][0].astype(np.float64),
                "ElectronNumber": self.electron_number
            }
        else:
            self.electron_number = "Unknown"
            self.energy = None

    def write_as_file(self, path=os.getcwd(), filename='structure_new.as', save_kind = 'final'):
        '''- write the as file
        - path: the path to save the file
        - file: the name of the file
        - save_kind: the kind of structure to save, 'initial', 'final' or the index of the structure in the RelaxTraj
        '''
        latt_form = ' %19.6f'
        cform = ' %19.16f'
        element_form = ' %2s'
        if save_kind == 'initial':
            self.atoms = self.atoms_initial
        elif save_kind == 'final':
            self.atoms = self.atoms_final
        elif isinstance(save_kind, int):
            self.atoms = self.RelaxTraj[save_kind]
        with open(path + '\\' + filename, 'w') as f:
            f.write("Total number of atoms" + "\n")
            f.write(str(len(self.atoms)) + "\n")
            f.write("Lattice" + "\n")
            for i in range(3):
                f.write(latt_form % self.atoms.cell[i][0] + " " + latt_form % self.atoms.cell[i][1] + " " + latt_form % self.atoms.cell[i][2] + "\n")
            f.write("Cartesian Mag")
            if self.atoms[0].tag != -1:
                f.write(" Fix_x Fix_y Fix_z" + "\n")
            else:
                f.write("\n")
            for i in range(len(self.atoms)):
                f.write(element_form % self.atoms[i].symbol + " ")
                f.write(cform % self.atoms[i].position[0] + " " + cform % self.atoms[i].position[1] + " " + cform % self.atoms[i].position[2] + " ")
                f.write(cform % self.atoms[i].magmom)
                if self.atoms[0].tag != -1:
                    f.write(" " + tran(self.atoms[i].tag) + "\n")
                else:
                    f.write("\n")

class pchargeLoad:
    """def __init__(self, filename = "pcharge.h5"):
    filename;
    atoms;
    chgInfo
    bandIndex
    grid
    """
    
    filename = None
    atoms = None
    chgInfo = {}
    bandIndex = {}
    grid = None
    # position = []
    # elements = []
    def __init__(self, filename = "pcharge.h5"):
        self.filename = filename
        self.atoms = self.__structureLoad()
        self.__pchargeLoad()
    def __structureLoad(self):
        data = h5py.File(self.filename, 'r')
        CoordinateType = [x.decode('UTF-8') for x in data["AtomInfo"]["CoordinateType"]]
        CoordinateType = "".join(CoordinateType)
        Elements = [x.decode('UTF-8') for x in data["AtomInfo"]["Elements"]]
        temp1 = [] # 将Elements中的符号转换为正常的格式
        temp2 = []
        for letter in Elements:
            if letter != ";":
                temp2.append(letter)
            elif letter == ";":
                temp1.append("".join(temp2))
                temp2 = []
        temp1.append("".join(temp2))
        Elements = temp1
        # self.elements = Elements
        # 读取结构信息
        if "Fix" in data["AtomInfo"]:
            Fix = np.array([x.astype(np.float64) for x in data["AtomInfo"]["Fix"]])
            Fix = Fix.reshape(-1, 3)
            Fix = Fix[:, 0]*100 + Fix[:, 1] * 10 + Fix[:, 2]
        else:
            Fix = np.zeros(len(data["AtomInfo"]["Position"]), dtype=int) - 1 # 如果文件中没有Fix信息，则所有固定值设为-1
        if "Mag" in data["AtomInfo"]:
            Mag = [x.astype(np.float64) for x in data["AtomInfo"]["Mag"]]
        else:
            Mag = np.zeros(len(Elements), dtype=np.float64).tolist()
        Position = [x.astype(np.float64) for x in data["AtomInfo"]["Position"]]
        Position = np.array(Position)
        Position = Position.reshape(-1, 3)
        # self.position = Position
        Lattice = [x.astype(np.float64) for x in data["AtomInfo"]["Lattice"]]
        Lattice = np.array(Lattice)
        Lattice = Lattice.reshape(-1, 3)
        # 构造原子类
        atoms = Atoms(cell = Lattice, pbc = True)
        if CoordinateType == "Cartesian":
            pass
        else:
            Position = Position.dot(Lattice)
        for i in range(len(Elements)):
            atom = Atom(symbol = Elements[i], 
                        position = Position[i], 
                        mass = molar_mass(Elements[i]),
                        magmom = Mag[i], tag = Fix[i])
            atoms.append(atom)
        return atoms
    def __pchargeLoad(self):
        data = h5py.File(self.filename, 'r')
        self.grid = data['AtomInfo']['Grid'][:]
        for i in data["Pcharge"]:
            try:
                float(i)
                self.chgInfo[i] = data["Pcharge"][i]['TotalCharge'][:]
                self.bandIndex[i] = data["Pcharge"][i]['BandIndex'][:]
            except:
                pass

class scfLoad:
    """def __init__(self, filename = "pcharge.h5"):
    filename;
    atoms;
    chgInfo
    bandIndex
    grid
    """
    filename = None
    atoms = None
    bandEnergies_Sorted_spin1 = {}
    bandEnergies_Sorted_spin2 = {}
    bandsIndex_Sorted_spin1 = {}
    bandsIndex_Sorted_spin2 = {}
    bandOccupation_Sorted_spin1 = {}
    bandOccupation_Sorted_spin2 = {}
    grid = None
    kpointsNumber = None
    bandsNumber = None
    # position = []
    # elements = []
    def __init__(self, filename = "scf.h5"):
        self.filename = filename
        self.atoms = self.__structureLoad()
        self.__bandLoad()
    def __structureLoad(self):
        data = h5py.File(self.filename, 'r')
        CoordinateType = [x.decode('UTF-8') for x in data["AtomInfo"]["CoordinateType"]]
        CoordinateType = "".join(CoordinateType)
        Elements = [x.decode('UTF-8') for x in data["AtomInfo"]["Elements"]]
        temp1 = [] # 将Elements中的符号转换为正常的格式
        temp2 = []
        for letter in Elements:
            if letter != ";":
                temp2.append(letter)
            elif letter == ";":
                temp1.append("".join(temp2))
                temp2 = []
        temp1.append("".join(temp2))
        Elements = temp1
        # self.elements = Elements
        # 读取结构信息
        if "Fix" in data["AtomInfo"]:
            Fix = np.array([x.astype(np.float64) for x in data["AtomInfo"]["Fix"]])
            Fix = Fix.reshape(-1, 3)
            Fix = Fix[:, 0]*100 + Fix[:, 1] * 10 + Fix[:, 2]
        else:
            Fix = np.zeros(len(data["AtomInfo"]["Position"]), dtype=int) - 1 # 如果文件中没有Fix信息，则所有固定值设为-1
        if "Mag" in data["AtomInfo"]:
            Mag = [x.astype(np.float64) for x in data["AtomInfo"]["Mag"]]
        else:
            Mag = np.zeros(len(Elements), dtype=np.float64).tolist()
        Position = [x.astype(np.float64) for x in data["AtomInfo"]["Position"]]
        Position = np.array(Position)
        Position = Position.reshape(-1, 3)
        # self.position = Position
        Lattice = [x.astype(np.float64) for x in data["AtomInfo"]["Lattice"]]
        Lattice = np.array(Lattice)
        Lattice = Lattice.reshape(-1, 3)
        # 构造原子类
        atoms = Atoms(cell = Lattice, pbc = True)
        if CoordinateType == "Cartesian":
            pass
        else:
            Position = Position.dot(Lattice)
        for i in range(len(Elements)):
            atom = Atom(symbol = Elements[i], 
                        position = Position[i], 
                        mass = molar_mass(Elements[i]),
                        magmom = Mag[i], tag = Fix[i])
            atoms.append(atom)
        #/AtomInfo/Grid
        self.grid = [int(x) for x in data["AtomInfo"]["Grid"]]
        return atoms
    def __bandLoad(self):
        # 首先读取spin1的数据
        data = h5py.File(self.filename, 'r')
        temp = np.array(data["Eigenvalue"]["Spin1"]["BandEnergies"])
        self.bandsNumber = temp.shape[0]
        self.kpointsNumber = temp.shape[1]
        bandEnergies_spin1 = {}
        for i in range(self.kpointsNumber):
            bandEnergies_spin1[i+1] = temp[:, i]
            self.bandsIndex_Sorted_spin1[i+1] = np.argsort(bandEnergies_spin1[i+1])+1
            self.bandEnergies_Sorted_spin1[i+1] = bandEnergies_spin1[i+1][self.bandsIndex_Sorted_spin1[i+1]-1]
            bandsOccupation_spin1 = np.array(data["Eigenvalue"]["Spin1"]["Occupation"])
            self.bandOccupation_Sorted_spin1[i+1] = bandsOccupation_spin1[:, i][self.bandsIndex_Sorted_spin1[i+1]-1]
        # 读取spin2的数据
        bandEnergies_spin2 = {}
        if "Spin2" in data["Eigenvalue"]:
            temp = np.array(data["Eigenvalue"]["Spin2"]["BandEnergies"])
            for i in range(self.kpointsNumber):
                bandEnergies_spin2[i+1] = temp[:, i]
                self.bandsIndex_Sorted_spin2[i+1] = np.argsort(bandEnergies_spin2[i+1])+1
                self.bandEnergies_Sorted_spin2[i+1] = bandEnergies_spin2[i+1][self.bandsIndex_Sorted_spin2[i+1]-1]
                bandsOccupation_spin2 = np.array(data["Eigenvalue"]["Spin2"]["Occupation"])
                self.bandOccupation_Sorted_spin2[i+1] = bandsOccupation_spin2[:, i][self.bandsIndex_Sorted_spin2[i+1]-1]
        else:
            self.bandEnergies_spin2 = None
            self.bandsIndex_Sorted_spin2 = None
            self.bandOccupation_Sorted_spin2 = None
    def write_bands_occupation(self):
        latt_form = ' %19.6f'
        with open("bands_occupation_spin1.txt", "w") as f:
            f.write("kpointsNumber = %d\n" % self.kpointsNumber)
            for i in range(self.kpointsNumber):
                f.write(" %19s" % "Index/KPoints{}".format(str(i+1)))
                f.write(" %19s" % "Occupation")
                f.write(" %19s" % "Energy")
                f.write("%5s" % "|")
            f.write("\n")
            for i in range(self.bandsNumber):
                for j in range(self.kpointsNumber):
                    f.write(' %19s' % self.bandsIndex_Sorted_spin1[j+1][i])
                    f.write(latt_form % self.bandOccupation_Sorted_spin1[j+1][i])
                    f.write(latt_form % self.bandEnergies_Sorted_spin1[j+1][i])
                    f.write("%5s" % "|")
                f.write("\n")
        if self.bandEnergies_Sorted_spin2 != None:
            with open("bands_occupation_spin2.txt", "w") as f:
                f.write("kpointsNumber = %d\n" % self.kpointsNumber)
                for i in range(self.kpointsNumber):
                    f.write(" %19s" % "Index/KPoints{}".format(str(i+1)))
                    f.write(" %19s" % "Occupation")
                    f.write(" %19s" % "Energy")
                    f.write("%5s" % "|")
                f.write("\n")
                for i in range(self.bandsNumber):
                    for j in range(self.kpointsNumber):
                        f.write(' %19s' % self.bandsIndex_Sorted_spin2[j+1][i])
                        f.write(latt_form % self.bandOccupation_Sorted_spin2[j+1][i])
                        f.write(latt_form % self.bandEnergies_Sorted_spin2[j+1][i])
                        f.write("%5s" % "|")
                    f.write("\n")

