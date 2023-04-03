# %%
import numpy as np
import scipy.constants as sc
from ase.thermochemistry import IdealGasThermo
from ase import Atoms, Atom
from ase.visualize import view
import ase.units as units
from pymatgen.core.structure import Molecule
from pymatgen.symmetry.analyzer import PointGroupAnalyzer
import json
import os

# %%
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


# 首先是读取frequency.json文件，构建热力学对象
class thermochemistry():
    def __init__(self, file_path='frequency.json', geometry='nonlinear'):
        self.atoms = None
        self.symmetry_ckeck = None
        self.symmetrynumber = None
        self.symmetry = None
        self.gibbs_energy = None
        self.entropy = None
        self.enthalpy = None
        self.internal_energy = None
        self.zpe = 0
        self.freq = None
        self.freq_list = None
        self.geometry = geometry
        self.__read_data(file_path)

    def __read_data(self, file_path):
        # 读取原子信息
        with open(file_path, 'r') as f:
            data = json.load(f)
        self.atoms = Atoms(
            cell=np.array(data["AtomInfo"]["Lattice"]).reshape(3, 3).tolist(),
            pbc=True)
        for a in data["AtomInfo"]["Atoms"]:
            temp = Atom(symbol=a['Element'], position=a['Position'], mass = molar_mass(a['Element']))
            self.atoms.append(temp)
        # 检查对称性
        mol = Molecule(species=self.atoms.get_chemical_symbols(), 
                    coords=self.atoms.get_positions())
        self.symmetry_ckeck = PointGroupAnalyzer(mol)
        self.symmetrynumber = self.symmetry_ckeck.get_rotational_symmetry_number()
        #len(self.symmetry_ckeck.get_symmetry_operations())是错误的，应该是旋转对称数
        self.symmetry = self.symmetry_ckeck.get_pointgroup()

        # 读取频率信息
        self.freq = np.array(data['FrequencyInfo'])
        self.freq_list = []
        for i in self.freq:
            self.freq_list.append(i['eigenvalues']['frequency'])
        # 将单位从THz转为eV
        self.freq_list = np.array(self.freq_list) * sc.tera * sc.h * 1 / sc.e
        # For Gas, the last six are translation and rotation
        '''因为线性分子的振动自由度为(3N-5)。发现有三个虚频和两个非常小的频率，
        这是平动和转动在振动自由度上的投影，直接忽略即可(VASPKIT自动判断并忽略)。
        忽略最小的5个(线型分子)或6个(非线型分子)振动频率，并不是直接忽略了平动和转动的
        贡献。而是通过平动和转动的配分函数另外计算其对热力学量的贡献。其中平动熵是气体
        分子熵的主要贡献。'''
        # sort the frequency and return the index
        index = np.argsort(self.freq_list)

        if self.geometry == 'linear':
            self.freq_list = np.delete(self.freq_list, index[:5])
            # self.freq_list = self.freq_list[:-5]
        elif self.geometry == 'nonlinear':
            self.freq_list = np.delete(self.freq_list, index[:6])
            # self.freq_list = self.freq_list[:-6]
        else:
            raise ValueError('geometry must be linear or nonlinear')

    def thermochemistry_caculate(self, T, spin=0, pressure = 0):
        '''Calculate the thermochemistry of adsorption molecule, the nuit is eV'''
        # 计算分子的热力学参数
        for energy in self.freq_list:
            self.zpe += 0.5 * energy
        self.thermo = IdealGasThermo(vib_energies=self.freq_list,
                        atoms=self.atoms,
                        symmetrynumber=self.symmetrynumber,
                        geometry=self.geometry,
                        spin=spin)
        pressure = pressure * sc.bar
        self.gibbs_energy = self.thermo.get_gibbs_energy(T, pressure=pressure, verbose=False)
        self.entropy = self.thermo.get_entropy(T, pressure=pressure, verbose=False)
        self.enthalpy = self.thermo.get_enthalpy(T, verbose=False)
        # U = H - pV = H - nRT
        self.internal_energy = self.enthalpy - sc.R * T / sc.N_A * units.J

# %%
print(""" +-------------------------- Warm Tips --------------------------+
   Included Vibrations, Translation, Rotation & Electron contributions.       
 GAS molecules should not be with any fix.       """)

geometry = int(input(''' Please input the geometry of the molecule: 
 (0: 'linear', 1: 'nonlinear')
 WARNING! The 'monatomic' is not available!\n'''))
if geometry == 0:
    geometry = 'linear'
elif geometry == 1:
    geometry = 'nonlinear'
else:
    print('Wrong input!')
    exit()

print('-->> (01) Reading frequency.json File...', end=' ')

if not os.path.exists('frequency.json'):
    print('Failed!')
    print('the file "frequency.json" does not exist!')
    exit()
else:
    thermo = thermochemistry(geometry=geometry)
    print('Done!')

print("""  -->> (02) Analyzing Molecular Symmetry Information...Done!
 Molecular Symmetry is: {}
 rotational symmetry number is: {}
 +---------------------------------------------------------------+""".format(thermo.symmetry, thermo.symmetrynumber))

tempuerature = float(input(' Please input the temperature (K): \n'))
pressure = float(input(' Please input the pressure (Atm): \n'))
spin = float(input(''' Please input the total electronic spin: 
 (the total electronic spin. (0 for molecules in which all electrons are paired, 
 0.5 for a free radical with a single unpaired electron, 1.0 for a triplet with 
 two unpaired electrons, such as O_2.))\n'''))

print('-->> (03) Calculating Thermochemistry...', end=' ')
thermo.thermochemistry_caculate(T=tempuerature, 
                                spin=spin, 
                                pressure=pressure)
print('Done!\n')
print(""" +---------------------------------------------------------------+""")
print(""" U(T) = ZPE + Delta_U(0->T)
 H(T) = U(T) + PV = ZPE + Delta_U(0->T) + PV
 G(T) = H(T) - TS = ZPE + Delta_U(0->T) + PV - TS
""")
string = """ Temperature (T)           :      {} K
 Pressure (P)              :      {} Pa
 Zero-point energy E_ZPE   :      {} eV
 Thermal correction to U(T):      {} eV
 Thermal correction to H(T):      {} eV
 Thermal correction to G(T):      {} eV
 Entropy S                 :      {} eV/K
 Entropy contribution T*S  :      {} eV
""".format(tempuerature, pressure * sc.atm, thermo.zpe, thermo.internal_energy , thermo.enthalpy, 
        thermo.gibbs_energy, thermo.entropy, thermo.entropy * tempuerature)
print(string)

with open('thermoCorrection.txt', 'w') as f:
    f.write(""" Molecular Symmetry: {}
 rotational symmetry number: {}
 +---------------------------------------------------------------+\n""".format(thermo.symmetry, thermo.symmetrynumber))
    f.write(string)

# %%



