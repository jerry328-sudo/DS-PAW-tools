# %%
import numpy as np
import scipy.constants as sc
from ase.thermochemistry import IdealGasThermo
from ase.visualize import view
import ase.units as units
from pymatgen.core.structure import Molecule
from pymatgen.symmetry.analyzer import PointGroupAnalyzer
import os
import sys
current_path = os.path.dirname(os.path.abspath(__file__))
parent_path = os.path.dirname(current_path)
sys.path.append(parent_path)
import functions.aseFunction as af



# 首先是读取frequency.json文件，构建热力学对象
class thermochemistry():
    atoms = None
    symmetry_ckeck = None
    symmetrynumber = None
    symmetry = None
    gibbs_energy = None
    entropy = None
    enthalpy = None
    internal_energy = None
    zpe = 0
    freq = None
    freq_list = None
    geometry = None
    def __init__(self, file_path='frequency.h5', geometry='nonlinear'):
        self.geometry = geometry
        self.__read_data(file_path)

    def __read_data(self, file_path):
        # 读取原子信息
        self.atoms = af.freLoad(file_path)
        # 检查对称性
        mol = Molecule(species=self.atoms.atoms.get_chemical_symbols(),
                    coords=self.atoms.atoms.get_positions())
        self.symmetry_ckeck = PointGroupAnalyzer(mol)
        self.symmetrynumber = self.symmetry_ckeck.get_rotational_symmetry_number()
        #len(self.symmetry_ckeck.get_symmetry_operations())是错误的，应该是旋转对称数
        self.symmetry = self.symmetry_ckeck.get_pointgroup()

        # 读取频率信息
        self.freq = self.atoms.Frequency
        self.freq_list = self.atoms.Frequency
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
                        atoms=self.atoms.atoms,
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

print('-->> (01) Reading frequency.h5 File...', end=' ')

if not os.path.exists('frequency.h5'):
    print('Failed!')
    print('the file "frequency.h5" does not exist!')
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



