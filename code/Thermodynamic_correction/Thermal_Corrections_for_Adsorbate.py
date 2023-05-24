# %%
# 吸附分子配分函数计算
import numpy as np
import scipy.constants as sc
from ase.thermochemistry import HarmonicThermo
import os
import sys
current_path = os.path.dirname(os.path.abspath(__file__))
parent_path = os.path.dirname(current_path)
sys.path.append(parent_path)
import functions.aseFunction as af


class thermochemistry():
    def __init__(self, filename='frequency.h5'):
        '''- Read the frequency file frequency.json
        if the user assumes the pV term (in H = U + pV) is zero 
        this can also be interpreted as the Gibbs free energy
        
        - Neglect PV contribution to translation for adsorbed molecules.
        To avoid abnormal entropy contribution,
        frequencies less than 50 cm-1 are set to 50 cm-1.'''
        freInfo = af.freLoad(filename)
        freInfo.write_jmol() # 生成jmol文件
        self.freq_list = []
        for i in freInfo.Frequency:
        # 将小于50cm-1的频率设置为50cm-1
            if i < 1.5:
                self.freq_list.append(1.5)
            else:
                self.freq_list.append(i)
        # self.freq_list = freInfo.Frequency
        # 将单位从THz转为eV
        self.freq_list = np.array(self.freq_list) * sc.tera * sc.h * 1 / sc.e
        self.gibbs_energy = None
        self.internal_energy = None
        self.entropy = None
        self.zpe = 0

    def thermochemistry_caculate(self, T):
        '''Calculate the thermochemistry of adsorption molecule, the nuit is eV'''
        # 计算分子的热力学参数
        self.thermo = HarmonicThermo(self.freq_list)
        self.gibbs_energy = self.thermo.get_helmholtz_energy(
            temperature=T, verbose=False)
        self.internal_energy = self.thermo.get_internal_energy(
            temperature=T, verbose=False)
        self.entropy = self.thermo.get_entropy(temperature=T, verbose=False)
        for energy in self.freq_list:
            self.zpe += 0.5 * energy


# %%
if not os.path.exists('frequency.json'):
    print('the file frequency.json does not exist!')
    exit()
else:
    print('-->> (01) Reading frequency.json File...')
    thermo = thermochemistry()

tempurature = input('''Please Enter The Temperature (K): 
 
 ------------>>
''')

thermo.thermochemistry_caculate(T = float(tempurature))

string = """
 U(T) = H(T) = ZPE + Delat_U(0->T)
 G(T) = H(T) - TS = ZPE + Delat_U(0->T) - TS
 
 Temperature (T): {} K
 Zero-point energy E_ZPE   :       {} eV
 Thermal correction to U(T):       {} eV
 Thermal correction to H(T):       {} eV
 Thermal correction to G(T):       {} eV
 Entropy S                 :       {} eV/K
 Entropy contribution T*S  :       {} eV
""".format(tempurature, thermo.zpe, thermo.internal_energy, thermo.internal_energy, 
           thermo.gibbs_energy, thermo.entropy, thermo.entropy * float(tempurature))

print(string)

with open('Thermal_Corrections_for_Adsorbate.txt', 'w') as f:
    f.write(string)
