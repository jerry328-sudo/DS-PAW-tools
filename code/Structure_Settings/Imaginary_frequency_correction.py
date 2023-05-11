# -*- coding: utf-8 -*-
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

data = af.freLoad()
l_position = data.imFreLocate  #虚频振动方向部分在OUTCAR中的起始行数
correction_factor = 0.1
inputs = input("Please input the correction factor: (default: 0.1)")
if inputs != '':
    correction_factor = float(inputs)
else:
    pass
print("The correction factor is: ", str(correction_factor))

if data.Unit[0] == 'THz':
    pass
else:
    raise ValueError('Unit of frequency is not THz')

# 思路就是把最大虚频对应的特征向量乘以校正因子，然后直接加到原子坐标上。
max_imfre_index = l_position
for i in range(len(data.Frequency)):
    if data.FrequencyType[i] == "f/i":
        if data.Frequency[i] > data.Frequency[max_imfre_index]:
            max_imfre_index = i
        else:
            pass
    else:
        pass
data.atoms.positions += data.EigenVectors[max_imfre_index] * correction_factor

data.write_as_file(file="structure_corr.as")