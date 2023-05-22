#!/usr/bin/env python
import os
import sys
current_path = os.path.dirname(os.path.abspath(__file__))
parent_path = os.path.dirname(current_path)
sys.path.append(parent_path)
import functions.aseFunction as af
import numpy as np

# default accuracy = 0.04
# default accuracy factor = 1
accuracy = input('Please enter the accuracy of kpoints (default: 0.04): ')
factor = input('Please enter the accuracy factor of kpoints (default: 1): ')
if accuracy == '':
    accuracy = 0.04
if factor == '':
    factor = 1
accuracy = float(accuracy) / float(factor)

structure = af.read_as()
# 获取晶胞的晶格常数
lattice = structure.cell.cellpar()
# 获取倒易晶格的晶格常数
reciprocal_lattice = structure.cell.reciprocal().cellpar()
# 将倒易晶格常数除以0.04后四舍五入,如果为0则取1
kmesh = np.round(reciprocal_lattice[:3] / accuracy)
kmesh[kmesh == 0] = 1
kmesh = kmesh.astype(int)
print(kmesh)