import sys
sys.path.append("f:\\code\\DS-PAW\\DS-PAW-tools\\code")
import functions.function_normal as fn
import numpy as np

cell = fn.Read_As_File()

print('Tne number of atoms is: ', str(cell.n_atoms))
print('The atom type is: ', str(cell.atom_type))

atom_type = input('Please input the atom type you want to replace: ')
atom_type = atom_type.capitalize()
while atom_type not in cell.atom_type.keys():
    atom_type = input('The atom type is not in the as file, pleasse input again: ')
    atom_type = atom_type.capitalize()

atom_doped = input('Please input the atom type you want to replace with: ')
atom_doped = atom_doped.capitalize()

# 当报错时循环执行
indicate = 1
while indicate == 1:
    try:
        n_doped = int(input('Input random doping quantity (enter "-1" to enter proportion): '))
        indicate = 0
    except:
        print('Please input an integer!')


if n_doped == -1:
    n_doped = float(input('''Input random doping proportion (0-1): 
eg. the atom you want to replace is Cu, and you want to replace 10% of Cu with Ag, then you should input 0.1: '''))
    n_doped = cell.atom_type[atom_type]*n_doped
    if n_doped % 1 != 0:
        print('The number of atoms to be doped is: ', str(n_doped))
        n_doped = input(
            'The proportion is not an integer, please re-enter (the number): ')

n_doped = int(n_doped)
atom_num = np.random.choice(cell.atom_type[atom_type], n_doped, replace=False)
ia = 0
for i in range(cell.target+1, len(cell.I)):
    if cell.Ic[i][0] == atom_type:
        if ia in atom_num:
            cell.Ic[i][0] = atom_doped
        ia += 1

cell.write_as_file()
