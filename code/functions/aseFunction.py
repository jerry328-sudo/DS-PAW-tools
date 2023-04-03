# %%
from ase import Atoms, Atom
from ase.visualize import view
import re

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