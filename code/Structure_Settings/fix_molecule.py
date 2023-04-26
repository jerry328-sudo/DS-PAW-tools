import os
import sys
current_path = os.path.dirname(os.path.abspath(__file__))
parent_path = os.path.dirname(current_path)
sys.path.append(parent_path)
import functions.function_normal as fn

dir = os.getcwd()
# determine whether the file 'structure.as' exists
print('reading structure.as ...', end=' ')
if not os.path.exists('structure.as'):
    print('structure.as does not exist!')
    exit()
else:
    structure = fn.Read_As_File(path=dir)
    print('done')

# 将固定坐标转换为数字


def direction_to_number(direction):
    if direction == 'x':
        direction = 0
    elif direction == 'y':
        direction = 1
    elif direction == 'z':
        direction = 2
    return direction


# %% 首先确定固定分子的方式
# Select a fixed method (by element or by position)
type = input('Select a fixed method (by element (0) or by position (1)): ')
while type not in ['0', '1']:
    type = input(
        'Select a fixed method (by element (0) or by position (1), -1 is exit): ')
    if type == '-1':
        exit()

count = 0

# %% 如果固定分子的方式为按元素
if type == '0':
    # Select the element of the fixed atom
    elements = []
    element = input('Please select the element of the fixed atom:(eg. ni co) ')
    element = element.split()
    for i in element:
        elements.append(i.capitalize())
    print(' ')
    print('The fixed atoms are: ', elements)
    # 选择固定的坐标轴
    directions = []
    direction = input(
        'Please select the direction of the fixed atom (x, y, z): ')
    direction = direction.split()
    for i in direction:
        directions.append(i.lower())
        # directions.append(direction_to_number(i.lower()))
    print(' ')
    print('The fixed directions are: ', directions)

    # 开始固定原子坐标
    judge = 0  # 判断是否有固定的原子, 0表示没有, 1表示有
    for i in range(len(structure.Ic)):
        if i == structure.target:
            if 'Fix' in structure.I[i]:  # 如果已经有固定的原子, 则不再固定
                judge = 1
            else:
                judge = 0
                structure.Ic.append('Fix_x Fix_y Fix_z')
        if i > structure.target:
            if judge == 0:  # 首先是先前没有固定过原子的情况
                if structure.Ic[i][0] in elements:
                    count += 1
                    if 'x' in directions:
                        structure.Ic[i].append('T')
                    else:
                        structure.Ic[i].append('F')
                    if 'y' in directions:
                        structure.Ic[i].append('T')
                    else:
                        structure.Ic[i].append('F')
                    if 'z' in directions:
                        structure.Ic[i].append('T')
                    else:
                        structure.Ic[i].append('F')
            if judge == 1:  # 其次是先前已经固定过原子的情况
                if structure.Ic[i][0] in elements:
                    count += 1
                    if 'x' in directions:
                        structure.Ic[i][-3] = 'T'
                    if 'y' in directions:
                        structure.Ic[i][-2] = 'T'
                    if 'z' in directions:
                        structure.Ic[i][-1] = 'T'
# %% 如果固定分子的方式为按位置
elif type == '1':
    # 选择判断用的坐标轴
    coordinate = input(
        'Please select the coordinate of the fixed atom (x, y, z): ')
    while coordinate not in ['x', 'y', 'z']:
        coordinate = input(
            'Please select the coordinate of the fixed atom (x, y, z): ')
    coordinate = coordinate.lower()
    if coordinate == 'x':
        coordinate = 1
    elif coordinate == 'y':
        coordinate = 2
    elif coordinate == 'z':
        coordinate = 3
    # 设置固定原子的范围
    range_a = input('Please select the range of the fixed atom (eg. 0 1): ')
    range_a = range_a.split()
    range_a = [float(i) for i in range_a]
    # 选择固定的坐标轴
    directions = []
    direction = input(
        'Please select the direction of the fixed atom (x, y, z): ')
    direction = direction.split()
    for i in direction:
        directions.append(i.lower())

    # 开始固定原子坐标
    judge = 0  # 判断是否有固定的原子, 0表示没有, 1表示有
    for i in range(len(structure.Ic)):
        if i == structure.target:
            if 'Fix' in structure.I[i]:  # 如果已经有固定的原子, 则不再固定
                judge = 1
            else:
                judge = 0
                structure.Ic.append('Fix_x Fix_y Fix_z')
        if i > structure.target:
            coor = float(structure.Ic[i][coordinate])
            if judge == 0:
                if range_a[0] <= coor and coor <= range_a[1]:
                    count += 1
                    if 'x' in directions:
                        structure.Ic[i].append('T')
                    else:
                        structure.Ic[i].append('F')
                    if 'y' in directions:
                        structure.Ic[i].append('T')
                    else:
                        structure.Ic[i].append('F')
                    if 'z' in directions:
                        structure.Ic[i].append('T')
                    else:
                        structure.Ic[i].append('F')
            if judge == 1:
                if range_a[0] <= coor and coor <= range_a[1]:
                    count += 1
                    if 'x' in directions:
                        structure.Ic[i][-3] = 'T'
                    if 'y' in directions:
                        structure.Ic[i][-2] = 'T'
                    if 'z' in directions:
                        structure.Ic[i][-1] = 'T'

# %% 保存固定原子的坐标
filename = input('Enter the name of the output file (Enter 0 to use the original name): ')
if filename == '0':
    filename = 'structure'
structure.write_as_file(path = dir, file = '{}.as'.format(filename))

with open('fix_molecule.txt', 'w') as f:
    f.write('The number of fixed atoms is: {}\r'.format(count))