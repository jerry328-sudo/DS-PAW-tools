import os
import sys
sys.path.append("f:\\code\\DS-PAW\\DS-PAW-tools\\code")
import functions.function_normal as fn

methode = input('Select a fixed method (by element (0) or by position (1)): ')
while methode not in ['0', '1']:
    methode = input('Select a fixed method (by element (0) or by position (1), -1 is exit): ')
    if methode == '-1':
        exit()
methode = int(methode)
if methode == 0:
    elements = []
    element = '0'
    while True:
        element = input('Please select the element of the fixed atom (enter -1 to end): ')
        element = element.split()
        for lenele in range(len(element)):
            elements.append(element[lenele].capitalize())
        if element == ['-1']:
            break
elif methode == 1:
    direction = input('Please select the direction of the fixed atom (x, y, z): ')
    direction = direction.lower()

    while direction not in ['x', 'y', 'z']:
        print('Please enter x, y, or z!')
        direction = input('Please select the direction of the fixed atom (x, y, z): ')
        direction = direction.lower()
    direct = direction
    print(' ')

    if direction == 'x':
        direction = 0
    elif direction == 'y':
        direction = 1
    elif direction == 'z':
        direction = 2
    thick = input('Enter the atom below which coordinate to fix: ')
    thick = float(thick)
filename = input('Enter the name of the output file (Enter 0 to use the original name): ')
filename_check = filename


directory = os.getcwd()
# Iterate over all the files in the directory
for filename in os.listdir(directory):
    # check whether the file is a directory
    if not os.path.isdir(directory + '\\' + filename):
        continue
    # check whether filename can transform to int
    try:
        int(filename)
    except ValueError:
        continue
    path = directory + '\\' + filename
    os.chdir(path)
    structure = fn.Read_As_File(path = path, file = 'structure{}.as'.format(filename))
    if filename_check == '0':
        filename = 'structure{}.as'.format(filename)
    if methode == 0:
        structure.set_fix_atom(methode=methode+1, fixed_element=elements, path=path, filename=filename)
    elif methode == 1:
        structure.set_fix_atom(methode=methode+1, direction=direction, thick=thick, path=path, filename=filename)