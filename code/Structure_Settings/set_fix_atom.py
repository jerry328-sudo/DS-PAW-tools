import os
import sys
sys.path.append("f:\\code\\DS-PAW\\DS-PAW-tools\\code")
import functions.function_normal as fn

# get the directory where the py file is located
# dir = os.path.dirname(os.path.realpath(__file__))
# os.chdir(dir)
dir = os.getcwd()

# determine whether the file 'structure.as' exists
print('reading structure.as ...', end=' ')
if not os.path.exists('structure.as'):
    print('structure.as does not exist!')
    exit()
else:
    structure = fn.Read_As_File(path = dir)
    print('done')


# Select a fixed method (by element or by position)
type = input('Select a fixed method (by element (0) or by position (1)): ')
while type not in ['0', '1']:
    type = input('Select a fixed method (by element (0) or by position (1), -1 is exit): ')
    if type == '-1':
        exit()

if type == '1':
    # Select the direction of the fixed atom
    direction = input('Please select the direction of the fixed atom (x, y, z): ')
    direction = direction.lower()

    while direction not in ['x', 'y', 'z']:
        print('Please enter x, y, or z!')
        direction = input('Please select the direction of the fixed atom (x, y, z): ')
        direction = direction.lower()
    direct = direction
    print(' ')

    if direction == 'x':
        direction = 1
    elif direction == 'y':
        direction = 2
    elif direction == 'z':
        direction = 3

    dcount = 0
    thick = []
    target = 0

    for ic in structure.Ic:
        if 'Cartesian' in ic:
            target = 1
        if target == 1:
            if 'Cartesian' not in ic:
                temp = round(float(ic[direction]), 3)
                if temp not in thick:
                    thick.append(temp)
    thick = sorted(thick, reverse=True)
    if len(thick) == 1:
        print('There is only one layer of atoms!')
        exit()
    elif len(thick) > 50:
        temp = input('There are more than 50 layers of atoms, do you want to continue? (y/n): ')
        temp = temp.lower()
        if temp == 'n':
            exit()
        elif temp == 'y':
            pass
        else:
            if dcount == 0:
                print('Please select the correct direction!')
                dcount += 1
            else:
                pass
            pass

    print('{} coordinate of the atoms:'.format(direct)+'\n')

    for thicks in thick:
        print(thicks)

    print('  ')

    dire = input('Enter the atom below which coordinate to fix: ')
    dire = float(dire)

    target = 0
    temp = ['F', 'F', 'F']
    tempt =['T', 'T', 'T']
    count_free_atom = 0
    count_fix_atom = 0
    count_atom = 0


    for i in range(len(structure.Ic)):
        if 'Cartesian' in structure.Ic[i]:
            if 'Fix' not in structure.I[i]:
                structure.Ic[i].append('Fix_x Fix_y Fix_z')
                target = 1
            else:
                target = 2
        if target == 1:
            if 'Cartesian' not in structure.Ic[i]:
                count_atom += 1
                if float(structure.Ic[i][direction]) <= dire:
                    structure.Ic[i].append(' '.join(tempt))
                    count_fix_atom += 1
                else:
                    structure.Ic[i].append(' '.join(temp))
                    count_free_atom += 1
        if target == 2:
            if 'Cartesian' not in structure.Ic[i]:
                count_atom += 1
                if float(structure.Ic[i][direction]) <= dire:
                    del structure.Ic[i][-3:]
                    structure.Ic[i].append(' '.join(tempt))
                    count_fix_atom += 1
                else:
                    count_free_atom += 1

    filename = input('Enter the name of the output file (Enter 0 to use the original name): ')
    if filename == '0':
        filename = 'structure'

    with open('{}.as'.format(filename), 'w') as f:
        for i in range(len(structure.Ic)):
            f.write(' '.join(structure.Ic[i])+'\n')

    with open('atom_count.txt', 'w') as f:
        f.write('Total number of atoms: {}\n'.format(count_atom))
        f.write('Number of fixed atoms: {}\n'.format(count_fix_atom))
        f.write('Number of free atoms: {}\n'.format(count_free_atom))
elif type == '0':
    # Select the element of the fixed atom
    print('The atom type of the structure: ')
    print(structure.atom_type)
    elements = []
    element = '0'
    while True:
        element = input('Please select the element of the fixed atom (enter -1 to end): ')
        element = element.split()
        for lenele in range(len(element)):
            elements.append(element[lenele].capitalize())
        if element == ['-1']:
            break

    # Judge whether the file has been fixed
    if 'Fix' in structure.I[structure.target]:
        target = 2
    else:
        target = 1

    count_free_atom = 0
    count_fix_atom = 0
    count_atom = 0

    if target == 1:
        for i in range(len(structure.Ic)):
            if 'Cartesian' in structure.Ic[i]:
                structure.Ic[i].append('Fix_x Fix_y Fix_z')
            if i > structure.target:
                count_atom += 1
                if structure.Ic[i][0] in elements:
                    structure.Ic[i].append('T T T')
                    count_fix_atom += 1
                else:
                    structure.Ic[i].append('F F F')
                    count_free_atom += 1
    elif target == 2:
        for i in range(len(structure.Ic)):
            if i > structure.target:
                count_atom += 1
                if structure.Ic[i][0] in elements:
                    del structure.Ic[i][-3:]
                    structure.Ic[i].append('T T T')
                    count_fix_atom += 1
                else:
                    count_free_atom += 1
    
    if count_atom != structure.n_atoms:
        raise ValueError('The number of atoms is not correct!')


    filename = input('Enter the name of the output file (Enter 0 to use the original name): ')
    if filename == '0':
        filename = 'structure'
    structure.write_as_file(path = dir, file = '{}.as'.format(filename))

    with open('atom_count.txt', 'w') as f:
        f.write('Total number of atoms: {}\n'.format(count_atom))
        f.write('Number of fixed atoms: {}\n'.format(count_fix_atom))
        f.write('Number of free atoms: {}\n'.format(count_free_atom))