import numpy as np
import os
import json

def tran(fix):
    if fix == 'T':
        return 'F'
    else:
        return 'T'

# get the directory where the py file is located
# dir = os.path.dirname(os.path.realpath(__file__))
# os.chdir(dir)
dir = os.getcwd()

# load DS-PAW file
with open('frequency.json', 'r') as f:
    freq = json.load(f)
frelen = len(freq['FrequencyInfo'])
with open('structure.as', 'r') as f:
    I = f.readlines()

# Determine whether the file folder exists
if not os.path.exists('vaspfiles'):

    # make a new directory in the current directory
    os.mkdir(dir + '/vaspfiles')

    # set the work path as the new directory
    os.chdir(dir + '/vaspfiles')
else:
    os.chdir(dir + '/vaspfiles')


# make a new file
with open('OUTCAR', 'w') as f:
    f.write(''' Eigenvectors and eigenvalues of the dynamical matrix
 ----------------------------------------------------
 

''')
    for i in range(frelen):
        # 零点能单位换算：1(THz) = 6.28(2PiTHz) = 33.356 cm-1 = 4.135meV
        a = freq['FrequencyInfo'][i]['eigenvalues']['frequency']
        p = np.around(a * 2 * np.pi, decimals=6)
        cm = np.around(a * 33.356, decimals=6)
        mev = np.around(a * 4.135, decimals=6)
        a = np.around(a, decimals=6)

        if freq['FrequencyInfo'][i]['eigenvalues']['frequencyTpye'] == 'f':
            f.write('   '+str(i+1)+' '+freq['FrequencyInfo'][i]['eigenvalues']['frequencyTpye']+'  '+'='+'%+12s'
                %str(a.tolist())+' '+freq['FrequencyInfo'][i]['eigenvalues']['unit']+'%+13s'%str(p.tolist())
                +' 2PiTHz'+'%+13s'%str(cm.tolist())+' cm-1'+'%+13s'%str(mev.tolist())+' meV'+'\n'+'\n')
        else:
            f.write('   '+str(i+1)+' '+freq['FrequencyInfo'][i]['eigenvalues']['frequencyTpye']+'='+'%+12s'
                %str(a.tolist())+' '+freq['FrequencyInfo'][i]['eigenvalues']['unit']+'%+13s'%str(p.tolist())
                +' 2PiTHz'+'%+13s'%str(cm.tolist())+' cm-1'+'%+13s'%str(mev.tolist())+' meV'+'\n'+'\n')

with open('CONTCAR', 'w') as f:
    f.write('POSCAR generated by tran_as_POSCAR.py\n')
    f.write('1.0\n')
    len_I = len(I)

    for i in range(len_I):
        if 'atoms' in I[i]:
            a_tot = float(I[i+1].split()[0])
        if 'Lattice' in I[i]:
            al = I[i+1].split()
            al = ['%.16f' % float(x) for x in al[0:3]]
            bl = I[i+2].split()
            bl = ['%.16f' % float(x) for x in bl[0:3]]
            cl = I[i+3].split()
            cl = ['%.16f' % float(x) for x in cl[0:3]]
        if 'Cartesian' in I[i]:
            if 'Fix' in I[i]:
                # It = I[i].split()
                pos = {'atom': [], 'cart': [], 'Fix': []}
                for j in range(i+1, i+1+int(a_tot)):
                    Ic = I[j].split()
                    pos['atom'].append(Ic[0])
                    pos['cart'].append(['%.16f' % float(x) for x in Ic[1:4]])
                    pos['Fix'].append([tran(x) for x in Ic[4:7]])
            else:
                pos = {'atom': [], 'cart': []}
                for j in range(i+1, i+1+int(a_tot)):
                    Ic = I[j].split()
                    pos['atom'].append(Ic[0])
                    pos['cart'].append(['%.16f' % float(x) for x in Ic[1:4]])

    f.write('%+22s' % al[0]+'%+22s' % al[1]+'%+22s' % al[2]+'\n')
    f.write('%+22s' % bl[0]+'%+22s' % bl[1]+'%+22s' % bl[2]+'\n')
    f.write('%+22s' % cl[0]+'%+22s' % cl[1]+'%+22s' % cl[2]+'\n')

    # Count the number of atoms in the order of 'pos' list
    atom_num = {}
    list_atom = []
    list_atom_value = []
    # count = -1
    count1 = -1
    for i in pos['atom']:
        # count += 1
        if i not in list_atom:
            crit = i
            count1 += 1
            list_atom.append(i)
            list_atom_value.append(1)
        elif i != crit: 
            crit = i
            count1 += 1
            list_atom.append(i)
            list_atom_value.append(1)
        else:
            list_atom_value[count1] += 1

    for i in list_atom:
        f.write('%+5s' % i)
    f.write('\n')
    for i in list_atom_value:
        f.write('%+6s' % i)
    f.write('\n')
    f.write('''Cartesian
''')
    if 'Fix' in pos:
        for i in range(len(pos['atom'])):
            f.write('%+20s' % pos['cart'][i][0]+'%+20s' % pos['cart'][i][1]+'%+20s' % pos['cart'][i]
                    [2]+'%+4s' % pos['Fix'][i][0]+'%+4s' % pos['Fix'][i][1]+'%+4s' % pos['Fix'][i][2]+'\n')
    else:
        for i in range(len(pos['atom'])):
            f.write('%+20s' % pos['cart'][i][0]+'%+20s' % pos['cart'][i][1]+'%+20s' % pos['cart'][i][2]+'\n')