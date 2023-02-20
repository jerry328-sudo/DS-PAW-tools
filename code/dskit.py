import os
import functions.function_normal as fn
# Create a UI for selecting functions
print('Please select the function you want to use:')
print('1. Thermodynamic correction')
print('2. Structure setting')
print('3. NEB tools')
print('9. others')

# Select the function
try:
    choice = int(
        input('Please enter the number of the function you want to use: '))
    if choice == 1:
        choice = int(input('''You have selected the function of thermodynamic correction.

Please enter the following function you want to use:
101. Convert DS frequency calculation file to vasp file
102. Thermal Corrections for Adsorbate 
103. Thermal Corrections for Gas \n'''))

    elif choice == 2:
        choice = int(input('''You have selected the function of structure setting.

Please enter the following function you want to use:
201. Set fixed atoms
202. Randomly dope atom
203. Set fixed atoms for NEB operation 
204. fix molecule\n'''))

    elif choice == 3:
        choice = int(input('''You have selected the function of NEB tools.

Please enter the following function you want to use:
301. Plot NEB reaction coordinate \n'''))
        
    elif choice == 9:
        choice = int(input('''You have selected the function of others.

Please enter the following function you want to use:
901. Plot the 'Max force' trace \n'''))
        
    else:
        print('Please enter the correct number!')
except ValueError:
    print('Please enter the correct number!')

dir = os.path.dirname(os.path.realpath(__file__))
# print(dir)
if choice == 101:
    # call the py file under the 'other tools' folder
    os.system('python ' + dir + "\\Thermodynamic_correction\\tran_fre_vasp.py")
elif choice == 102:
    os.system('python ' + dir + "\\Thermodynamic_correction\\Thermal_Corrections_for_Adsorbate.py")
elif choice == 103:
    os.system('python ' + dir + "\\Thermodynamic_correction\\Thermal_Corrections_for_Gas.py")
elif choice == 201:
    os.system('python ' + dir + "\\Structure_Settings\\set_fix_atom.py")
elif choice == 202:
    os.system('python ' + dir + "\\Structure_Settings\\Randomly_doped.py")
elif choice == 203:
    os.system('python ' + dir + "\\Structure_Settings\\set_fix_atom_NEB.py")
elif choice == 204:
    os.system('python ' + dir + "\\Structure_Settings\\fix_molecule.py")
elif choice == 301:
    os.system('python ' + dir + "\\NEB\\NEB_plot.py")
elif choice == 901:
    os.system('python ' + dir + "\\other_tools\\maxforce_convergence_plot.py")



