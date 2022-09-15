import matplotlib.pyplot as plt
from tkinter import filedialog
from tkinter import *
import os

print('The current working directory is: ' + '\n' + os.getcwd())
dir = input('Please enter the directory where the "Maxforce.log" file is located(press any key to continue): \n'
            'If the file is in the current directory, please enter "." \n')
if dir == '.':
    dir = os.getcwd()
else:
    root = Tk()
    root.withdraw()
    dir = filedialog.askdirectory()

os.chdir(dir)
print('The current working directory is: ' + '\n' + os.getcwd() + '\n')
# Determine whether the file folder exists
print('reading Maxforce.log ...', end=' ')
if not os.path.exists('Maxforce.log'):
    print('Maxforce.log does not exist!')
    exit()
else:
    print('done')
# Read the Maxforce.log file
Ic = []
with open('Maxforce.log', 'r') as f:
    I = f.readlines()
    for i in range(len(I)):
        Ic.append(I[i].split())

trace = []
for ic in Ic:
    trace.append(float(ic[2]))

plt.plot(trace)
plt.xlabel('Step')
plt.ylabel('Maxforce')
plt.show()
