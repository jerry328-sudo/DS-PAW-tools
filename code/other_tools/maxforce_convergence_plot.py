import matplotlib.pyplot as plt
from tkinter import filedialog
from tkinter import *
import os
import numpy as np

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

# chose the part of the trace to plot
p1 = int(input('Please enter the starting point of the trace to plot: \n'))
p2 = int(input('Please enter the end point of the trace to plot (enter "-1" to not set the end point): \n'))
if p2 == -1:
    p2 = len(trace)

while p1 < 0 or p2 > len(trace) or p1 > p2:
    print('The starting point and end point you entered are out of range, please re-enter!')
    p1 = int(input('Please enter the starting point of the trace to plot (enter "-1" to exit): \n'))
    if p1 == -1:
        exit()
    p2 = int(input('Please enter the end point of the trace to plot (enter "-1" to not set the end point): \n'))
if p2 == -1:
    p2 = len(trace)
if p1 != 0:
    p1 -= 1
trace = trace[p1:p2]

# Set the value of x-axis
x = np.arange(p1+1, p2+1, 1)

# Plot the trace
plt.plot(x, trace, 'r')
plt.xlabel('Iteration')
plt.ylabel('Maxforce')
plt.title('Maxforce Convergence Plot')
plt.show()

