#!/usr/bin/env python
import os
import sys
current_path = os.path.dirname(os.path.abspath(__file__))
parent_path = os.path.dirname(current_path)
sys.path.append(parent_path)
import aseFunction as af
from tkinter import filedialog
from tkinter import Tk

# %%
path_root = os.getcwd()
dirList = os.listdir(path_root)
for dirName in dirList:
    if "NEB" in dirName:
        print(dirName)
nebDir = input("Please input the name of the NEB directory: ")

# %%
# 读取初始磁矩
root = Tk()
root.withdraw()
mag_file = filedialog.askopenfilename(initialdir = path_root)

# %%
mag_structure = af.read_as(mag_file)
mag = mag_structure.get_initial_magnetic_moments()

# %%
# 读取结构文件
dir = path_root + "\\" + nebDir + "\\DS-PAW-NEB"
for dirName in os.listdir(dir):
    if os.path.isdir(dir + "\\" + dirName):
        count = 0
        subDir = dir + "\\" + dirName
        for fileName in os.listdir(subDir):
            if fileName.endswith("as"):
                structureFile = subDir + "\\" + fileName
                count += 1
                if count > 1:
                    print(dirName)
                    raise ValueError("More than one structure file in the directory")
        structure = af.read_as(file = structureFile)
        structure.set_initial_magnetic_moments(mag)
        # change the name of the structure file
        os.rename(structureFile, subDir + "\\structure_old.as")
        af.write_as_file(path = subDir, structure = structure, filename = "structure" + dirName + ".as")
        print("Already written: " + "structure" + dirName + ".as")

# %%



