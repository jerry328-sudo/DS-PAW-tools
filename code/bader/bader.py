# %%
import json
from ase import Atom, Atoms
from ase.visualize import view
import numpy as np
from tkinter import filedialog
import tkinter as tk
from ase.visualize.ngl import NGLDisplay as ngl

# root = tk.Tk()
# root.withdraw()
# file_path = filedialog.askopenfilename(filetypes=[("json", "*.json")])
# print(file_path)


# %%
# with open(file_path, "r") as f:
with open("./bader.json", "r") as f:
    data = json.load(f)

# %%
cell=np.mat(data["AtomInfo"]["Lattice"]).reshape(3, 3)

# %%
atoms = Atoms(
            cell=np.array(data["AtomInfo"]["Lattice"]).reshape(3, 3).tolist(),
            pbc=True)
i=0
for a in data["AtomInfo"]["Atoms"]:
    p=a['Position']*cell
    p=p.tolist()
    temp = Atom(symbol=a['Element'], position=p[0], charge=data["BaderInfo"]["Charge"][i])
    atoms.append(temp)
    i+=1

# %%
view(atoms,block=True)
# ngl(atoms)
# input("Press Enter to continue...")