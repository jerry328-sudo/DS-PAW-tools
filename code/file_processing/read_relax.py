import os
import sys
current_path = os.path.dirname(os.path.abspath(__file__))
parent_path = os.path.dirname(current_path)
sys.path.append(parent_path)
import functions.aseFunction as af
import json
import ase.io

structure = af.relaxLoad()

structure.write_as_file()

with open('energy.json','w',encoding='utf-8') as f:
    f.write(json.dumps(structure.energy,ensure_ascii=False,indent=4))

ase.io.write("relax.xyz", structure.RelaxTraj, format='xyz')

print("""relaxed structure: structure_new.as
energy: energy.json
""")
print(json.dumps(structure.energy,ensure_ascii=False,indent=4))

# press any key to exit
os.system('pause')
