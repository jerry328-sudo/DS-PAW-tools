# %%
import sys
sys.path.append("e:\\code\\DS-PAW\\DS-PAW-tools\\code")
import functions.function_normal as fn
import matplotlib.pyplot as plt
import os

path = os.getcwd()
interpolation_kind = input('''Enter interpolation kind: 
1. 'zero'
2. 'slinear'
3. 'quadratic'
4. 'cubic'(when the number of data points is more than 4, means if data_clean is true, you should not choose this option)
5. enter function interp1d's interp_kinds
''')
if interpolation_kind == '1':
    interpolation_kind = 'zero'
elif interpolation_kind == '2':
    interpolation_kind = 'slinear'
elif interpolation_kind == '3':
    interpolation_kind = 'quadratic'
elif interpolation_kind == '4':
    interpolation_kind = 'cubic'
else:
    pass
data_clean = input("Whether to clean the data (only keep the highest point) (y/n): ")
if data_clean == 'y':
    data_clean = True
else:
    data_clean = False

data_save = input("Whether to save the data (y/n): ")
if data_save == 'y':
    data_save = True
else:
    data_save = False
figure_plot = input("Whether to plot the figure (y/n): ")
if figure_plot == 'y':
    figure_plot = True
else:
    figure_plot = False

# %% MAIN CODE

fn.plot_neb_barrier(neb_json=path + "\\neb.json", interp_kinds=interpolation_kind, data_clean=data_clean, data_save=data_save)
plt.xlabel('Reaction Coordinate', fontsize=16)
plt.ylabel('Energy (eV)', fontsize=16)
plt.tick_params(labelsize=16)
# plt.savefig("neb_reaction_coordinate.png", img_format="png")
if figure_plot:
    plt.show()