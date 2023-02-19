import numpy as np
from scipy import constants as const
import scipy.integrate as integrate
import scipy as sp

def molar_mass(element):
    """
    Calculation of molar mass
    - param element: Element
    - return: Molar mass
    """
    mass_list = {'H': 1.00794, 'He': 4.002602, 'Li': 6.941, 'Be': 9.012182, 'B': 10.811, 
                'C': 12.0107, 'N': 14.0067, 'O': 15.9994, 'F': 18.9984032, 'Ne': 20.1797, 
                'Na': 22.98976928, 'Mg': 24.3050, 'Al': 26.9815386, 'Si': 28.0855, 
                'P': 30.973762, 'S': 32.065, 'Cl': 35.453, 'Ar': 39.948, 'K': 39.0983, 
                'Ca': 40.078, 'Sc': 44.955912, 'Ti': 47.867, 'V': 50.9415, 'Cr': 51.9961, 
                'Mn': 54.938045, 'Fe': 55.845, 'Co': 58.933195, 'Ni': 58.6934, 'Cu': 63.546, 
                'Zn': 65.38, 'Ga': 69.723, 'Ge': 72.64, 'As': 74.92160, 'Se': 78.96, 
                'Br': 79.904, 'Kr': 83.798, 'Rb': 85.4678, 'Sr': 87.62, 'Y': 88.90585, 
                'Zr': 91.224, 'Nb': 92.90638, 'Mo': 95.96, 'Tc': 98, 'Ru': 101.07, 
                'Rh': 102.90550, 'Pd': 106.42, 'Ag': 107.8682, 'Cd': 112.411, 'In': 114.818, 
                'Sn': 118.710, 'Sb': 121.760, 'Te': 127.60, 'I': 126.90447, 'Xe': 131.293, 
                'Cs': 132.9054519, 'Ba': 137.327, 'La': 138.90547, 'Ce': 140.116, 
                'Pr': 140.90765, 'Nd': 144.242, 'Pm': 145, 'Sm': 150.36, 'Eu': 151.964,
                'Gd': 157.25, 'Tb': 158.92535, 'Dy': 162.500, 'Ho': 164.93032, 'Er': 167.259,
                'Tm': 168.93421, 'Yb': 173.054, 'Lu': 174.9668, 'Hf': 178.49, 'Ta': 180.94788,
                'W': 183.84, 'Re': 186.207, 'Os': 190.23, 'Ir': 192.217, 'Pt': 195.084,
                'Au': 196.966569, 'Hg': 200.59, 'Tl': 204.3833, 'Pb': 207.2, 'Bi': 208.98040,
                'Po': 209, 'At': 210, 'Rn': 222, 'Fr': 223, 'Ra': 226, 'Ac': 227, 'Th': 232.03806,
                'Pa': 231.03588, 'U': 238.02891, 'Np': 237, 'Pu': 244, 'Am': 243, 'Cm': 247,
                'Bk': 247, 'Cf': 251, 'Es': 252, 'Fm': 257, 'Md': 258, 'No': 259, 'Lr': 262,
                'Rf': 261, 'Db': 262, 'Sg': 266, 'Bh': 264, 'Hs': 277, 'Mt': 268, 'Ds': 281,
                'Rg': 272, 'Cn': 285, 'Nh': 284, 'Fl': 289, 'Mc': 288, 'Lv': 293, 'Ts': 294,
                'Og': 294}
    return mass_list[element]

# caculate the moment of inertia of  molecule
# 计算转动惯量
def moment_of_inertia(molecule):
    """
    Calculation of moment of inertia
    - param molecule: Molecule
    - return: Moment of inertia
    """
    # read the molecule file
    with open(molecule, 'r') as f:
        lines = f.readlines()
    # get the number of atoms
    natoms = int(lines[0])
    # get the coordinates of atoms
    coordinates = np.zeros((natoms, 3))
    for i in range(natoms):
        coordinates[i, :] = np.array(lines[i + 2].split()[1:4], dtype=float)
    # get the element of atoms
    elements = []
    for i in range(natoms):
        elements.append(lines[i + 2].split()[0])
    # get the mass of atoms
    mass = np.zeros(natoms)
    for i in range(natoms):
        mass[i] = molar_mass(elements[i])
    # calculate the moment of inertia
    I = 0
    for i in range(natoms):
        for j in range(natoms):
            if i == j:
                continue
            else:
                I += mass[i] * mass[j] * np.linalg.norm(coordinates[i, :] - coordinates[j, :]) ** 2
    return I

# caculate the Spin multiplicity
# 计算自旋多重度
def spin_multiplicity(molecule):
    """
    Calculation of spin multiplicity
    - param molecule: Molecule
    - return: Spin multiplicity
    """
    # read the molecule file
    with open(molecule, 'r') as f:
        lines = f.readlines()
    # get the number of atoms
    natoms = int(lines[0])
    # get the coordinates of atoms
    coordinates = np.zeros((natoms, 3))
    for i in range(natoms):
        coordinates[i, :] = np.array(lines[i + 2].split()[1:4], dtype=float)
    # get the element of atoms
    elements = []
    for i in range(natoms):
        elements.append(lines[i + 2].split()[0])
    # get the mass of atoms
    mass = np.zeros(natoms)
    for i in range(natoms):
        mass[i] = molar_mass(elements[i])
    # calculate the moment of inertia
    I = 0
    for i in range(natoms):
        for j in range(natoms):
            if i == j:
                continue
            else:
                I += mass[i] * mass[j] * np.linalg.norm(coordinates[i, :] - coordinates[j, :]) ** 2
    # calculate the spin multiplicity
    S = 1 # spin of electron S是自旋角动量
    while True:
        if (2 * S + 1) ** 2 * I < 2 * const.physical_constants['Bohr radius'][0] ** 2:
            S += 1
        else:
            break
    return 2 * S + 1


# Calculation of translational partition function
def partition_function_translational(T, V, N, m):
    """
    Calculation of translational partition function
    - param T: Temperature in K
    - param V: Volume in m^3
    - param N: Number of particles
    - return: Translational partition function
    """
    k = const.k  # Boltzmann constant: 1.38064852e-23 J/K
    h = const.h  # Planck constant: 6.62607004e-34 J*s
    qt = (2 * np.pi * m * k * T / h ** 2) ** (3 * N / 2) * V ** N
    return qt

# Calculation of rotational partition function
def partition_function_rotational(T, N, m, I):
    """
    Calculation of rotational partition function
    - param T: Temperature in K
    - param N: Number of particles
    - param m: Mass of particle in kg
    - param I: Moment of inertia in kg*m^2
    - return: Rotational partition function
    """
    k = const.k  # Boltzmann constant: 1.38064852e-23 J/K
    qr = (8 * np.pi ** 2 * k * T / I) ** (N / 2) * (np.pi * I / m / k / T) ** (3 * N / 2)
    return qr

# Calculation of vibrational partition function
def partition_function_vibrational(T, N, m, w):
    """
    Calculation of vibrational partition function
    - param T: Temperature in K
    - param N: Number of particles
    - param m: Mass of particle in kg
    - param w: Frequency in Hz
    - return: Vibrational partition function
    """
    k = const.k  # Boltzmann constant: 1.38064852e-23 J/K
    qv = (2 * np.pi * m * k * T / w) ** (3 * N / 2)
    return qv

# Calculation of partition function
def partition_function(T, V, N, m, I, w):
    """
    Calculation of partition function
    - param T: Temperature in K
    - param V: Volume in m^3
    - param N: Number of particles
    - param m: Mass of particle in kg
    - param I: Moment of inertia in kg*m^2
    - param w: Frequency in Hz
    - return: Partition function
    """
    qt = partition_function_translational(T, V, N, m)
    qr = partition_function_rotational(T, N, m, I)
    qv = partition_function_vibrational(T, N, m, w)
    q = qt * qr * qv
    return q

# Calculation of thermodynamic energy U
def thermodynamic_energy(T, V, N, m, I, w):
    """
    Calculation of thermodynamic energy U
    - param T: Temperature in K
    - param V: Volume in m^3
    - param N: Number of particles
    - param m: Mass of particle in kg
    - param I: Moment of inertia in kg*m^2
    - param w: Frequency in Hz
    - return: Thermodynamic energy U
    """
    k = const.k  # Boltzmann constant: 1.38064852e-23 J/K
    q = partition_function(T, V, N, m, I, w)
    u = k * T * np.log(q)
    return u