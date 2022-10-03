import os


# read the as file class

class Read_As_File:
    '''- a function to deal with the as file
    - Read_As_File(file = 'structure.as', path = os.getcwd())
    '''
    def __init__(self, file = 'structure.as', path = os.getcwd()):
        self.file = path+'\\'+file
        self.Ic = []
        self.I = []
        self.n_atoms = 0
        self.target = 6
        self.atom_type = {}
        self.__read_as_file()

    def __read_as_file(self):
        '''- read the as file'''
        with open(self.file, 'r') as f:
            self.I = f.readlines()
            for i in range(len(self.I)):
                self.Ic.append(self.I[i].split())

        # count the number of atoms
        self.n_atoms = int(self.Ic[1][0])
        self.target = len(self.I)
        for i in range(len(self.I)):
            if 'Cartesian' in self.Ic[i]:
                self.target = i
            if i > self.target:
                if self.Ic[i][0] not in self.atom_type.keys():
                    self.atom_type[self.Ic[i][0]] = 1
                else:
                    self.atom_type[self.Ic[i][0]] += 1
    
    def write_as_file(self, file = 'structure.as'):
        '''- write the as file'''
        with open(file, 'w') as f:
            for i in range(len(self.Ic)):
                f.write(' '.join(self.Ic[i])+'\n')