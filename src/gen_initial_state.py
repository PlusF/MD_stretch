import numpy as np
from Params import Params
from constant import *


class FCCMaker:
    def __init__(self, lattice_constant, num_lattice, kind, temperature):
        first = np.array([0, 0, 0])
        second = np.array([ 0.5,  0.5, 0])
        third = np.array([ 0.5, 0,  0.5])
        fourth = np.array([0,  0.5,  0.5])
        self.my_atoms_ = [first, second, third, fourth]
        self.lattice_constant_ = lattice_constant
        self.num_lattice_ = num_lattice
        self.kind_ = kind
        self.temperature_ = temperature

        self.params_ = Params()

    def make(self):
        self.lattice_ = []
        self.vel_ = []
        for ix in range(self.num_lattice_[0]):
            for iy in range(self.num_lattice_[1]):
                for iz in range(self.num_lattice_[2]):
                    org = np.array([ix, iy, iz]) * self.lattice_constant_
                    for atom in self.my_atoms_:
                        self.lattice_.append(org + atom * self.lattice_constant_)
                        self.vel_.append(self.get_vel_from_temp())

    def to_xyz(self, filename):
        with open(filename, 'w') as f:
            for atom, vel in zip(self.lattice_, self.vel_):
                f.write(f'{" ".join([self.kind_, *atom.astype(str), *vel.astype(str)])}\n')

    def get_vel_from_temp(self):
        m = self.params_.dict_params_[self.kind_]['mass']
        rand = np.array([(sum([np.random.rand() for _ in range(12)]) - 6) for _ in range(3)])
        return np.sqrt(kB * self.temperature_ / m) * rand


def main():
    lattice_constant = 4.05e-10  # m
    num_lattice = [2, 2, 2]  # x, y, z方向の格子の数 
    temperature = 100

    maker = FCCMaker(lattice_constant, num_lattice, 'Al', temperature)
    maker.make()
    maker.to_xyz('./data/Al_fcc_8.in')


if __name__ == '__main__':
    main()
