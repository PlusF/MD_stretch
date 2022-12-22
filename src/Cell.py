import inspect
from Params import Params
from CaseData import CaseData
from Atom import Atom
import numpy as np


class Cell:
    def __init__(self, casedata: CaseData):
        self.casedata_ = casedata
        self.d1_ = np.array([0, 0, 0])
        self.d2_ = np.array([0, 0, 0])
        self.cell_index_1d_ = 0
        self.cell_index_3d_ = np.array([0, 0, 0])
        self.atom_list_ = []
        self.n_atoms_ = 0
        self.up_ = self.uk_ = 0

        self.params_ = Params()

    def contains(self, a: Atom):
        if np.all((self.d1_ <= a.r_) | (a.r_ < self.d2_)):
            return True
        return False

    def add_atom(self, a: Atom):
        self.atom_list_.append(a)
        self.n_atoms_ += 1

    def remove_atom(self, a: Atom):
        self.atom_list_.remove(a)
        self.n_atoms_ -= 1

    def clear_force(self):
        for atom in self.atom_list_:
            atom.a_ = np.zeros(3).astype(float)

    def calc_force(self):
        for i in range(self.n_atoms_):
            for j in range(i + 1, self.n_atoms_):
                displacement = self.atom_list_[i].r_ - self.atom_list_[j].r_
                r = np.linalg.norm(displacement)
                if r > self.casedata_.cutoff_:
                    continue
                acc = self.params_.Morse(r) * (displacement / r)
                self.atom_list_[i].a_ -= acc
                self.atom_list_[j].a_ += acc

    def calc_force_with_surrounding(self, offset):
        for i in range(self.n_atoms_):
            for j in range(self.n_atoms_):
                displacement = self.atom_list_[i].r_ - (self.atom_list_[j].r_ + offset)
                r = np.linalg.norm(displacement)
                if r > self.casedata_.cutoff_:
                    continue
                acc = self.params_.Morse(r) * (displacement / r)
                self.atom_list_[i].a_ -= acc

    def calc_force_and_up(self):
        self.up_ = 0
        for i in range(self.n_atoms_):
            for j in range(i + 1, self.n_atoms_):
                displacement = self.atom_list_[i].r_ - self.atom_list_[j].r_
                r = np.linalg.norm(displacement)
                if r > self.casedata_.cutoff_:
                    continue
                acc, phi = self.params_.Morse_calc_up(r)
                acc *= displacement / r
                self.atom_list_[i].a_ -= acc
                self.atom_list_[j].a_ += acc
                self.up_ += phi

    def calc_force_and_up_with_surrounding(self, offset):
        for i in range(self.n_atoms_):
            for j in range(self.n_atoms_):
                displacement = self.atom_list_[i].r_ - (self.atom_list_[j].r_ + offset)
                r = np.linalg.norm(displacement)
                if r > self.casedata_.cutoff_:
                    continue
                acc, phi = self.params_.Morse_calc_up(r)
                acc *= displacement / r
                self.atom_list_[i].a_ -= acc
                self.up_ += phi * 0.5

    def update_velocity_half(self):
        for atom in self.atom_list_:
            atom.v_ += atom.a_ * self.casedata_.dt_ / 2

    def update_velocity_half_and_calc_uk(self):
        self.uk_ = 0
        for atom in self.atom_list_:
            atom.v_ += atom.a_ * self.casedata_.dt_ / 2
            self.uk_ += 0.5 * self.params_.dict_params_['Al']['mass'] * np.dot(atom.v_, atom.v_)

    def update_position(self):
        for atom in self.atom_list_:
            atom.r_ += atom.v_ * self.casedata_.dt_

    def migrate_and_reflect(self):
        for atom in self.atom_list_:
            for i in range(3):
                if self.casedata_.periodic_[i]:  # migrate
                    if atom.r_[i] < self.d1_[i] - self.casedata_.margin_[i]:
                        atom.r_[i] += self.d2_[i] - self.d1_[i]
                    elif atom.r_[i] >= self.d2_[i] + self.casedata_.margin_[i]:
                        atom.r_[i] -= self.d2_[i] - self.d1_[i]
                else:  # reflect
                    if atom.r_[i] < self.d1_[i]:
                        atom.v_[i] *= -1
                        atom.r_[i] = 2 * self.d1_[i] - atom.r_[i]
                    if atom.r_[i] >= self.d2_[i]:
                        atom.v_[i] *= -1
                        atom.r_[i] = 2 * self.d2_[i] - atom.r_[i]

    def relax(self):
        for atom in self.atom_list_:
            if np.dot(atom.v_, atom.a_) < 0:
                atom.v_ = np.zeros(3).astype(float)

    def get_sum_vel(self):
        return sum([np.linalg.norm(atom.v_) for atom in self.atom_list_])

    def update_box_size(self):
        self.d1_ = np.array([0, 0, 0])
        self.d2_ = np.array(self.casedata_.box_size_)

    def stretch(self):
        self.update_box_size()
        for atom in self.atom_list_:
            atom.r_ *= self.casedata_.stretch_eps_

    def get_trajectory(self):
        trajectory_str = ''
        for atom in self.atom_list_:
            trajectory_str += atom.out_traj()
        return trajectory_str

    def get_restart(self):
        restart_str = ''
        for atom in self.atom_list_:
            restart_str += atom.out_restart()
        return restart_str


def main():
    pass


if __name__ == '__main__':
    main()
