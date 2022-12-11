import inspect
from Params import Params
from CaseData import CaseData
from Atom import Atom
import numpy as np


class Cell:
    def __init__(self, casedata: CaseData, atom_list):
        self.casedata_ = casedata
        self.d1_ = np.array([0, 0, 0])
        self.d2_ = np.array(casedata.box_size_)
        self.atom_list_ = atom_list
        self.n_atoms_ = len(self.atom_list_)
        self.up_ = self.uk_ = 0

        self.params_ = Params()

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

    def migrate(self):
        for atom in self.atom_list_:
            for i in range(3):
                if not self.casedata_.periodic_[i]:
                    continue
                if atom.r_[i] < self.d1_[i] - self.casedata_.margin_[i]:
                    atom.r_[i] += self.d2_[i] - self.d1_[i]
                elif atom.r_[i] >= self.d2_[i] + self.casedata_.margin_[i]:
                    atom.r_[i] -= self.d2_[i] - self.d1_[i]

    def relax(self):
        for atom in self.atom_list_:
            if np.dot(atom.v_, atom.a_) < 0:
                atom.v_ = np.zeros(3).astype(float)

    def is_relaxed(self):
        avg_vel = 0
        for atom in self.atom_list_:
            avg_vel += np.linalg.norm(atom.v_)
        avg_vel /= self.n_atoms_

        if avg_vel < 1:  # 閾値は適切か？
            return True
        else:
            return False

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


def test_clear_force(cell: Cell):
    cell.calc_force()
    for atom in cell.atom_list_:
        if np.linalg.norm(atom.a_) == 0:
            print(f'Failed in {inspect.currentframe().f_code.co_name}. There are some isolated atoms.')
            return
    cell.clear_force()
    for atom in cell.atom_list_:
        if np.linalg.norm(atom.a_) != 0:
            print(f'Failed in {inspect.currentframe().f_code.co_name}. clear_force() is not properly working.')
            return
    print(f'Succeed in {inspect.currentframe().f_code.co_name}')


def test_calc_force_and_up(cell: Cell):
    cell.calc_force_and_up()


def test_calc_force_and_up_with_surrounding():
    casedata = CaseData('./data/case0.json')
    atom_list = [
        Atom('Al', [0, 0, 0], [-100, 0, 0]),
        Atom('Al', [0, 8e-10, 0], [0, 0, 0]),
        Atom('Al', [0, 0, 4e-10], [0, 0, 0]),
    ]
    cell = Cell(casedata, atom_list)


def update_velocity_half_and_calc_uk(cell: Cell):
    pass


def test_update_position(cell: Cell):
    pass


def test_migrate():
    print('Testing migration: ', end='')
    casedata = CaseData('./data/case0.json')
    atom_list = [
        Atom('Al', [0, 0, 0], [-1000, 0, 0])
    ]
    cell = Cell(casedata, atom_list)
    for _ in range(10):
        print(cell.atom_list_[0].r_[0], end=', ')
        cell.update_position()
        cell.migrate()
    print('...')


def main():
    casedata = CaseData('./data/case0.json')
    atom_list = [
        Atom('Al', [0, 0, 0], [-100, 0, 0]),
        Atom('Al', [0, 4e-10, 0], [0, 0, 0]),
        Atom('Al', [0, 0, 4e-10], [0, 0, 0]),
    ]
    cell = Cell(casedata, atom_list)
    test_clear_force(cell)
    test_calc_force_and_up_with_surrounding()
    test_migrate()


if __name__ == '__main__':
    main()
