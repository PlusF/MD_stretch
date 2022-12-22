import os
from Atom import Atom
from Cell import Cell
from CaseData import CaseData
from tqdm import tqdm
from constant import *
import numpy as np


class Driver:
    def __init__(self, casedata: CaseData):
        self.casedata_ = casedata
        self.lc_ = self.casedata_.box_size_ / self.casedata_.n_cell_  # 1セルのサイズ
        self.cell_list_ = []
        self.n_atoms_ = 0
        self.current_step_ = 0
        self.counter_ = 0  # for checking relaxation

        self.init_cell()
        self.load_atom()

    def init_cell(self):
        for icx in range(-1, self.casedata_.n_cell_[0] + 1):
            for icy in range(-1, self.casedata_.n_cell_[1] + 1):
                for icz in range(-1, self.casedata_.n_cell_[2] + 1):
                    cell = Cell(self.casedata_)
                    cell.d1_ = np.array([self.lc_[0] * icx, self.lc_[1] * icy, self.lc_[2] * icz])
                    cell.d2_ = np.array([self.lc_[0] * (icx + 1), self.lc_[1] * (icy + 1), self.lc_[2] * (icz + 1)])
                    cell.cell_index_1d_ = self.get_1d_from_3d_cell_index([icx, icy, icz])
                    cell.cell_index_3d_ = np.array([icx, icy, icz])
                    self.cell_list_.append(cell)  # TODO: cell_indexをどうふるか（surroundingの扱い方）

    def migrate_inside(self, a: Atom):
        for i in range(3):
            if a.r_[i] < 0:
                a.r_[i] += self.casedata_.box_size_[i]
            elif a.r_[i] >= self.casedata_.box_size_[i]:
                a.r_[i] -= self.casedata_.box_size_[i]

    def get_cell_index_from_coord(self, coord):
        ic = (coord // self.lc_).astype(int)
        if np.all((np.array([0, 0, 0]) <= ic) & (ic < self.casedata_.n_cell_)):
            return self.get_1d_from_3d_cell_index(ic)
        else:
            raise ValueError(f'Atom out of range: {coord}')

    def get_1d_from_3d_cell_index(self, ic3):
        ic1 = ic3[0] * self.casedata_.n_cell_[1] * self.casedata_.n_cell_[2] + ic3[1] * self.casedata_.n_cell_[2] + ic3[2]
        return ic1

    def get_3d_from_1d_cell_index(self, ic1):
        icx = ic1 // (self.casedata_.n_cell_[1] * self.casedata_.n_cell_[2])
        icy = (ic1 - icx * (self.casedata_.n_cell_[1] * self.casedata_.n_cell_[2])) // self.casedata_.n_cell_[2]
        icz = ic1 % self.casedata_.n_cell_[2]
        return np.array([icx, icy, icz]).astype(int)

    def load_atom(self):
        # restartのとき
        if self.casedata_.restart_:
            if not os.path.exists(self.casedata_.restart_file_):
                raise ValueError("No restart file found. You must set the value 'restart' false.")
            with open(self.casedata_.restart_file_, 'r') as f:
                for i, line in enumerate(f):
                    line.strip('\n')
                    if i == 0:  # 前回のステップ数が記録されている
                        self.current_step_ = int(line.split()[-1])
                        self.casedata_.n_loop_ += self.current_step_
                        continue
                    kind, x, y, z, vx, vy, vz = line.split()
                    a = Atom(kind, [x, y, z], [vx, vy, vz])
                    self.migrate_inside(a)  # restartの場合、marginにある原子をmigrateしておく必要がある
                    ic = self.get_cell_index_from_coord(a.r_)
                    self.cell_list_[ic].add_atom(a)
                    self.n_atoms_ += 1
            
        else:
            with open(self.casedata_.in_file_, 'r') as f:
                for line in f:
                    line.strip('\n')
                    kind, x, y, z, vx, vy, vz = line.split()
                    a = Atom(kind, [x, y, z], [vx, vy, vz])
                    ic = self.get_cell_index_from_coord(a.r_)
                    self.cell_list_[ic].add_atom(a)
                    self.n_atoms_ += 1

        with open(self.casedata_.restart_file_, 'w'):
            pass

    def run(self):
        for _ in tqdm(range(self.casedata_.n_loop_ - self.current_step_)):
            if self.current_step_ % self.casedata_.interval_ == 0:
                self.do_step_with_output()
                self.write_trajectory()
                self.write_energy()
                self.write_cell_state()
            else:
                self.do_step()
            self.current_step_ += 1

        self.write_restart()

    def do_step(self):
        self.calc_force()
        self.update_velocity_half()
        self.update_position()
        self.migrate_and_reflect()
        self.calc_force()
        self.update_velocity_half()
        if self.casedata_.relax_:
            self.relax()
        if self.casedata_.need_stretch_:
            self.stretch()

    def do_step_with_output(self):
        self.calc_force()
        self.update_velocity_half()
        self.update_position()
        self.migrate_and_reflect()
        self.calc_force_and_up()
        self.update_velocity_half_and_calc_uk()
        if self.casedata_.relax_:
            self.relax()
        if self.casedata_.need_stretch_:
            self.stretch()

    def copy_to_surrounding(self):
        # TODO: IMPLEMENT ME
        pass

    def calc_force(self):
        for cell in self.cell_list_:
            cell.clear_force()
            cell.calc_force()
            # periodic B.C.
            for offset in self.casedata_.offsets_:  # TODO: 複数セルに対応できるようsurroundingを定義
                cell.calc_force_with_surrounding(offset)

    def calc_force_and_up(self):
        for cell in self.cell_list_:
            cell.clear_force()
            cell.calc_force_and_up()
            # periodic B.C.
            for offset in self.casedata_.offsets_:  # TODO: 複数セルに対応できるようsurroundingを定義
                    cell.calc_force_and_up_with_surrounding(offset)

    def update_velocity_half(self):
        for cell in self.cell_list_:
            cell.update_velocity_half()

    def update_velocity_half_and_calc_uk(self):
        for cell in self.cell_list_:
            cell.update_velocity_half_and_calc_uk()

    def update_position(self):
        for cell in self.cell_list_:
            cell.update_position()

    def migrate_and_reflect(self):
        for cell in self.cell_list_:
            cell.migrate_and_reflect()

    def relax(self):
        for cell in self.cell_list_:
            cell.relax()

    def stretch(self):
        self.counter_ += 1
        if self.counter_ < 50:  # ある程度時間が経たないと速度が出ず、緩和したかどうか判定できない
            return
        avg_vel = 0
        for cell in self.cell_list_:
            avg_vel += cell.get_sum_vel()
        avg_vel /= self.n_atoms_

        if avg_vel < 1:  # well relaxed
            new_box_size = self.casedata_.box_size_ * self.casedata_.stretch_eps_
            print('%05d (%03d) %.03e -> %.03e' % (self.current_step_, self.counter_, self.casedata_.box_size_[0], new_box_size[0]))
            self.casedata_.set_box_size(new_box_size)
            for cell in self.cell_list_:
                cell.stretch()
            self.counter_ = 0

    def write_trajectory(self):
        trajectory_str = ''
        trajectory_str += f'{self.n_atoms_}\n'
        trajectory_str += f'#step {self.current_step_}\n'
        for cell in self.cell_list_:
            trajectory_str += cell.get_trajectory()

        with open(self.casedata_.out_file_traj_, 'a') as f:
            f.write(trajectory_str)

    def write_energy(self):
        up = uk = 0
        for cell in self.cell_list_:
            up += cell.up_
            uk += cell.uk_
        temperature = uk / (1.5 * self.n_atoms_ * kB)

        with open(self.casedata_.out_file_energy_, 'a') as f:
            f.write(f'{self.current_step_} {up} {uk} {temperature}\n')

    def write_cell_state(self):
        cell_str = f'{self.current_step_} '
        d1 = np.array([0, 0, 0])
        d2 = self.casedata_.box_size_ * 1e10
        cell_str += f'{d2[0] - d1[0]} 0 0 '
        cell_str += f'0 {d2[1] - d1[1]} 0 '
        cell_str += f'0 0 {d2[2] - d1[2]} '
        cell_str += ' '.join(d1.astype(str))
        cell_str += '\n'
        with open(self.casedata_.out_file_cell_, 'a') as f:
            f.write(cell_str)

    def write_restart(self):
        restart_str = ''
        restart_str += f'#step {self.current_step_}\n'
        for cell in self.cell_list_:
            restart_str += cell.get_restart()

        with open(self.casedata_.restart_file_, 'a') as f:
            f.write(restart_str)


def test_periodic(dr: Driver):
    dr.run()


def main():
    casedata = CaseData('./data/case2.json')
    dr = Driver(casedata)
    # test_periodic(dr)
    print(dr.get_1d_from_3d_cell_index([0, 0, 0]))
    print(dr.get_1d_from_3d_cell_index([0, 0, 1]))
    print(dr.get_1d_from_3d_cell_index([0, 1, 1]))
    print(dr.get_1d_from_3d_cell_index([1, 1, 1]))
    for i in range(60):
        print(dr.get_3d_from_1d_cell_index(i))

if __name__ == '__main__':
    main()
