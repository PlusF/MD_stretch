import os
from Atom import Atom
from Cell import Cell
from CaseData import CaseData
from tqdm import tqdm
from constant import *


class Driver:
    def __init__(self, casedata: CaseData):
        self.casedata_ = casedata
        self.cell_ = None
        self.n_atoms_ = 0
        self.current_step_ = 0
        self.counter_ = 0  # for checking relaxation
        self.load()

    def load(self):
        atom_list = []

        # restartのとき
        if self.casedata_.restart_:
            if not os.path.exists(self.casedata_.restart_file_):
                raise ValueError("No restart file found. You must set the value 'restart' false.")
            with open(self.casedata_.restart_file_, 'r') as f:
                for i, line in enumerate(f):
                    line.strip('\n')
                    if i == 0:
                        self.current_step_ = int(line.split()[-1])
                        self.casedata_.n_loop_ += self.current_step_
                        continue
                    kind, x, y, z, vx, vy, vz = line.split()
                    a = Atom(kind, [x, y, z], [vx, vy, vz])
                    atom_list.append(a)
                    self.n_atoms_ += 1
            with open(self.casedata_.restart_file_, 'w') as f:
                pass
            
        else:
            with open(self.casedata_.in_file_, 'r') as f:
                for line in f:
                    line.strip('\n')
                    kind, x, y, z, vx, vy, vz = line.split()
                    a = Atom(kind, [x, y, z], [vx, vy, vz])
                    atom_list.append(a)
                    self.n_atoms_ += 1

        self.cell_ = Cell(self.casedata_, atom_list)

    def run(self):
        for i in tqdm(range(self.casedata_.n_loop_ - self.current_step_)):
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
        self.migrate()
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
        self.migrate()
        self.calc_force_and_up()
        self.update_velocity_half_and_calc_uk()
        if self.casedata_.relax_:
            self.relax()
        if self.casedata_.need_stretch_:
            self.stretch()

    def calc_force(self):
        self.cell_.clear_force()
        self.cell_.calc_force()
        # periodic B.C.
        for offset in self.casedata_.offsets_:
            self.cell_.calc_force_with_surrounding(offset)

    def calc_force_and_up(self):
        self.cell_.clear_force()
        self.cell_.calc_force_and_up()
        # periodic B.C.
        for offset in self.casedata_.offsets_:
            self.cell_.calc_force_and_up_with_surrounding(offset)

    def update_velocity_half(self):
        self.cell_.update_velocity_half()

    def update_velocity_half_and_calc_uk(self):
        self.cell_.update_velocity_half_and_calc_uk()

    def update_position(self):
        self.cell_.update_position()

    def migrate(self):
        self.cell_.migrate()

    def relax(self):
        self.cell_.relax()

    def stretch(self):
        self.counter_ += 1
        if self.counter_ < 50:  # ある程度時間が経たないと速度が出ず、緩和したかどうか判定できない
            return
        if self.cell_.is_relaxed():
            new_box_size = self.casedata_.box_size_ * self.casedata_.stretch_eps_
            # print('%05d (%03d) %.03e -> %.03e' % (self.current_step_, self.counter_, self.casedata_.box_size_[0], new_box_size[0]))
            self.casedata_.set_box_size(new_box_size)
            self.cell_.stretch()
            self.counter_ = 0

    def write_trajectory(self):
        trajectory_str = ''
        trajectory_str += f'{self.n_atoms_}\n'
        trajectory_str += f'#step {self.current_step_}\n'
        trajectory_str += self.cell_.get_trajectory()

        with open(self.casedata_.out_file_traj_, 'a') as f:
            f.write(trajectory_str)

    def write_energy(self):
        up = uk = 0
        up += self.cell_.up_
        uk += self.cell_.uk_
        temperature = uk / (1.5 * self.n_atoms_ * kB)

        with open(self.casedata_.out_file_energy_, 'a') as f:
            f.write(f'{self.current_step_} {up} {uk} {temperature}\n')

    def write_cell_state(self):
        cell_str = f'{self.current_step_} '
        d1 = self.cell_.d1_ * 1e10
        d2 = self.cell_.d2_ * 1e10
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
        restart_str += self.cell_.get_restart()

        with open(self.casedata_.restart_file_, 'a') as f:
            f.write(restart_str)


def test_periodic(dr: Driver):
    dr.run()


def main():
    casedata = CaseData('./data/case0.json')
    dr = Driver(casedata)
    test_periodic(dr)


if __name__ == '__main__':
    main()
