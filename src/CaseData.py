import os
import json
import numpy as np


class CaseData:
    def __init__(self, filename):
        if not os.path.exists('./log'):
            os.mkdir('./log')

        with open(filename, 'r') as f:
            case_dict = json.load(f)

        self.in_file_ = case_dict["in_file"]
        self.restart_ = case_dict["restart"]
        self.restart_file_ = case_dict["restart_file"]
        self.out_file_traj_ = case_dict["out_file_traj"]
        self.out_file_energy_ = case_dict["out_file_energy"]
        self.out_file_cell_ = case_dict["out_file_cell"]

        self.dt_ = case_dict["dt"]  # time step
        self.cutoff_ = case_dict["cutoff"]  # cut off radius
        self.n_loop_ = case_dict["n_loop"]  # number of loop iterations
        self.interval_ = case_dict["interval"]  # interval between output
        self.periodic_ = case_dict["periodic"]  # B.C x, y, z
        self.n_cell_ = np.array(case_dict["n_cell"])

        self.margin_ = self.offsets_ = None
        if self.restart_:
            with open(self.out_file_cell_, 'r') as f:
                _, x, _, _, _, y, _, _, _, z, _, _, _ = f.readlines()[-1].split()  # cellファイルの最後の出力を読む
                box_size = np.array([x, y, z]).astype(float) * 1e-10
            self.set_box_size(box_size)
        else:
            # reset the output file
            with open(self.out_file_traj_, 'w'):
                pass
            with open(self.out_file_energy_, 'w') as f:
                f.write('#step up uk temperature\n')
            with open(self.out_file_cell_, 'w') as f:
                f.write('#$LABELS step a_x a_y a_z b_x b_y b_z c_x c_y c_z o_x o_y o_z\n')

            self.set_box_size(np.array(case_dict["box_size"]).astype(float))

        self.relax_ = case_dict["relax"]  # structure relaxation
        self.stretch_eps_ = np.array(case_dict["stretch_eps"]).astype(float)  # distortion of stretch per step
        self.need_stretch_ = True if np.linalg.norm(self.stretch_eps_) > 0 else False
        self.stretch_eps_ += np.array([1, 1, 1]).astype(float)

    def set_box_size(self, box_size):
        self.box_size_ = np.array(box_size)
        self.margin_ = self.box_size_ / self.n_cell_ * 0.05
        # periodic B.C.を計算する際に適用するオフセット
        self.set_offsets()

    def set_offsets(self):
        d = [-1, 0, 1]
        combi = [[i, j, k] for i in d for j in d for k in d]
        combi.remove([0, 0, 0])
        combi = np.array(combi).astype(float)

        # 指定された境界条件に応じて、力を計算する際のオフセットを設定 (periodicがFalseの軸が0であるようなオフセット)
        self.offsets_ = combi[np.sum(self.periodic_ | ((1 - np.abs(combi)).astype(bool)), axis=1) == 3]

        self.offsets_ *= self.box_size_

    def check_cell_size_and_cutoff(self):
        # TODO: IMPLEMENT ME
        pass


def test_offsets(casedata: CaseData):
    pass


def main():
    casedata = CaseData('./data/case0.json')
    test_offsets(casedata)


if __name__ == '__main__':
    main()
