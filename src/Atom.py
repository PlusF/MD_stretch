import numpy as np


class Atom:
    def __init__(self, kind, r, v):
        self.kind_ = kind
        self.r_ = np.array(r).astype(np.double)
        self.v_ = np.array(v).astype(np.double)
        self.a_ = np.array([0, 0, 0]).astype(np.double)

    def __str__(self):
        return f'Atom({self.kind_}: {self.r_}, {self.v_})'

    def out_traj(self):
        return f'{" ".join([self.kind_, *(self.r_ * 1e10).astype(str), *self.v_.astype(str)])}\n'

    def out_restart(self):
        return f'{" ".join([self.kind_, *(self.r_).astype(str), *self.v_.astype(str)])}\n'