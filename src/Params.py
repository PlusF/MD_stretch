import numpy as np
import matplotlib.pyplot as plt


class Params:
    def __init__(self):
        self.dict_params_ = {
            'Al': {
                'mass': 1.67e-27 * 27,
                'epsilon': 0.27 * 1.6e-19,  # eV -> J
                'alpha': 1.16 * 1e10,  # Å-1 -> m-1
                'r0': 3.25 * 1e-10  # Å -> m
                }
        }
        self.init_params()

    def init_params(self):
        for key, value in self.dict_params_.items():
            self.dict_params_[key]['minus_alpha'] = -value['alpha']
            self.dict_params_[key]['minus_2_alpha'] = -2 * value['alpha']
            self.dict_params_[key]['minus_2_alpha_eps'] = -2 * value['alpha'] * value['epsilon']

    def Morse(self, r):
        # TODO: Al以外の原子のとき
        params = self.dict_params_['Al']
        force = params['minus_2_alpha_eps'] * (np.exp(params['minus_2_alpha'] * (r - params['r0'])) - np.exp(params['minus_alpha'] * (r - params['r0'])))
        return force / params['mass']

    def Morse_calc_up(self, r):
        # TODO: Al以外の原子のとき
        params = self.dict_params_['Al']
        force = params['minus_2_alpha_eps'] * (np.exp(params['minus_2_alpha'] * (r - params['r0'])) - np.exp(params['minus_alpha'] * (r - params['r0'])))
        phi = params['epsilon'] * (np.exp(params['minus_2_alpha'] * (r - params['r0'])) - 2 * np.exp(params['minus_alpha'] * (r - params['r0'])))
        return force / params['mass'], phi


def test_Params():
    params = Params()
    r = np.linspace(2e-10, 7e-10, 100)
    _, u = params.Morse_calc_up(r)
    plt.plot(r, u, color='k')
    min_u = u.min()
    min_r = r[np.where(u==min_u)][0]
    plt.plot(min_r, min_u, marker='o', color='r')
    plt.text(min_r * 1.3, min_u, 'U(%.2e) = %.2e' % (min_r, min_u))
    plt.show()


def main():
    test_Params()


if __name__ == '__main__':
    main()