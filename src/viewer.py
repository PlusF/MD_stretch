import numpy as np
import matplotlib.pyplot as plt

def load_energy(filename):
    step_list = []
    up_list = []
    uk_list = []
    temperature_list = []
    with open(filename, 'r') as f:
        for line in f:
            line.strip('\n')
            step, up, uk, temperature = line.split()
            if step == '#step':
                continue
            step_list.append(int(step))
            up_list.append(float(up))
            uk_list.append(float(uk))
            temperature_list.append(float(temperature))

    return np.array(step_list), np.array(up_list), np.array(uk_list), np.array(temperature_list)


def calc_stress_strain(up_list, lx_list):
    stride = 50
    sigma_list = []
    epslion_list = []
    for i, lx in enumerate(lx_list):
        if i < stride or i >= len(up_list) - stride:
            continue
        sigma = (up_list[i + stride] - up_list[i - stride]) / (2 * stride)
        epslion = (lx - lx_list[0]) / lx_list[0]
        sigma_list.append(sigma)
        epslion_list.append(epslion)

    return sigma_list, epslion_list


def load_cell(filename):
    box_sizes = []
    with open(filename, 'r') as f:
        for line in f:
            line.strip('\n')
            if line[0] == '#':
                continue
            _, x, _, _, _, y, _, _, _, z, _, _, _ = line.split()
            box_size = [x, y, z]
            box_sizes.append(box_size)

    return np.array(box_sizes).astype(float) * 1e-10


def main():
    step_list, up_list, uk_list, temperature_list = load_energy('./log/energy_8.out')

    box_sizes = load_cell('./log/cell_8.out')
    lx_list = box_sizes[:, 0]

    # sigma_list, epslion_list = calc_stress_strain(up_list, lx_list)

    fig, axes = plt.subplots(2, 2, figsize=(10, 7))
    axes[0, 0].plot(step_list, up_list, color='r', label='Up')
    axes[0, 0].plot(step_list, uk_list, color='b', label='Uk')
    axes[0, 0].plot(step_list, up_list + uk_list, color='k', label='Total')
    axes[0, 0].set_xlabel('Step')
    axes[0, 0].set_ylabel('Energy')
    axes[0, 0].legend()
    axes[0, 1].plot(step_list, temperature_list, color='k')
    axes[0, 1].set_xlabel('Step')
    axes[0, 1].set_ylabel('Temperature [K]')
    axes[1, 0].plot(step_list, lx_list, color='b', label='Length of x')
    axes[1, 0].set_xlabel('Step')
    axes[1, 0].set_ylabel('Lx')
    axes[1, 1].plot(lx_list, up_list, color='k')
    axes[1, 1].set_xlabel('Lx')
    axes[1, 1].set_ylabel('Potential Energy')
    plt.show()


if __name__ == '__main__':
    main()
    