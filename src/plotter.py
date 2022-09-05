import sys

import numpy as np
import scipy.io as spio
from cycler import cycler
from matplotlib import pyplot as plt

# takes in filenames as sys.argv

# can do rudimentary line plots for absolute, percentage difference, or simple difference

AU2ANG = 1 / 0.5291772103


def read(filename):
    extension = filename.split('.')[-1]
    if extension == 'mat':
        # save a mat file
        print('Detected .mat file %s' % filename)
        data = spio.loadmat(filename)
        q = data['q'][0]
        res = data['I'][0]

    elif extension == 'npy':
        print('Detected .npy file %s' % filename)
        data = np.load(filename)
        q = data[:, 0]
        res = data[:, 1]

    else:
        print('Detected .txt file %s' % filename)
        data = np.genfromtxt(filename)
        q = data[:, 0]
        res = data[:, 1]

    return q, res


def perc_diff(res, ref_data, exc_frac):
    return exc_frac * (res - ref_data) / ref_data


def check_arrays_the_same(q1, q2):
    if len(q1) != len(q2):
        print(len(q1), len(q2))
        return False
    for i in range(len(q1)):
        if not np.isclose(q1[i], q2[i]):
            print(q1[i], q2[i])
            return False
    return True


def plotter(files):
    fig, ax = plt.subplots(figsize=(4, 4))
    cyc = cycler(color=[
        '#332288', '#88ccee', '#44aa99', '#117733', '#999933', '#ddcc77',
        '#cc6677', '#882255', '#aa4499'
    ])
    ax.set_prop_cycle(cyc)

    qmin = np.inf
    qmax = -np.inf
    resmin = np.inf
    resmax = -np.inf

    print('Detected %i files to plot' % len(files))

    for i, j in enumerate(files):
        print('Reading %s' % j)
        q, res = read(j)
        print(q, res)

        qmin = min(np.min(q), qmin)
        qmax = max(np.max(q), qmax)

        ax.set_xlabel('q / $\AA$')

        if pd or diff:
            if i == 0:
                print('Setting reference data to %s' % j)
                ref_q = q
                ref_data = res
                ref_j = j
            else:
                # check that q is the same for all graphs
                if (not check_arrays_the_same(q, ref_q)) and (pd or diff):
                    print(
                        'Trying to compare two arrays with different qs, failing'
                    )
                    sys.exit()
                if pd:
                    resmin = min(np.min(perc_diff(res, ref_data, 1.)), resmin)
                    resmax = max(np.max(perc_diff(res, ref_data, 1.)), resmin)
                else:
                    resmin = min(np.min(res - ref_data), resmin)
                    resmax = max(np.max(res - ref_data), resmin)
        else:
            resmin = min(np.min(res), resmin)
            resmax = max(np.max(res), resmax)

        if pd:
            ax.set_ylabel('Percentage difference')
            ax.set_title('Relative to %s' % ref_j)
            if i != 0:
                ax.plot(q * AU2ANG, perc_diff(res, ref_data, 1.), label=j)

        elif diff:
            ax.set_ylabel('Absolute difference')
            ax.set_title('Relative to %s' % ref_j)
            if i != 0:
                ax.plot(q * AU2ANG, res - ref_data, label=j)

        else:
            ax.set_ylabel('Absolute scattering')
            ax.plot(q * AU2ANG, res, label=j)

    ax.legend()
    ax.set_ylim([resmin, resmax])
    ax.set_xlim([qmin, qmax])

    plt.show()


def main():
    plotter(sys.argv[1:])


if __name__ == '__main__':
    diff = False
    pd = True

    if pd and diff:
        print(
            'Requested percentage difference and difference, just plotting percentage difference'
        )
        diff = False

    main()
