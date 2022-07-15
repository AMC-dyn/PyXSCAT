import numpy as np
import scipy.io as spio


def writer(filename, q, data):
    extension = filename.split('.')[-1]

    if extension == 'mat':
        # save a mat file
        print('Detected that .mat file requested, saving to %s' % filename)
        spio.savemat(filename, {'q': q, 'I': data})

    elif extension == 'npy':
        print('Detected that .npy file requested, saving to %s' % filename)
        np.save(filename, np.c_(q, data))

    else:
        print('Saving file as .txt into %s' % filename)
        np.savetxt(filename, np.c_(q, data), fmt='%.18e')


