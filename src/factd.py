import numpy as np
from scipy.special import gamma


# factd(n), the double factorial n!!
# The variable n may be complex valued and can have any size.
# The function uses the complex Gamma routine from Scipy.
def factd(n):
    siz = n.size  # Note: Why was siz defined as an array, [siz]?
    # n = n[:]  # Note: Why that?
    p = np.cos(np.multiply(np.pi, n)) - 1
    f = np.multiply(np.multiply(np.power(2, np.divide(-p + n + n, 4)), np.power(np.pi, np.divide(p, 4))),
                    gamma(1 + np.divide(n, 2)))
    p = np.argwhere(np.around(n).all() == n and np.imag(n).all() == 0 and np.real(n).all() >= -1)
    if p.size != 0:  # Note: Is that correct? ~isempty(p)
        f[p] = np.around(f[p])
    p = np.argwhere(np.around(np.divide(n, 2)).all() == np.divide(n, 2).all() and np.imag(n).all() == 0 and np.real(n).all() < -1)
    if p.size != 0:  # Note: Is that correct? ~isempty(p)
        f[p] = np.inf
    f = np.reshape(f, (siz, 1))
    return f  # Note: Is f what has to be returned?
