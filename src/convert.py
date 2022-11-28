import numpy as np
from numba import njit

# Hide compilation warnings
from numba.core.errors import NumbaDeprecationWarning, NumbaPendingDeprecationWarning
import warnings
warnings.simplefilter('ignore', category=NumbaDeprecationWarning)
warnings.simplefilter('ignore', category=NumbaPendingDeprecationWarning)

@njit
def orbital_to_cartesian(ell: list[float] | np.ndarray[float], GM: float) -> tuple[list[float], list[float]]:
    """
    Converts orbital elements to cartesian coordinates
    :param ell: orbital elements [a, e, i, Ω, ω, M]
    :param GM: gravitational parameter
    :return: cartesian coordinates and speeds ([x, y, z], [vx, vy, vz])
    """

    a = ell[0]
    k = ell[1]
    h = ell[2]
    q = ell[3]
    p = ell[4]
    l = ell[5]

    n = np.sqrt(GM / a ** 3)
    fle = l-k*np.sin(l)+h*np.cos(l)
    corf = (l - fle + k * np.sin(fle) - h * np.cos(fle)) / (1 - k * np.cos(fle) - h * np.sin(fle))
    shift = 1e-8

    while abs(corf) > shift:
        fle += corf
        corf = (l - fle + k * np.sin(fle) - h * np.cos(fle)) / (1 - k * np.cos(fle) - h * np.sin(fle))
        shift *= 1.1

    lf = -k * np.sin(fle) + h * np.cos(fle)
    sam1 = -k * np.cos(fle) - h * np.sin(fle)
    asr = 1 / (1 + sam1)

    phi = np.sqrt(1 - k*k - h*h)
    psi = 1/(1+phi)

    x1 = a * (np.cos(fle) - k -psi * h * lf)
    y1 = a * (np.sin(fle) - h + psi * k * lf)

    vx1 = n * asr * a * (-np.sin(fle) - psi * h * sam1)
    vy1 = n * asr * a * (np.cos(fle) + psi * k * sam1)

    cis2 = 2 * np.sqrt(1 - p*p - q*q)

    tp = 1 - 2 * p*p
    tq = 1 - 2 * q*q
    dg = 2 * p * q
    
    pos = np.empty(3)
    vit = np.empty(3)

    pos[0]=x1*tp+y1*dg
    pos[1]=x1*dg+y1*tq
    pos[2]=(-x1*p+y1*q)*cis2
    vit[0]=vx1*tp+vy1*dg
    vit[1]=vx1*dg+vy1*tq
    vit[2]=(-vx1*p+vy1*q)*cis2

    return pos, vit

if __name__ == "__main__":
    print(orbital_to_cartesian([1, 0, 0, 0, 0, 0], 1))