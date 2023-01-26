import numpy as np
from numba import njit # optional (on-the-fly compilation)

# Hide compilation warnings
from numba.core.errors import NumbaDeprecationWarning, NumbaPendingDeprecationWarning
import warnings
warnings.simplefilter('ignore', category=NumbaDeprecationWarning)
warnings.simplefilter('ignore', category=NumbaPendingDeprecationWarning)

@njit # optional (on-the-fly compilation)
def orbital_to_cartesian(ell: list[float] | np.ndarray[float], GM: float) -> tuple[list[float], list[float]]:
    # Remove these type hints if you are running python < 3.10
    """
    Converts orbital elements to cartesian coordinates
    :param ell: orbital elements [a, e, i, Ω, ω, λ]
    :param GM: gravitational parameter
    :return: cartesian coordinates and speeds ([x, y, z], [vx, vy, vz])
    """

    # Parsing orbital elements to original function parameters
    a = ell[0]
    e = ell[1]
    i = ell[2]
    omega = ell[3]
    varpi = ell[4]
    l = ell[5]

    k = e*np.cos(varpi)
    h = e*np.sin(varpi)
    q = np.sin(i/2)*np.cos(omega)
    p = np.sin(i/2)*np.sin(omega)

    # The original function converted to python (not my code, don't ask me how it works)
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


def cartesian_to_orbital(pos: list[float] | np.ndarray[float], vit: list[float] | np.ndarray[float], GM: float) -> list[float] | np.ndarray[float]:
    # Remove these type hints if you are running python < 3.10
    """
    Converts cartesian coordinates to orbital elements
    :param pos: cartesian coordinates [x, y, z]
    :param vit: cartesian speeds [vx, vy, vz]
    :param GM: gravitational parameter
    :return: orbital elements [a, e, i, Ω, ω, λ]
    """

    # The original function converted to python (not my code, don't ask me how it works)
    pos = np.array(pos)
    vit = np.array(vit)

    rayon = np.sqrt(np.sum(pos**2))
    v2 = np.sum(vit**2)
    a = GM * rayon / (2 * GM - rayon * v2)
    gx = pos[1] * vit[2] - pos[2] * vit[1]
    gy = pos[2] * vit[0] - pos[0] * vit[2]
    gz = pos[0] * vit[1] - pos[1] * vit[0]
    gg = np.sqrt(gx*gx + gy*gy + gz*gz)
    # print(gz, gg)
    cis2 = np.sqrt(0.5 * (1 + gz / gg))
    # print(gg, cis2)
    q = -gy / (2 * gg * cis2)
    p = gx / (2 * gg * cis2)

    tp = 1 - 2 * p*p
    tq = 1 - 2 * q*q
    dg = 2 * p * q

    x1 = tp * pos[0] + dg * pos[1] - 2 * p * cis2 * pos[2]
    y1 = dg * pos[0] + tq * pos[1] + 2 * q * cis2 * pos[2]
    vx1 = tp * vit[0] + dg * vit[1] - 2 * p * cis2 * vit[2]
    vy1 = dg * vit[0] + tq * vit[1] + 2 * q * cis2 * vit[2]

    k = gg * vy1 / GM - x1 / rayon
    h = -gg * vx1 / GM - y1 / rayon

    psi = 1/(1+np.sqrt(1-k*k-h*h))

    ach = 1 - psi * h*h
    ack = 1 - psi * k*k

    adg = psi * h * k
    det = ach * ack - adg*adg
    sm1 = x1/a + k
    sm2 = y1/a + h

    cf = (sm1*ack - sm2*adg) / det
    sf = (ach*sm2 - adg*sm1) / det

    fle = np.arctan2(sf, cf)

    l = fle - k * sf + h * cf

    # Parsing the original function outputs to the orbital elements
    e = np.sqrt(k*k - h*h)
    i = 2 * np.arcsin(np.sqrt(p**2 + q**2))
    omega = np.arctan2(p, q)
    varpi = np.arctan2(k, h)

    return np.array([a, e, i, omega, varpi, l])

if __name__ == "__main__":
    vit, pos = orbital_to_cartesian([1, 0, 0, 0, 1, 0], 1)
    GM = 1
    ell = cartesian_to_orbital(pos,vit, GM)
    for i in ell:
        print(i, end="\t")
    print()