from control.xferfcn import *
from numpy import (sqrt)
from scipy import signal


def oustapp(*args):

    """
    oustapp(G, wb, wh, N, method), Obtain a integer-order approximation of a fractional-order system.

    Args:
        G - fotf object
        wb, wh - lower and higher frequency bounds
        N - approximation order
        method - 'oust' for Oustaloup's method (default) or
                 'ref' for the refined Oustaloup method

        defaults: wb = 0.001, wh = 1000, N=5, method='oust')

        See also: fsparam

    """
    #Default method to use is args is less than 5
    if len(args) < 5:
        method = 'outs'
    else:
        method = str(args[4]).lower()
        if method != 'outs' or method!= 'ref':
            raise Warning('OUSTAPP:BadMethod', "Method must be 'oust' or 'ref'. Using 'oust' as default")

    s = tf('s')

    #get FOTF Object parameters
    if isinstance(args[0], FOTransFunc):
        [a, na, b, nb, iodelay] = fotfparam(args[0])

    # Order of approximation is set to default of 3 is not given
    if len(args) < 4:
        N = 5
    #setting default upper frequency bound if not given to 3
    if len(args) < 3:
        wh = 1000

    # setting default lower frequency bound if not given to 3
    if len(args) < 2:
        wb = 0.001

    # raise error is no input
    if len(args) < 1:
        raise ValueError('oustapp: NotEnoughInputArguments', 'Not enough input arguments')

    method
def oustafod(*args):
    """

    param r:   signifies s^r

    param N:   order

    param wb:  frequency lower band

    param wh:  frequency upper band

    return: lti object


    oustfod(r,N,wb,wh): computes the Oustaloup filter approximation of a
    fractional-order operator s^r of the order N and valid in the
    frequency range of (wb, wh). The function returns a ZPK object
    containing the continuous-time filter. The following equation
    is used to construct the filter:

             N
           -----
    Gp =    | |    (s+w'_k)/(s+wk)
            | |
            k= -N
    
    wk  = (b*wh/d)^((r+2k)/(2N+1))
    w'_k = (d*wb/b)^((r-2k)/(2N+1)).
    """

    if len(args) < 4:
        raise ValueError('oustafod: Not enough input arguments')
    elif len(args) == 4:
        r = args[0]
        N = args[1]
        wb = args[2]
        wh = args[3]

    wu = sqrt(wh / wb)

    w_kz = [wb*(wu**((kz + N + 0.5 - (0.5 * r)) / (2 * N + 1))) for kz in range(1,N+1,1)] #Zeros
    w_kp = [wb*(wu**((kp + N + 0.5 + (0.5 * r)) / (2 * N + 1))) for kp in range(1,N+1,1)]  #Poles
    K = wh ** r  # gain

    tff = signal.lti(w_kz, w_kp, K)
    dtff = signal.dlti(w_kz, w_kp, K)
    return [tff,dtff]


def new_fod(*args):
    """
    new_fod(r, N, wb, wh, b, d) creates a ZPK object with a refined
    Oustaloup filter approximation of a fractional-order operator
    s^r of order N and valid within frequency range (wb, wh).

    s^r = (d*wh/b)^r * (ds^2 + b*wh*s)/(d*(1-r)*s^2+b*wh*s+d*r)*Gp
    where
             N
           -----
    Gp =    | |    (s+w'_k)
            | |  ----------
            | |    (s+wk)
    k= -N

    wk  = (b*wh/d)^((r+2k)/(2N+1))
    w'k = (d*wb/b)^((r-2k)/(2N+1)).

    Should parameters b and d be omitted, the default values will be
    used: b=10, d=9.

    :param r:   signifies s^r
    :param N:   order
    :param wb:  frequency lower band
    :param wh:  frequency upper band
    :param b: 
    :param d: 
    :return: Fractional Order Object
    """

    if len(args)< 4:
        raise ValueError('new_fod: Not enough input arguments')
    elif len(args) == 4:
        r = args[0]
        N = args[1]
        wb = args[2]
        wh = args[3]
        b = 10
        d = 9
    elif len(args) == 6:
        r = args[0]
        N = args[1]
        wb = args[2]
        wh = args[3]
        b = args[4]
        d = args[5]

    if r == 0:
        return signal.lti(0)
    else:
        mu = wh / wb

        w_kz = [wb * (mu ** ((k + N + 0.5 - (0.5 * r)) / (2 * N + 1))) for k in range(1, N + 1, 1)]
        w_kp = [wb * (mu ** ((k + N + 0.5 + (0.5 * r)) / (2 * N + 1))) for k in range(1, N + 1, 1)]

        K = pow((d * wh / b), r)
        tff = signal.lti((w_kz,w_kp, K)*tf([d, b * wh, 0], [d * (1 - r), b * wh, d * r]))
        dtff = signal.dlti((w_kz,w_kp, K)*tf([d, b * wh, 0], [d * (1 - r), b * wh, d * r]))
        return [tff, dtff]

def fotfparam(FOTransFunc):
    return