from control.xferfcn import *
from numpy import sqrt
import numpy as np
from scipy import signal
from fotf import FOTransFunc


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
        method = 'oust'
    else:
        method = args[4].lower()
        if method != 'oust' and method != 'ref':
            raise Warning('OUSTAPP.oustapp:BadMethod', "Method must be 'oust' or 'ref'. Using 'oust' as default")

    #s = tf([1],[1]) #check this

    #get FOTF Object parameters
    if isinstance(args[0], FOTransFunc):
        [num, nnum, den, nden, iodelay] = fotfparam(args[0])

    # Order of approximation is set to default of 3 is not given
    if len(args) < 4:
        N = 5
    else:
        N = args[3]
    #setting default upper frequency bound if not given to 3
    if len(args) < 3:
        wh = 1000
    else:
        wh = args[2]

    # setting default lower frequency bound if not given to 3
    if len(args) < 2:
        wb = 0.001
    else:
        wb = args[1]

    # raise error is no input
    if len(args) < 1:
        raise ValueError('oustapp: NotEnoughInputArguments', 'Not enough input arguments')

    method = method.lower()

    #Get fotf paramenters



    #Go through zero array
    zeroPoly = approxfotf(num, nnum, wb, wh, N, method)
    polePoly = approxfotf(den, nden, wb, wh, N, method)



    # Convert to ZPK model
    fractf = zeroPoly / polePoly
    #fractzpk = signal.tf2zpk(zeroPoly,polePoly)
    return fractf


def approxfotf(num, nnum, wb, wh, N, method='oust'):
    zeroPoly = None
    for i in range(len(num[0][0])):
        thisExp = nnum[0][0][i]
        intPart = int(nnum[0][0][i])
        fracPart = float(nnum[0][0][i] - intPart)
        toadd = [num[0][0][i]]
        toaddDen = [1]

        for j in range(intPart):
            toadd.append(0)
            #toaddDen.append(0)
        seprated = tf(toadd, toaddDen)
        if fracPart != 0:
            if method == 'oust':
                approx = oustafod(fracPart, N, wb, wh)
            else:
                approxnum, approxden = tf(new_fod(fracPart, N, wb, wh))

        approxSeptf = np.array(seprated) * approx

        if i == 0:
            zeroPoly = approxSeptf
        else:
            zeroPoly = zeroPoly + approxSeptf

        # Go through Poles array
    return zeroPoly


def oustafod(*args):
    """

    param r:   signifies s^r

    param N:   order

    param wb:  frequency lower band

    param wh:  frequency upper band

    return: transfer Function OBject


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

    w_kz = [wb*(wu**((kz + N + 0.5 - (0.5 * r)) / (2 * N + 1))) for kz in range(-N,N+1,1)] #Zeros
    w_kp = [wb*(wu**((kp + N + 0.5 + (0.5 * r)) / (2 * N + 1))) for kp in range(-N,N+1,1)]  #Poles
    K = wh ** r  # gain

    return signal.zpk2tf(w_kz, w_kp, K )

# def new_fod(*args):
#     """
#     new_fod(r, N, wb, wh, b, d) creates a ZPK object with a refined
#     Oustaloup filter approximation of a fractional-order operator
#     s^r of order N and valid within frequency range (wb, wh).
#
#     s^r = (d*wh/b)^r * (ds^2 + b*wh*s)/(d*(1-r)*s^2+b*wh*s+d*r)*Gp
#     where
#              N
#            -----
#     Gp =    | |    (s+w'_k)
#             | |  ----------
#             | |    (s+wk)
#     k= -N
#
#     wk  = (b*wh/d)^((r+2k)/(2N+1))
#     w'k = (d*wb/b)^((r-2k)/(2N+1)).
#
#     Should parameters b and d be omitted, the default values will be
#     used: b=10, d=9.
#
#     :param r:   signifies s^r
#     :param N:   order
#     :param wb:  frequency lower band
#     :param wh:  frequency upper band
#     :param b:
#     :param d:
#     :return: Fractional Order Object
#     """
#
#     if len(args)< 4:
#         raise ValueError('new_fod: Not enough input arguments')
#     elif len(args) == 4:
#         r = args[0]
#         N = args[1]
#         wb = args[2]
#         wh = args[3]
#         b = 10
#         d = 9
#     elif len(args) == 6:
#         r = args[0]
#         N = args[1]
#         wb = args[2]
#         wh = args[3]
#         b = args[4]
#         d = args[5]
#
#     if r == 0:
#         return signal.lti(0)
#     else:
#         mu = wh / wb
#
#         w_kz = [wb * (mu ** ((k + N + 0.5 - (0.5 * r)) / (2 * N + 1))) for k in range(1, N + 1, 1)]
#         w_kp = [wb * (mu ** ((k + N + 0.5 + (0.5 * r)) / (2 * N + 1))) for k in range(1, N + 1, 1)]
#
#         K = pow((d * wh / b), r)
#         tff = signal.lti((w_kz,w_kp, K)*tf([d, b * wh, 0], [d * (1 - r), b * wh, d * r]))
#         dtff = signal.dlti((w_kz,w_kp, K)*tf([d, b * wh, 0], [d * (1 - r), b * wh, d * r]))
#         return [tff, dtff]

def fotfparam(fotObject):
    if isinstance(fotObject, FOTransFunc):
        return [fotObject.num, fotObject.nnum, fotObject.den, fotObject.nden, fotObject.dt]
    else:
        raise ValueError("oustapp.fotfparam, Object is not of type FOTransFunc")

def fix_s(b):
    """
    Extract integer part from a given number
    :param b: List of real numbers
    :return a: List of integer numbers
    """

    a=[]  #output array
    if isinstance(b,list):
        for number in b:
            a = a.append(int(number)) #extrct only the integer part
        return a
    else:
        raise ValueError("OUSTAPP.fix_s: input should be a list with real numbers not a".format(type(b)))

def test():
    from fotf import fotf
    from oustapp import oustapp, approxfotf
    nnum = [1.25, 0.]
    num = [1., 2.]
    den = [9., 8., 7.]
    nden = [2.5, 0.25, 0.]
    g = fotf(num, nnum, den, nden)
    x = oustapp(g, 0.01, 100, 4,'oust')
    print(x)

