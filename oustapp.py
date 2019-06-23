from math import *
import numpy as np
from scipy import signal
from fotf import fotfparam, FOTransFunc
from fotf import *
from control.matlab import *
from matplotlib import pyplot as plt

__all__  = ['oustapp', 'test','testbode']

def oustapp(*args):

    """
    Obtains an integer-order approximation of a fractional-order transfer function object.

    Usage:
        oustapp(G, wb, wh, N, method)

    Params:
        G - fotf object

        wb - lower bound frequency

        wh - higher bound frequency

        N - approximation order

        method -    'oust' for Oustaloup's method (default) or
                    'ref' for the refined Oustaloup method

    Defaults:
        wb = 0.001, wh = 1000, N=5, method='oust')

    Returns:
        zpk()
    """

    #Default method to use is args is less than 5
    if len(args) < 5:
        method = 'oust'
    else:
        method = args[4].lower()
        if method != 'oust' and method != 'ref':
            raise Warning('OUSTAPP.oustapp:BadMethod', "Method must be 'oust' or 'ref'. Using 'oust' as default")

    #s = tf([1],[1]) #check this

    # Order of approximation is set to default of 5 is not given
    if len(args) < 4:
        N = 5
    else:
        N = args[3]
    #setting default upper frequency bound if not given to 3
    if len(args) < 3:
        wh = 10000
    else:
        wh = args[2]

    # setting default lower frequency bound if not given to 3
    if len(args) < 2:
        wb = 0.0001
    else:
        wb = args[1]

    # raise error is no input
    if len(args) < 1:
        raise ValueError ('oustapp: NotEnoughInputArguments', 'Not enough input arguments')

    method = method.lower()

    #Get fotf paramenters
    if isinstance(args[0], FOTransFunc):
        [num, nnum, den, nden, dt] = fotfparam(args[0])

    #Go through zero array
    zeroPoly = _approxfotf(num, nnum, wb, wh, N, method)
    polePoly = _approxfotf(den, nden, wb, wh, N, method)



    # Convert to ZPK model
    # zeroZ, zeroP, zeroK = tf2zpk(zeroPoly.num[0][0], zeroPoly.den[0][0])
    # poleZ, poleP, poleK = tf2zpk(polePoly.num[0][0], polePoly.den[0][0])
    #
    # ZPKfrac = [zeroZ/poleZ, zeroP/poleP, zeroK / poleK]
    fractf = zeroPoly/ polePoly
    fractf.num[0][0] = fractf.num[0][0]/fractf.den[0][0][0]     #deviding by the first coefficient of the denominator to be similar to matlab
    fractf.den[0][0] = fractf.den[0][0] / fractf.den[0][0][0]   #deviding by the first coefficient of the denominator to be similar to matlab


    if dt >=0:
        fractf.dt = dt
    return fractf


def _approxfotf(num, nnum, wb, wh, N, method='oust'):
    """
    Approximates an FOTF  Object to a transfer function object of integer order
    :param num: coefficient of Fractional Polynomial
    :type num:  numpy.ndarray
    :param nnum: exponents of Transfer Function polynomial
    :type nnum:  numpy.ndarray
    :param wb:  Lower bound Frequency
    :type wb:   int or float
    :param wh:  Upper bound Frequency
    :type wh:   int or float
    :param N:   The desired order of the approximation must be greater than 0
    :type N: int
    :param method: Desired method 'oust' or 'ref'
    :return:
    """
    zeroPoly = 0
    for i in range(num.size):
        thisExp = nnum[i]
        intPart = int(nnum[i])
        fracPart = float(nnum[i] - intPart)
        toadd = (tf([1,0],1)**intPart)*num[i]
        toaddDen = [1]


        if fracPart != 0:
            if method == 'oust':
                toadd *=  _oustafod(fracPart, N, wb, wh)
            else:
                toadd *= tf(new_fod(fracPart, N, wb, wh))

        zeroPoly += toadd

        # approxSept2 = np.polymul(approx,np.array(seperated.num))
        # zeroPoly2+= approxSept2
    # Go through Poles array
    return zeroPoly

def _oustafod(r,N,wb,wh):
    """
    _oustfod(r,N,wb,wh): computes the Oustaloup filter approximation of a
    fractional-order operator s^r of the order N and valid in the frequency
    range (wb, wh). The function returns a ZPK object containing the continuous-time filter.
    The following equation is used to construct the filter:

            N
           -----
    Gp =    | |    (s+w'_k)/(s+wk)
            | |
            k= -N

    wk  = (b*wh/d)**((r+2k)/(2N+1))
    w'_k = (d*wb/b)**((r-2k)/(2N+1))

    :param r:   an exponent coefficient e.g. s^r, where (0 < r < 1)
    :type r:    float
    :param N:   The desired order of the approximation must be greater than 0
    :type N:    int
    :param wb:  frequency lower band
    :type wb:   float
    :param wh:  frequency upper band
    :type wh:   float
    :returns : Transfer Function object signal.lti()

    """
    if isinstance(r,float) and int(r)== 0:
        pass
    else:
        raise ValueError("oustapp._oustafod: r, must be in range (0 < r < 1)")

    if isinstance(N, int) and N > 0:
        pass
    else:
        raise ValueError("N must be an integer greater than 0")

    if isinstance(wh,(float,int)) and isinstance(wb,(float,int)) and wh > wb:
        pass
    else:
        raise ValueError("wh must be greater than wb")

    wu = wh/wb
    wb = -1*wb  #The minus sign is important. Was the cause of many debugging issues when compare with matlab
    w_kz = [(wu**((kz + N + 0.5*(1-r)) / (2 * N + 1)))* wb for kz in range(-N,N+1,1)]  #Zeros
    w_kp = [(wu**((kp + N + 0.5*(1+r)) / (2 * N + 1)))* wb for kp in range(-N,N+1,1)]  #Poles
    K = wh ** r  # gain
    ttf = zpk2tf(w_kz, w_kp, K )

    return tf(ttf[0],ttf[1])

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
        return tff

def test():
    from fotf import fotf
    g1 = newfotf(1., '14994s^{1.31}+6009.5s^{0.97}+1.69', 0)
    g2 = newfotf(1., '0.8s^{2.2}+0.5s^{0.9}+1', 0)
    g3 = newfotf('-2s^{0.63}+4', '2s^{3.501}+3.8s^{2.42}+2.6s^{1.798}+2.5s^{1.31}+1.5', 0)
    x1 = g1.oustapp( 0.0001, 10000, 5, 'oust')
    x2 = oustapp(g1, 0.0001, 10000, 5, 'oust')
    x3 = oustapp(g3, 0.0001, 10000, 5, 'oust')

    print(x1)
    print(x2)
    print(x3)
    w = np.linspace(0.01, 6.284, 1000)
    t = linspace(0, 30, 300)

    (Y1, T1) = step(x1, t)
    (Y2, T2) = step(x2, t)
    (Y3, T3) = step(x3, t)

    #plot approximation 3
    plt.figure()
    plt.plot(T1, Y1)
    plt.grid(True, axis='both', which='both')
    plt.title("1/14994s^{1.31}+6009.5s^{0.97}+1.69")
    plt.show()

    #plot approximation 2
    plt.figure()
    plt.plot(T2, Y2)
    plt.grid(True, axis='both', which='both')
    plt.title("1/0.8s^{2.2}+0.5s^{0.9}+1")
    plt.show()

    #plote approximation 3
    plt.figure()
    plt.plot(T3, Y3)
    plt.grid(True, axis='both', which='both')
    plt.title('-2s^{0.63}+4 / 2s^{3.501}+3.8s^{2.42}+2.6s^{1.798}+2.5s^{1.31}+1.5')
    plt.show()

def testbode():
    from scipy import signal
    from matplotlib import pyplot as plt
    from math import pi
    import numpy as np
    k = 1.0 / (4000.0 * pi)
    system = signal.TransferFunction([1], [k, 1])
    f = np.logspace(-1, 10, 1000) #generate 1000 values in log scale means 10**-1 to 10**2
    kk = 2 * pi * f #angular frequency
    w, mag, phase = signal.bode(system,kk) #phase is in degrees
    plt.figure()
    plt.semilogx(w, mag)  # Bode magnitude plot, w is xvalue, mag is yvalue
    plt.show()
    plt.figure()
    plt.semilogx(w, phase)  # Bode phase plot
    plt.show()
    #plt(w,phase,'r-^'), plt(w,phase,'b-o', label=phase), plt.plot(x, y)
    #scatter(x,y),subplot(2,1,1), legend(['phase','w'])
    #xlabel('Frequency'), ylabel('Phase'),title('Frequency Response'),grid()
    #clf(). close(), close('all'), colorbar(), hist(), axis('tight'), tight_layout()
    #from mayavi import mlab

    #ax1 = subploy(2,2,1)
    #subplot (2,2,2, sharex = ax1, sharey = ax1)
    #cumsum(array([1,2,3]))






