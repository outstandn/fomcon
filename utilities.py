import pandas as pd
from enum import Enum
from datetime import datetime
import numpy as np
from fotf import *
from scipy.optimize import least_squares,leastsq, curve_fit, shgo, dual_annealing, basinhopping, differential_evolution, Bounds
from control import matlab
from matplotlib import pyplot as plt
__all__ = ['optType', 'optAlgo', 'optFix', 'opt', 'fid', 'test']

def test():
    typset = optType.grunwaldLetnikov
    algset = optAlgo.TrustRegionReflective
    fixset = optFix.Exp
    guessset = newfotf('2s^{0.63}+4', '2s^{3.501}+3.8s^{2.42}+2.6s^{1.798}+2.5s^{1.31}+1.5', 0)
    polyfixset = [0, 0]
    optiset = opt(guessset,typset,algset,fixset,polyfixset)
    result = fid('PROC1.xlsx','PROC2.xlsx',optiset,[[0,20],[0,10]])
    print(result)




class optType(Enum):
    oustaloop='oust'
    grunwaldLetnikov = 'gl'

class optAlgo(Enum):
    LevenbergMarquardt = 'lm'
    TrustRegionReflective = 'trf'
    DogBox = 'dogbox'
    RobustLoss = 'cauchy'


class optFix(Enum):
    Free = 'n'
    Coeff = 'c'
    Exp = 'e'

class opt():

    def __init__(self, initialGuess, optiType, optiAlg, optiFix, polyFix, funcFix=None):
        """
        Usage:              opt(initialGuess, optiType, optiAlg, optiFix, polyFix)

        :param initialGuess: initial guessed fractional order transfer function
        :type initialGuess: str or list or FOTransFunc
        :param optiType:
        :param optiAlg:
        :param optiFix:
        :param polyFix:     a vector with two values: [BFIX; AFIX], where BFIX and AFIX can be 1, in which case the corresponding
                            polynomial is fixed during identification, or 0. Note that in case BFIX = AFIX = 1 the initial model
                            will be immediately returned with no identification conducted. Default: [0; 0]
        :type polyFix:      list or numpy.ndarray or tuple
        :param funcFix:

        """
        if isinstance(initialGuess, str):
            self.G = newfotf(initialGuess)
        elif isinstance(initialGuess, (list, int, float)):
            self.G = fotf(initialGuess)
        elif isinstance(initialGuess,FOTransFunc):
            self.G = initialGuess
        else:
            raise ValueError("utilities.opt: initialGuess should be of type 'str' or 'list' or 'FOTransFunc'")

        if isinstance(optiType, optType):
            self.type = optiType
        else:
            raise ValueError("utilities.opt: 2nd parameter should be of type optType")

        if isinstance(optiAlg, optAlgo):
            self.alg = optiAlg
        else:
            raise ValueError("utilities.opt: 3rd parameter should be of type optAlgo")

        if isinstance(optiFix, optFix):
            self._optiFix = optiFix
        else:
            raise ValueError("opt.fix: 4th parameter should be of type optAlgo")

        if isinstance(polyFix, (list, np.ndarray)):
            self._polyFix = np.array(polyFix)
        else:
            raise ValueError("utilities.opt: 5th Parameter should be of type list or numpy.ndarray")

        self._funcFix = funcFix
    @property
    def funcFix(self):
        return self._funcFix

    @funcFix.setter
    def funcFix(self,value):
        if isinstance(value,np.ndarray):
            self._funcFix = value
        else:
            raise ValueError("opt.funcFix: 'numpy.ndarray' type is allowed")

    @property
    def polyFix(self):
        return self._polyFix

    @polyFix.setter
    def polyFix(self, value):
        if isinstance(value, (list, np.ndarray, tuple)):
            self._polyFix = value
        else:
            raise ValueError("opt.polyFix: polyFix should be of type list or tuple or numpy.ndarray")

    @property
    def optiFix(self):
        return self._optiFix

    @optiFix.setter
    def optiFix(self, value):
        if isinstance(value, optFix):
            self._optiFix = value
        else:
            raise ValueError("opt.optiFix: 4th parameter should be of type optAlgo")









def _fracidfun(x0, y, u, t, opti):
    """
    Compute fractional time-domain identification cost function.

    Usage: z = _fracidfun(x0,y,u,t,opt)

          opt.type - 'gl', 'oust'

          opt.wb, opt.wh, opt.N - parameters to be used
                                  with 'oust'
          opt.fix - identification type, 'n' for free identification,
                    'c' to fix coefficients, 'e' to fix exponents
          opt.fixpoly - a matrix [zeroPolyLength; polePolyLength],
                        if length is zero, the polynomial is considered
                        fixed and derived from initial model. Both lengths
                        cannot be zero at once.
          opt.G - initial model
          opt.alg - identification algorithm (1: TRR, 2: LM)
    :param x0:  random guesses of x0 used to test in objective function
    :type x0:   numpy.ndarray
    :param y:   Output From Data
    :type y:    numpy.ndarray
    :param u:   Input from data
    :type u:    numpy.ndarray
    :param t:   time in seconds
    :type t:    numpy.ndarray
    :param opti:optimise options
    :type opti: opt
    :return err:y-y_id
    :type err: numpy.ndarray
    """

    # Cannot identify if both polynomials fixed
    if isinstance(opti.funcFix, (list,np.ndarray, tuple)):
        # Check fixpoly
        if opti.funcFix[0] == 0 and opti.funcFix[1] == 0:
            raise ValueError('FRACIDFUN:BothPolynomialsFixed: Cannot identify model because both polynomials are set to be fixed')
        else:
            numSize = int(opti.funcFix[0])
            denSize = int(opti.funcFix[1])
    else:
        raise ValueError('FRACIDFUN: wrong type for opti.polyFix: list, tuple, ndarray allowed')


    # Get initial model parameters
    [inum, innum, iden, inden,idt] = fotfparam(opti.G)
    opt = opti.optiFix

    # Pole polynomial is fixed
    if numSize > 0 and denSize == 0:
        den = iden
        nden = inden
        [num, nnum] = _fracidfun_getpolyparam(opt, numSize, x0, inum, innum)

    # Zero polynomial is fixed
    elif numSize == 0 and denSize > 0:

        num = inum
        nnum = innum
        [den, nden] = _fracidfun_getpolyparam(opt, denSize, x0, iden, inden)

    # Free identification
    elif numSize > 0 and denSize > 0:
        [num, nnum, den, nden] = _fracidfun_getfotfparam(opt, numSize, denSize, x0, iden, inden, inum, innum)

    # Get identified fotf object
    G = FOTransFunc(num,nnum,den,nden,idt)

    # Build model based on type
    if opti.type ==optType.grunwaldLetnikov:
        y_id = lsim(G,u,t)
    elif opti.type == optType.oustaloop:
        G = G.oustapp()
        (y_id, t) = matlab.step(G,t)
    else:
        raise  ValueError("utilities._fracidfun: Unknown simulation type 'optType' specified!")
    return y - y_id

# Returns polynomial parameters based on desired identification type
def _fracidfun_getpolyparam(fix, vec_len, vec, ip, inp):
    if fix == optFix.Free:
        # Free identification
        p = vec[0: vec_len/2]
        np = vec[vec_len/2:]
    elif fix == optFix.Exp:
        # Fix exponents
        p = vec
        np = inp
    elif fix == optFix.Coeff:
        # Fix coefficients
        p = ip
        np = vec
    return [p,np]

# Returns both polynomial parameters based on desired identification type
def _fracidfun_getfotfparam(fix, numSize, denSize, vec, ia, ina, ib, inb):

    if fix == optFix.Free:
        b = vec[0:int(numSize/2)]
        nb = vec[int(numSize/2):numSize]
        a = vec[numSize:numSize+int(denSize/2)]
        na = vec[numSize+int(denSize/2):]
    elif fix == optFix.Exp:
        b = vec[0:numSize]
        nb = inb
        a = vec[numSize:]
        na = ina
    elif fix == optFix.Coeff:
        b = ib
        nb = vec[0:numSize]
        a = ia
        na = vec[numSize:]

    return [b, nb, a, na ]

def fid(idd, vidd, opti, limits=None ):
    """

    :param idd:     Name of Excel file (with .extension). File should have heading(y,u,t,dt).Should be in working directory of Source code
    :type idd:      str
    :param vidd:     Name of Verification Excel file (with .extension). File should have heading(y,u,t,dt).Should be in working directory of Source code
    :type vidd:      str

    :param limits:  a cell array with two vectors, which may also be empty, containing polynomial coefficient and
                    exponent limits in the form [[CMIN;CMAX],[EMIN;EMAX]]. Default: [] (i.e. no limits are imposed.)
    :type limits:   list,np.ndarray
    :param opti:    see opt()
    :type opti:     opt
    :return:    FOTransFunc()
    """

    if isinstance(idd,str):
        data = pd.read_excel(idd)
        y = np.array(data.y.values)
        u = np.array(data.u.values)
        t = np.array(data.t.values)
        # dt = t[1]-t[0]

        plt.figure()

        del (data)  # to save memory
        if y.size != u.size or u.size != t.size or t.size != y.size:
            raise IOError("utilies.fid: size of data idd are not the same. Kindly Fix your data")
    else:
        raise ValueError("utilities:fid: idd should be a string to an excel file (Extension should be include)")

    if isinstance(vidd,str):
        data = pd.read_excel(vidd)
        vy= np.array(data.y.values)
        vu = np.array(data.u.values)
        vt = np.array(data.t.values)
        # vdt = vt[1]-vt[0]

        del (data)  # to save memory

        if vy.size != vu.size or vu.size != vt.size or vt.size != vy.size:
            raise IOError("utilies.fid: size of data vidd are not the same. Kindly Fix your data")

    else:
        raise ValueError("utilities:fid: vidd should be a string to an excel file (Extension should be include)")


    #check limits
    if limits == None:
        clim = np.array([0,np.inf])
        elim = np.array([0,10])
    else:
        clim = np.array(limits[0])
        elim = np.array(limits[1])

        if clim[0] > clim[1]:
            raise Warning('Swapping limits')
            temp = clim[0]
            clim[1] = clim[0]
            clim[0] = temp

        if elim[0] > elim[1]:
            raise Warning('Swapping limits')
            temp = elim[0]
            elim[1] = elim[0]
            elim[0] = temp
        if elim[0] < 0:
            raise Warning('exponent lower bound must be >= 0')
            elim[0]=0
        if clim[0] < 0:
            raise Warning('Coefficent lower bound must be >= 0')
            clim[0] = 0

    xnum,xnnum,xden,xnden,xdt = fotfparam(opti.G)

    #Check polynomial fix options
    if opti.polyFix[0] == 1 and opti.polyFix[1] == 1:
        #No optimization was done
        return opt.G
    elif opti.polyFix[0] == 0 and opti.polyFix[1] == 1: #Poles polynomial is fixed i.e not optimized
        #check which to optimize coefficients, exponents or both
        if opti.optiFix == optFix.Free:
            #initial Guess
            x0 = np.concatenate([xnum,xnnum])
            #limits
            lb = np.concatenate([clim[0]*np.ones(x0.size/2), elim[0]*np.ones(x0.size/2)])
            ub = np.concatenate([clim[1] * np.ones(x0.size / 2), elim[1] * np.ones(x0.size / 2)])
        elif opti.optiFix == optFix.Exp:
            #Fix Exponent, optimize Coefficient
            x0 = xnum
            #limits
            lb = clim[0] * np.ones_like(x0)
            ub = clim[1] * np.ones_like(x0)
        elif opti.optiFix == optFix.Coeff:
            #Fix Coefficient, optimize Exponent
            x0 = xnnum
            # limits
            lb = elim[0] * np.ones_like(x0)
            ub = elim[1] * np.ones_like(x0)

        opti.funcFix = np.array([x0.size, 0])

    elif opti.polyFix[0] == 1 and opti.polyFix[1] == 0:#Zeroes polynomial is fixed i.e not optimized
        # check which to optimize coefficients, exponents or both
        if opti.optiFix == optFix.Free:
            # initial Guess
            x0 = np.concatenate([xden, xnden])
            # limits
            lb = np.concatenate([clim[0] * np.ones(x0.size / 2), elim[0] * np.ones(x0.size / 2)])
            ub = np.concatenate([clim[1] * np.ones(x0.size / 2), elim[1] * np.ones(x0.size / 2)])
        elif opti.optiFix == optFix.Exp:
            # Fix Exponent, optimize Coefficient
            x0 = xden
            # limits
            lb = clim[0] * np.ones_like(x0)
            ub = clim[1] * np.ones_like(x0)
        elif opti.optiFix == optFix.Coeff:
            # Fix Coefficient, optimize Exponent
            x0 = xnden
            # limits
            lb = elim[0] * np.ones_like(x0)
            ub = elim[1] * np.ones_like(x0)

        opti.funcFix = np.array([0, x0.size])

    else: #No fix, optimise Zero Polynomials and Poles Polynomials
        if opti.optiFix == optFix.Free:
            # initial Guess
            x0 = np.concatenate([xnum, xnnum, xden, xnden])
            # limits
            lb = np.concatenate([clim[0] * np.ones_like(xnum), elim[0] * np.ones_like(xnnum), clim[0] * np.ones_like(xden), elim[0] * np.ones_like(xnden)])
            ub = np.concatenate([clim[1] * np.ones_like(xnum), elim[1] * np.ones_like(xnnum), clim[1] * np.ones_like(xden), elim[1] * np.ones_like(xnden)])
        elif opti.optiFix == optFix.Exp:
            # Fix Exponent, optimize Coefficient
            x0 = xden
            # limits
            lb = clim[0] * np.ones_like(x0)
            ub = clim[1] * np.ones_like(x0)
        elif opti.optiFix == optFix.Coeff:
            # Fix Coefficient, optimize Exponent
            x0 = xnden
            # limits
            lb = elim[0] * np.ones_like(x0)
            ub = elim[1] * np.ones_like(x0)

        opti.funcFix = np.array([xnum.size*2, xden.size*2])

    print("Please wait,  System Identification in progress...")

    start = datetime.now()
    #Run the identification Algorithm
    if opti.alg == optAlgo.LevenbergMarquardt:
        #bounds will not be used
        res = least_squares(_fracidfun, x0, args=(y,u,t,opti), verbose=1, method='lm')
    elif opti.alg == optAlgo.TrustRegionReflective:
        res = least_squares(_fracidfun, x0, args=(y,u,t,opti), verbose=1, method='trf', bounds=(lb,ub))
    elif opti.alg == optAlgo.DogBox:
        res = least_squares(_fracidfun, x0, args=(y, u, t, opti), verbose=1, bounds=(lb,ub), method='dogbox')
    elif opti.alg == optAlgo.RobustLoss:
        res = least_squares(_fracidfun, x0, args=(y, u, t, opti), verbose=1, bounds=(lb,ub), loss = 'cauchy', f_scale = 0.1)

    print('Time: ',datetime.now() - start)

    #Display Resuls
    print(res.x,'\n\n', res.message, '\n\n', res.success, '\n\n', 'Number of Function Calls : {}'.format(res.nfev), '\n\n', 'Function Residuals: {}'.format(res.fun),'\n\n', 'Value at Solution: {}'.format(res.cost))

    if opti.polyFix[0] == 0 and opti.polyFix[1] == 1:  # Poles polynomial is fixed i.e not optimized
        # so update zeros with optimized result
        if opti.optiFix == optFix.Free:
            xnum = res.x[0:x.size/2]
            xnnum = res.x[x.size/2:]
        elif opti.optiFix == optFix.Exp:
            xnum = res.x
        elif opti.optiFix == optFix.Coeff:
            # Fix Coefficient, update Exponent
            xnnum = res.x
        opti.funcFix = np.array([x0.size, 0])

    elif opti.polyFix[0] == 1 and opti.polyFix[1] == 0:  # Zeroes polynomial is fixed i.e not optimized
        # update optimize coefficients in poles
        if opti.optiFix == optFix.Free:
            xden = res.x[0:x.size / 2]
            xnden = res.x[x.size / 2:]
        elif opti.optiFix == optFix.Exp:
            # Fix Exponent, so update optimize Coefficient
            xden = res.x
        elif opti.optiFix == optFix.Coeff:
            # Fix Coefficient, so update optimize Exponent
            xnden = res.x

    else:  # No fix, optimise Zero Polynomials and Poles Polynomials
        if opti.optiFix == optFix.Free:
            xnum = res.x[0:xnum.size]
            xnnum = res.x[xnum.size:xnum.size*2]
            xden = res.x[xnum.size*2:xden.size+(xnum.size*2)]
            xnden = res.x[xden.size+(xnum.size*2):]
        elif opti.optiFix == optFix.Exp:
            # Fix Exponent, update  optimized Coefficient
            xnum = res.x[0:xnum.size]
            xden = res.x[xnum.size :]

        elif opti.optiFix == optFix.Coeff:
            # Fix Coefficient, update optimized Exponent
            xnnum = res.x[0:xnum.size]
            xnden = res.x[xnum.size:]

    return FOTransFunc(xnum,xnnum,xden,xnden, xdt)





