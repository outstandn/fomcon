import pandas as pd
from enum import Enum
from datetime import datetime
import numpy as np
from fotf import *
from scipy.optimize import minimize, least_squares,leastsq, curve_fit, Bounds#, shgo, dual_annealing, basinhopping, differential_evolution
from control.matlab import lsim as controlsim
from matplotlib import pyplot as plt
__all__ = ['optMethod', 'optAlgo', 'optFix', 'opt', 'fid', 'optMethod', 'optAlgo', 'optFix', 'test']

def test():
    result = []
    counter = 1
    guessset = g3 = newfotf('-2s^{0.63}+4', '2s^{3.501}+3.8s^{2.42}+2.6s^{1.798}+2.5s^{1.31}+1.5', 0)

    for j in [optAlgo.LevenbergMarquardt]:
        # for k in [optFix.Exp]:
        for k in [optFix.Free, optFix.Coeff,optFix.Exp]:
            # for l in range(2):
            #     for m in range(2):
                    polyfixset = [0, 0]
                    optiset = opt(guessset, optMethod.grunwaldLetnikov, j, k, polyfixset)
                    print('{0}: Computing settings: {1}, optMethod.grunwaldLetnikov, {2}, {3}'.format(counter,  j, k, polyfixset))
                    res = fid('PROC1.xlsx', 'PROC2.xlsx', optiset, [[0, 20], [0, 10]],plot=[False, False], plotid=[False, False], cleanDelay = [True,2.5])
                    result.append(res)
                    print(res.G, "\n\n")
                    counter+=1

    guessset = g3 = newfotf('2s^{0.63}+4', '2s^{3.501}+3.8s^{2.42}+2.6s^{1.798}+2.5s^{1.31}+1.5', 0)
    for j in [optAlgo.TrustRegionReflective]:
        for k in [optFix.Free, optFix.Coeff,optFix.Exp]:
            polyfixset = [0, 0]
            optiset = opt(guessset, optMethod.grunwaldLetnikov, j, k, polyfixset)
            print('{0}: Computing settings: {1}, optMethod.grunwaldLetnikov, {2}, {3}'.format(counter,  j, k, polyfixset))
            res = fid('PROC1.xlsx', 'PROC2.xlsx', optiset, [[0, 20], [0, 10]],plot=[False, False], plotid=[False, False], cleanDelay = [True,2.5])
            result.append(res)
            print(res.G, "\n\n")
            counter+=1

    return result

    # typset = optMethod.grunwaldLetnikov
    # algset = optAlgo.RobustLoss
    # fixset = optFix.Free
    # guessset = newfotf('2s^{0.63}+4', '2s^{3.501}+3.8s^{2.42}+2.6s^{1.798}+2.5s^{1.31}+1.5', 0)
    # u = np.ones(1000)
    # t = np.linspace(0,30,1000)
    # y = lsim(guessset,u,t)
    # polyfixset = [0, 0]
    # optiset = opt(guessset,typset,algset,fixset,polyfixset)
    # result = fid('PROC1.xlsx','PROC2.xlsx',optiset,[[0,20],[0,10]], plot=[True,True])
    # return result




class optMethod(Enum):
    grunwaldLetnikov = 'gl'
    oustaloop='oust'


class optAlgo(Enum):
    LevenbergMarquardt = 'lm'
    TrustRegionReflective = 'trf'
    Softl1 = 'softl1'
    RobustLoss = 'cauchy'

class optFix(Enum):
    Free = 'n'
    Coeff = 'c'
    Exp = 'e'

class opt():

    def __init__(self, initialGuess, optiType, optiAlg, optiFix, polyFix, optidelay = False, funcFix=None):
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
        :param optidelay:   Should delay be optimized? if true
        :type optidelay:    bool or int

        """
        if isinstance(initialGuess, str):
            self.G = newfotf(initialGuess)
        elif isinstance(initialGuess, (list, int, float)):
            self.G = fotf(initialGuess)
        elif isinstance(initialGuess,FOTransFunc):
            self.G = initialGuess
        else:
            raise ValueError("utilities.opt: initialGuess should be of type 'str' or 'list' or 'FOTransFunc'")

        if isinstance(optiType, optMethod):
            self.type = optiType
        else:
            raise ValueError("utilities.opt: 2nd parameter should be of type optMethod")

        if isinstance(optiAlg, optAlgo):
            self.alg = optiAlg
        else:
            raise ValueError("utilities.opt: 3rd parameter should be of type optAlgo")

        if isinstance(optiFix, optFix):
            self._optiFix = optiFix
        else:
            raise ValueError("opt.fix: 4th parameter should be of type optAlgo")

        if optidelay == True or optidelay == 1:
            self.findDelay = True
        else:
            self.findDelay = False

        if isinstance(polyFix, (list, np.ndarray)):
            if polyFix[0] == True or 1:
                polyFix[0] = 1
            else:
                polyFix[0] = 0

            if polyFix[1] == True or 1:
                polyFix[0] = 1
            else:
                polyFix[0] = 0
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
            opti.findDelay = True
            raise ValueError('FRACIDFUN:BothPolynomialsFixed: Cannot identify model because both polynomials are set to be fixed')
        else:
            numSize = int(opti.funcFix[0])
            denSize = int(opti.funcFix[1])
    else:
        raise ValueError('FRACIDFUN: wrong type for opti.polyFix: list, tuple, ndarray allowed')


    # Get initial model parameters
    [inum, innum, iden, inden,delay] = fotfparam(opti.G)
    opt = opti.optiFix
    if opti.findDelay:
        delay = x0[-1]
        x0 = x0[:-1]

    # Pole polynomial is fixed
    if numSize > 0 and denSize == 0:
        #Update numerators (Zeroes)
        [inum, innum] = _fracidfun_getpolyparam(opt, numSize, x0, inum, innum)

    # Zero polynomial is fixed
    elif numSize == 0 and denSize > 0:
        #Update Denominator (Poles)
        [iden, inden] = _fracidfun_getpolyparam(opt, denSize, x0, iden, inden)

    # Free identification
    elif numSize > 0 and denSize > 0:
        # Update both (Zeros & Poles)
        [inum, innum, iden, inden] = _fracidfun_getfotfparam(opt, numSize, denSize, x0, iden, inden, inum, innum)

    # Get identified fotf object
    G = FOTransFunc(inum, innum, iden, inden,delay)

    # Build model based on type
    if opti.type == optMethod.grunwaldLetnikov:
        y_id = lsim(G,u,t)
    elif opti.type == optMethod.oustaloop:
        newG = G.oustapp()
        (y_id, t, x00) = controlsim(newG,u, t)
    else:
        raise  ValueError("utilities._fracidfun: Unknown simulation type 'optMethod' specified!")
    err = y - y_id
    return err

# Returns polynomial parameters based on desired identification type
def _fracidfun_getpolyparam(fix, vec_len, vec, ip, inp):
    if fix == optFix.Free:
        # Free identification
        p = vec[0: int(vec_len/2)]
        np = vec[int(vec_len/2):]
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

def fid(idd, vidd, opti, limits=None, plot = [False,False] , plotid = [True, True], cleanDelay = [True,2.5]):
    """

    :param idd:     Name of Excel file (with .extension). File should have heading(y,u,t).
                    Should be in working directory of Source code.This is the "Identification Data"
    :type idd:      str
    :param vidd:    Name of Excel file (with .extension). File should have heading(vy,vu,vt).
                    Should be in working directory of Source code. This is the "Verification Data"
    :type vidd:      str
    :param opti:    see opt()
    :type opti:     opt
    :param limits:  a cell array with two vectors, which may also be empty, containing polynomial coefficient and
                    exponent limits in the form [[CMIN;CMAX],[EMIN;EMAX]]. Default: None (i.e. no limits are imposed.)
    :type limits:   list,np.ndarray
    :param plot:    array of type bool with length=2. eg[True,True] where position 0 is to plot identification data
                    and position 1 is to plot validation data. Default [False,False]
    :type plot:     list, nummpy.ndarray
    :param plotid:  array of type bool with length=2. eg[True,True] where position 0 is to plot identified system output
                    and position 0 is to plot identified system output. Default [False, False]
    :type plotid:   list

    :param cleanDelay:  Default = [True,2.5]. Used to clean identification data with delays before identification. 2.5 is delay in seconds
    :type cleanDelay:   list

    :return :       fidOutput
    """

    EXP_LB = 0.001
    if isinstance(idd,str) and (idd[-4:] =='xlsx' or idd[-3:] =='xls'):
        data = pd.read_excel(idd)
        y = np.array(data.y.values)
        u = np.array(data.u.values)
        t = np.array(data.t.values)
        # dt = t[1]-t[0]
        del(data)  # to save memory

        if y.size != u.size or u.size != t.size or t.size != y.size:
            raise IOError("utilies.fid: size of data idd are not the same. Kindly Fix your data")
    else:
        raise ValueError("utilities:fid: idd should be a string to an excel file (Extension should be include)")





    #check plot
    if not isinstance(plot , (list,np.ndarray, tuple)):
        raise ValueError('plot can only be a list or tuple or numpy.array')
    if np.array(plot).size != 2:
        raise ValueError("plot must be have a length equal to 2")

    #Since no error, then you can visualise data
    if plot[0]:
        plt.figure(dpi=128)
        plt.subplot(2, 1, 1)
        plt.plot(t, y, 'b-')
        plt.title("Identification Data")
        plt.ylabel('output')
        plt.grid(True, axis='both', which='both')

        plt.subplot(2, 1, 2)
        plt.plot(t, u, 'r-')
        plt.xlabel('time (sec)')
        plt.ylabel('input')
        plt.grid(True, axis='both', which='both')
        plt.show()

    if cleanDelay[0]:
        opti.findDelay = False
        opti.G.dt = 0
        truncy = np.nonzero(cleanDelay[1] <= t)
        y = y[truncy]
        u = u[truncy]
        t = t[truncy]- cleanDelay[1]

    if plot[0]:
        plt.figure(dpi=128)
        plt.subplot(2, 1, 1)
        plt.plot(t, y, 'b-')
        plt.title("Identification Data without delay")
        plt.ylabel('output')
        plt.grid(True, axis='both', which='both')

        plt.subplot(2, 1, 2)
        plt.plot(t, u, 'r-')
        plt.xlabel('time (sec)')
        plt.ylabel('input')
        plt.grid(True, axis='both', which='both')
        plt.show()

    if isinstance(vidd,str) and (vidd[-4:] =='xlsx' or vidd[-3:] =='xls'):
        data = pd.read_excel(vidd)
        vy= np.array(data.y.values)
        vu = np.array(data.u.values)
        vt = np.array(data.t.values)
        # vdt = vt[1]-vt[0]
        del data  # to save memory

        if vy.size != vu.size or vu.size != vt.size or vt.size != vy.size:
            raise IOError("utilies.fid: size of data vidd are not the same. Kindly Fix your data")
    else:
        raise ValueError("utilities:fid: vidd should be a string to an excel file (Extension should be include)")



    # Since no error in validation data, then you can visualise data
    if plot[1]:
        plt.figure(dpi=128)
        plt.subplot(2, 1, 1)
        plt.plot(vt, vy, 'b-')
        plt.title("Verification Data")
        plt.ylabel('output')
        plt.grid(True, axis='both', which='both')

        plt.subplot(2, 1, 2)
        plt.plot(vt, vu, 'r-')
        plt.xlabel('time (sec)')
        plt.ylabel('input')
        plt.grid(True, axis='both', which='both')
        plt.show()


    #check limits
    # check that limit must be a  2 * 2 array
    if limits is not None and np.array(limits).shape != (2, 2):
        raise ValueError("limits must a 2 by 2 matrix")

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
        return fidOutput(opti.G, y, u, t, vy, vu, vt)
    elif opti.polyFix[0] == 0 and opti.polyFix[1] == 1: #Poles polynomial is fixed i.e not optimized
        #check which to optimize coefficients, exponents or both
        if opti.optiFix == optFix.Free:
            #initial Guess
            x0 = np.append(xnum,xnnum)
            #limits
            lb = np.append(clim[0]*np.ones_like(xnum), elim[0]*np.ones_like(xnnum))
            if elim[0] < EXP_LB:
                lb[xnum.size:-1] = EXP_LB
            lb[-1] = 0
            ub = np.append(clim[1] * np.ones_like(xnum), elim[1] * np.ones_like(xnnum))
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
            #so that last power of  's' is always zero
            if elim[0] < EXP_LB:
                lb[:-1] = EXP_LB
            lb[-1] = 0
            ub = elim[1] * np.ones_like(x0)

        opti.funcFix = np.array([x0.size, 0])

    elif opti.polyFix[0] == 1 and opti.polyFix[1] == 0:#Zeroes polynomial is fixed i.e not optimized
        # check which to optimize coefficients, exponents or both
        if opti.optiFix == optFix.Free:
            # initial Guess
            x0 = np.append(xden, xnden)
            # limits
            lb = np.append(clim[0] * np.ones_like(xden), elim[0] * np.ones_like(xnden))
            if elim[0] < EXP_LB:
                lb[xden.size:-1] = EXP_LB
            lb[-1] = 0
            ub = np.append(clim[1] * np.ones_like(xden), elim[1] * np.ones_like(xnden))
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
            if elim[0] < EXP_LB:
                lb[:-1] = EXP_LB
            lb[-1] = 0
            ub = elim[1] * np.ones_like(x0)

        opti.funcFix = np.array([0, x0.size])

    else: #No fix, optimise Zero Polynomials and Poles Polynomials
        if opti.optiFix == optFix.Free:
            # initial Guess
            x0 = np.concatenate([xnum, xnnum, xden, xnden])
            # limits
            lb = np.concatenate([clim[0] * np.ones_like(xnum), elim[0] * np.ones_like(xnnum), clim[0] * np.ones_like(xden), elim[0] * np.ones_like(xnden)])
            if elim[0] < EXP_LB:
                lb[xnum.size: (xnum.size * 2)-1] = EXP_LB
                lb[(xnum.size * 2) + xden.size : -1] = EXP_LB
            lb[(xnum.size * 2)-1] = 0
            lb[-1] = 0
            ub = np.concatenate([clim[1] * np.ones_like(xnum), elim[1] * np.ones_like(xnnum), clim[1] * np.ones_like(xden), elim[1] * np.ones_like(xnden)])

            opti.funcFix = np.array([xnum.size * 2, xden.size * 2])
        elif opti.optiFix == optFix.Exp:
            # Fix Exponent, optimize Coefficient
            x0 = np.append(xnum, xden)
            # limits
            lb = clim[0] * np.ones_like(x0)
            ub = clim[1] * np.ones_like(x0)
            opti.funcFix = np.array([xnum.size, xden.size])
        elif opti.optiFix == optFix.Coeff:
            # Fix Coefficient, optimize Exponent
            x0 = np.append(xnnum, xnden)
            # limits
            lb = elim[0] * np.ones_like(x0)
            if elim[0] < EXP_LB:
                lb[:xnnum.size-1] = EXP_LB
                lb[xnnum.size:-1] = EXP_LB
            lb[xnnum.size-1] = 0
            lb[-1] = 0
            ub = elim[1] * np.ones_like(x0)
            opti.funcFix = np.array([xnum.size, xden.size])

    # Account for delay optimization and bounds
    if opti.findDelay:
        x0 = np.append(x0, xdt)
        lb = np.append(lb,0)
        ub = np.append(ub,10) #delay should only be positive thus use upper bound

    print("Please wait,  System Identification in progress...")

    start = datetime.now()
    #Run the identification Algorithm
    if opti.alg == optAlgo.LevenbergMarquardt:
        #bounds will not be used
        print("Bounds will not be applied with Levenberg Marquardts optimization algorithm")
        res = least_squares(_fracidfun, x0, args=(y,u,t,opti), ftol=1e-10, verbose=2, method='lm',max_nfev = 1000)
    elif opti.alg == optAlgo.TrustRegionReflective:
        res = least_squares(_fracidfun, x0, args=(y,u,t,opti), verbose=2, method='trf', bounds=(lb,ub), ftol=1e-10,max_nfev = 1000)
    elif opti.alg == optAlgo.Softl1:
        res = least_squares(_fracidfun, x0, args=(y, u, t, opti), verbose=2, bounds=(lb,ub), loss = 'soft_l1',method='trf', max_nfev = 1000)
    elif opti.alg == optAlgo.RobustLoss:
        res = least_squares(_fracidfun, x0, args=(y, u, t, opti), verbose=2, bounds=(lb,ub), loss = 'cauchy', method='trf',f_scale = 0.1, max_nfev = 1000)

    print('Time: ',datetime.now() - start)

    #Display Resuls
    # print(res.x,'\n\n', res.message, '\n\n', res.success, '\n\n', 'Number of Function Calls : {}'.format(res.nfev), '\n\n', 'Function Residuals: {}'.format(res.fun),'\n\n', 'Value at Solution: {}'.format(res.cost))

    #Xtract Delay first
    if opti.findDelay:
        xdt = res.x[-1]
        res.x = res.x[:-1]


    if opti.polyFix[0] == 0 and opti.polyFix[1] == 1:  # Poles polynomial is fixed i.e not optimized
        # so update zeros with optimized result
        if opti.optiFix == optFix.Free:
            xnum = res.x[0:xnum.size]
            xnnum = res.x[xnum.size:]
        elif opti.optiFix == optFix.Exp:
            xnum = res.x
        elif opti.optiFix == optFix.Coeff:
            # Fix Coefficient, update Exponent
            xnnum = res.x
        opti.funcFix = np.array([x0.size, 0])

    elif opti.polyFix[0] == 1 and opti.polyFix[1] == 0:  # Zeroes polynomial is fixed i.e not optimized
        # update optimize coefficients in poles
        if opti.optiFix == optFix.Free:
            xden = res.x[0:xden.size]
            xnden = res.x[xden.size:]
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

    IdentifiedG = FOTransFunc(xnum, xnnum, xden, xnden, xdt)
    lsimG = lsim(IdentifiedG,u,t)
    lsimvG = lsim(IdentifiedG,vu,vt)


    if plotid[0]:

        # plot identified system output vs Data from initial system
        plt.figure(dpi=128)
        plt.subplot(2, 1, 1)
        plt.plot(t, y, 'b-', t, lsimG, 'g-')
        plt.title("Identification Data vs Identified System")
        plt.ylabel('output')
        plt.legend(['iddata', 'idsystem'], loc='upper right')
        plt.grid(True, axis='both', which='both')

        # Fitness measure
        erro= y - lsimG
        fitness = 100 * (1 - (np.linalg.norm(erro) / np.linalg.norm(y - np.mean(lsimG))))

        plt.subplot(2, 1, 2)
        plt.plot(t, y-lsimG, 'r-')
        plt.title("Identified System error. fitness: {}%".format(round(fitness,2)))
        plt.xlabel('time (sec)')
        plt.ylabel('error')
        plt.grid(True, axis='both', which='both')
        plt.show()


    if plotid[1]:
        # plot identified system output vs Data from Verification system
        plt.figure(dpi=128)
        plt.subplot(2, 1, 1)
        plt.plot(vt, vy, 'b-', vt, lsimvG, 'g-')
        plt.title("Verification Data vs Identified System")
        plt.ylabel('output')
        plt.legend(['vdata', 'idsystem'], loc='upper left')
        plt.grid(True, axis='both', which='both')

        # Fitness measure
        erro = vy - lsimvG
        fitness = 100 * (1 - (np.linalg.norm(erro) / np.linalg.norm(vy - np.mean(lsimvG))))

        plt.subplot(2, 1, 2)
        plt.plot(vt,  vy-lsimvG, 'r-')
        plt.title("validated System error. fitness: {}%".format(round(fitness,2)))
        plt.xlabel('time (sec)')
        plt.ylabel('error')
        plt.grid(True, axis='both', which='both')
        plt.show()

    return fidOutput(IdentifiedG, y, u, t, vy, vu, vt)


class fidOutput():

    def __init__(self,G,y,u,t, vy,vu,vt):
        """
        :param G:   Identified system
        :type G:    FOTransFunc
        :param y:   Identification data output
        :type y:    numpy.ndarray
        :param u:   Identification data input
        :type u:    numpy.ndarray
        :param t:   Identification data time variable
        :type t:    numpy.ndarray
        :param vy:  Verification data output
        :type vy:   numpy.ndarray
        :param vu:  Verification data input
        :type vu:   numpy.ndarray
        :param vt:  Verification data time variable
        :type vt:   numpy.ndarray
        """
        self.G = G
        self.y = y
        self.u = u
        self.t = t
        self.vy = vy
        self.vu = vu
        self.vt = vt