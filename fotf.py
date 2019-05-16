#!/usr/bin/env python3
"""
FOTF object class and method definition
derived from the FOTF/FOMCON Matlab toolbox And Control Toolbox Python
---------------------------------------
Ported to python by Tobechukwu Onyedi
Revision date: 12th April 2019
"""

# External function declarations
import numpy as np
from numpy import (angle, array, empty, finfo, ndarray, ones,
                   polyadd, polymul, polyval, roots, sqrt, zeros, squeeze, exp, pi,
                   where, delete, real, poly, nonzero)
import scipy as sp
from numpy.polynomial.polynomial import polyfromroots
from scipy.signal import lti, tf2zpk, zpk2tf, cont2discrete
from copy import deepcopy
from warnings import warn
import inspect
from itertools import chain
from control.matlab import *
from lti import LTI
from numpy import linalg as LA
from matplotlib import pyplot as plt
from matplotlib import axes

__all__ = ['FOTransFunc', 'fotf', 'ss2tf', 'tfdata', 'poly2str', 'str2poly', 'freqresp', 'dcgain', 'newfotf', 'lsim', 'step', 'impulse']
MIN_COMM_ORDER = 0.01


# noinspection SpellCheckingInspection
class FOTransFunc(LTI):

    def __init__(self, *args):
        """TransferFunction(num, den[, dt])

        Construct a transfer function.

        The default constructor is Fractional order TransferFunction(num, den, nnum,  nden), where num and
        den are lists of lists of arrays containing polynomial coefficients. nnum and nden
        are coefficients of the  power of s.

        To create a discrete time Fractional Order transfer funtion, use TransferFunction(num,
        den, dt,nnum, nden) where 'dt' is the sampling time (or True for unspecified
        sampling time).  To call the copy constructor, call
        TransferFunction(sys), where sys is a Fractional Order TransferFunction object
        (continuous or discrete).

        """
        _args = deepcopy(args)
        epsi = 0.01

        # internal parameters, used to manipulate initial values from user
        _num, _nnum, _den, _nden, _dt = None, None, None, None, None
        if len(_args) == (0 or None):
            pass
        elif len(_args) >= 4:
            [_num, _nnum, _den, _nden] = args[0:4]
            if len(_args) == 5 and isinstance(_args[4], (float,int)):
                if _args[4] > 0:
                    _dt = _args[4]
                else:
                    _dt = 0

        elif len(_args) == 1 and isinstance(_args[0],(float,int)):
            _num = _args[0]
            _nnum = 0
            _den = 1.
            _nden = 0
            _dt = 0

        elif len(_args) == 1 and (_args[0] == 's'):
            _num = 1.
            _nnum = 1.
            _den = 1.
            _nden = 0.
            _dt = 0

        elif len(_args) == 1 and isinstance((_args[0]), FOTransFunc):
            [_num, _nnum, _den, _nden] = _args[0:4]
            if len(_args) == 5:
                _dt = _args[4]

        elif len(_args) >= 2 and isinstance(_args[0], list) and isinstance(_args[1], list):
            _num = _args[0][0]
            _nnum = _args[0][1]
            _den = _args[1][0]
            _nden = _args[1][1]

            if len(_args) >= 3 and isinstance(_args[2], (float,int)):
                if _args[2] > 0:
                    _dt = _args[2]
            else:
                _dt = 0

        elif len(_args) >= 2 and isinstance(_args[0], str) and isinstance(_args[1], str):
            if len(_args)>= 3 and isinstance(_args[2], (float,int)):
                if _args[2] > 0:
                    _dt = _args[2]
                    bas = 'z'
            else:
                _dt = 0
                bas = 's'

            _num,_nnum = str2poly(_args[0], bas)
            _den,_nden = str2poly(_args[1], bas)
        else:
            raise ValueError("fotf.FOTransFunc: Needs 1, 2 , 3 ,4 or 5 string or int or list or ndarray arguments. received {}.".format(len(args)))

        _num = _clean_part(_num)
        _den = _clean_part(_den)
        _nnum = _clean_part(_nnum)
        _nden = _clean_part(_nden)

        inputs = len(_num[0])
        outputs = len(_num)

        # Make sure numerator and denominator matrices have consistent sizes
        if inputs != len(_den[0]):
            raise ValueError("The numerator has %i input(s), but the denominator has "
                             "%i\ninput(s)." % (inputs, len(_den[0])))
        if outputs != len(_den):
            raise ValueError("The numerator has %i output(s), but the denominator has "
                             "%i\noutput(s)." % (outputs, len(_den)))

        # Additional checks/updates on structure of the transfer function
        for i in range(outputs):
            # Make sure that each row has the same number of columns
            if len(_num[i]) != inputs:
                raise ValueError("Row 0 of the numerator matrix has %i elements, but row "
                                 "%i\nhas %i." % (inputs, i, len(_num[i])))
            if len(_den[i]) != inputs:
                raise ValueError("Row 0 of the denominator matrix has %i elements, but row "
                                 "%i\nhas %i." % (inputs, i, len(_den[i])))

            # Check for zeros in numerator or denominator
            # TODO: Right now these checks are only done during construction.
            # It might be worthwhile to think of a way to perform checks if the
            # user modifies the transfer function after construction.
            for j in range(inputs):
                # Check that we don't have any zero denominators.
                zeroden = True
                for k in _den[i][j]:
                    if k:
                        zeroden = False
                        break
                if zeroden:
                    raise ValueError("Input %i, output %i has a zero denominator."
                                     % (j + 1, i + 1))

                # If we have zero numerators, set the denominator to 1.
                zeronum = True
                for k in _num[i][j]:
                    if k:
                        zeronum = False
                        break
                if zeronum:
                    _den[i][j] = ones(1)

        LTI.__init__(self, inputs, outputs, dt=_dt)
        self.num = _num
        self.den = _den
        self.nnum = _nnum
        self.nden = _nden
        self.epsi = epsi
        # self.dt = _dt
        self._truncatecoeff()

    def __call__(self, s):
        """Evaluate the system's transfer function for a complex variable

        For a SISO transfer function, returns the value of the
        transfer function.  For a MIMO transfer fuction, returns a
        matrix of values evaluated at complex variable s."""

        if self.issiso():
            # return a scalar
            return self.horner(s)[0][0]
        else:
            # return a matrix
            return self.horner(s)

    def step(self, t, output=True, plot=True):
        """
        :param t: list, ndarray of in secs
        :type t: list, ndarray
        :param output: bool, used to determine if you want output
        :param plot: bool, used to indicate you want plot

        :return t, y: Time (t) and Linear Simulation (y)

        """
        if t is None:
            t = step_auto_range(self)
        elif not isinstance(t, ndarray):
            t = np.array(t)
        u = np.ones(t.size)
        y = lsim(self, u, t, plot=False)

        if output is True and plot is True:
            plt.figure()
            plt.plot(t, y)
            plt.title('Step response')
            plt.xlabel('Time [s]')
            plt.ylabel('Amplitude')
            plt.grid()
            plt.show()

            # Plot final value if present
            # Test DC gain
            myGain = dcgain(self)
            if np.isinf(myGain) or (np.abs(myGain) < np.finfo(float).resolution):
                pass
            else:
                plt.figure()
                plt.plot([t[0], t[-1]], [myGain, myGain], ':k')
                plt.title('Dc Gain')
                plt.xlabel('Time [s]')
                plt.ylabel('Amplitude')
                plt.grid()
                plt.show()
            return t, y
        elif output is False and plot is True:
            plt.figure()
            plt.plot(t, y)
            plt.title('Step response')
            plt.xlabel('Time [s]')
            plt.ylabel('Amplitude')
            plt.grid()
            plt.show()

            # Plot final value if present
            # Test DC gain
            myGain = dcgain(self)
            if np.isinf(myGain) or (np.abs(myGain) < np.finfo(float).resolution):
                pass
            else:
                plt.figure()
                plt.plot([t[0], t[-1]], [myGain, myGain], ':k')
                plt.title('Dc Gain')
                plt.xlabel('Time [s]')
                plt.ylabel('Amplitude')
                plt.grid()
                plt.show()
        elif output is True and plot is False:
            return y
        else:
            raise ValueError("fotf.step: Wrong input for keyword 'output' or 'plot', one must be True")

    def _truncatecoeff(self):
        """Remove extraneous zero coefficients from num and den.

        Check every element of the numerator and denominator matrices, and
        truncate leading zeros.  For instance, running self._truncatecoeff()
        will reduce self.num = [[[0, 0, 1, 2]]] to [[[1, 2]]].

        """

        # Beware: this is a shallow copy.  This should be okay.
        data = [self.num, self.den]
        for p in range(len(data)):
            for i in range(self.outputs):
                for j in range(self.inputs):
                    # Find the first nontrivial coefficient.
                    nonzero = None
                    for k in range(data[p][i][j].size):
                        if data[p][i][j][k]:
                            nonzero = k
                            break

                    if nonzero is None:
                        # The array is all zeros.
                        data[p][i][j] = zeros(1)
                    else:
                        # Truncate the trivial coefficients.
                        data[p][i][j] = data[p][i][j][nonzero:]
        [self.num, self.den] = data

        ndata = [self.nnum, self.nden]
        for p in range(len(data)):
            for i in range(self.outputs):
                for j in range(self.inputs):
                    # Find the first nontrivial coefficient.
                    nonzero = None
                    for k in range(data[p][i][j].size):
                        if data[p][i][j][k]:
                            nonzero = k
                            break

                    if nonzero is None:
                        # The array is all zeros.
                        data[p][i][j] = zeros(1)
                    else:
                        # Truncate the trivial coefficients.
                        data[p][i][j] = data[p][i][j][nonzero:]
        [self.nnum, self.nden] = ndata

    def isstable(self, doPlot=True):
        """

        :param doPlot: if set to 'True', will create a plot with relevant system poles.
                        Default is 'False'
        :type doPlot:   bool
        :return K:      K - bool, 'True' if system is stable, else 'False'
                        q - calculated commensurate-order
                        err - stability assessment error norm
                        apol - closest poles' absolute imaginary part value to the unstable region
        """
        comm_factor = MIN_COMM_ORDER ** -1
        if isinstance(self, FOTransFunc):
            b = self.den[0][0]
            nb = self.nden[0][0]
            nb1 = nb * comm_factor
            q = comm_order(self, 'den')
            newnb = np.array(nb1 / (comm_factor * q), dtype=np.int32)

            c = np.zeros(newnb[0] + 1)
            c[newnb] = b
            cslice = c[::-1]
            p = np.roots(cslice)

            if p is not None:
                absp = np.array(p * (np.abs(p) > np.finfo(float).resolution))

            err = np.linalg.norm(polyval(cslice, absp))
            apol = np.amin(np.abs(np.angle(absp)))
            K = apol > q * np.pi * 0.5

            # Check if drawing is requested
            if doPlot:

                # create new figure
                x = plt.figure(dpi=512)
                axes.Axes.set_autoscale_on(x, True)
                plt.plot(np.real(p), np.imag(p), 'x', 0, 0, '+')

                # Get and check x axis limit
                left, right = plt.xlim()
                # left = np.imag(p).min()
                # right = np.imag(p).max()
                if right <= 0:
                    right = abs(left)
                    plt.xlim(left, right)

                left = 0
                gpi = right * np.tan(q * np.pi * 0.5)

                x_fill = np.array([left, right, right, left])
                y_fill = np.array([0, gpi, -gpi, 0])
                plt.fill_between(x_fill, y_fill, color='red')
                plt.show()
            else:
                pass

        return [K, q, err, apol]

    def __str__(self, var=None):
        """String representation of the FRACTIONAL Order transfer function."""

        mimo = self.inputs > 1 or self.outputs > 1
        if var is None:
            # ! TODO: replace with standard calls to lti functions
            if self.dt is None or self.dt == 0:
                var = 's'
            else:
                var = 'z'
        outstr = ""

        for i in range(self.inputs):
            for j in range(self.outputs):
                if mimo:
                    outstr += "\nInput {0} to output {1}".format(i + 1, j + 1)

                # Convert the numerator and denominator polynomials to strings.
                numstr = poly2str(self.num[j][i], self.nnum[j][i], var=var)
                denstr = poly2str(self.den[j][i], self.nden[j][i], var=var)

                # Figure out the length of the separating line
                dashcount = max(len(numstr), len(denstr))
                dashes = '-' * dashcount

                # Center the numerator or denominator
                if len(numstr) < dashcount:
                    numstr = (' ' * int(round((dashcount - len(numstr)) / 2)) +
                              numstr)
                if len(denstr) < dashcount:
                    denstr = (' ' * int(round((dashcount - len(denstr)) / 2)) +
                              denstr)

                outstr += "\n" + numstr + "\n" + dashes + "\n" + denstr + "\n"

        # See if this is a discrete time system with specific sampling time
        if not (self.dt is None) and type(self.dt) != bool and self.dt > 0:
            # TODO: replace with standard calls to lti functions if possible
            outstr += "\ndt = " + self.dt.__str__() + "\n"

        return outstr

    # represent a string, makes display work for IPython
    def __repr__(self, ):
        return self.__str__()

    # still investigating on this
    def _repr_latex_(self, var=None):
        """LaTeX representation of the transfer function, for Jupyter notebook"""

        mimo = self.inputs > 1 or self.outputs > 1

        if var is None:
            # ! TODO: replace with standard calls to lti functions
            var = 's' if self.dt is None or self.dt == 0 else 'z'

        out = ['$$']

        if mimo:
            out.append(r"\begin{bmatrix}")

        for i in range(self.outputs):
            for j in range(self.inputs):
                # Convert the numerator and denominator polynomials to strings.
                numstr = poly2str(self.num[i][j], self.nnum, var=var)
                denstr = poly2str(self.den[i][j], self.nnum, var=var)

                out += [r"\frac{", numstr, "}{", denstr, "}"]

                if mimo and j < self.outputs - 1:
                    out.append("&")

            if mimo:
                out.append(r"\\")

        if mimo:
            out.append(r" \end{bmatrix}")

        # See if this is a discrete time system with specific sampling time
        if not (self.dt is None) and type(self.dt) != bool and self.dt > 0:
            out += ["\quad dt = ", str(self.dt)]

        out.append("$$")

        return ''.join(out)

    def __neg__(self):
        """Negate a Fractional order transfer function."""

        _num = deepcopy(self.num)
        for i in range(self.outputs):
            for j in range(self.inputs):
                _num[i][j] *= -1

        return FOTransFunc(_num, self.nnum, self.den, self.nden, self.dt)

    # Not done this yet
    def __add__(self, other):
        """Add two FOTransfunc objects ."""

        # Convert the second argument to a transfer function.
        if not isinstance(other, FOTransFunc):
            _other = fotf(other)
        else:
            _other = deepcopy(other)

        # Check that the input-output sizes are consistent.
        if self.inputs != other.inputs:
            raise ValueError("The first summand has %i input(s), but the second has %i."
                             % (self.inputs, other.inputs))
        if self.outputs != other.outputs:
            raise ValueError("The first summand has %i output(s), but the second has %i."
                             % (self.outputs, other.outputs))

        # Figure out the sampling time to use
        if self.dt == other.dt:
            aa = self.den[0][0]
            othera = _other.den[0][0]
            bb = self.num[0][0]
            otherb = _other.num[0][0]

            #change denominator shape to 2d
            aa = np.reshape(aa, (1, aa.shape[0]))
            othera = np.reshape(othera, (1, othera.shape[0]))
            # change numerator shape to 2d
            bb = np.reshape(bb, (1, bb.shape[0]))
            otherb = np.reshape(otherb, (1, otherb.shape[0]))

            #use Kron product
            a = sp.linalg.kron(aa, othera)
            b0 = sp.linalg.kron(aa, otherb)
            b1 = sp.linalg.kron(bb, othera)

            # revert shapes
            a = np.reshape(a, a.size)
            b0 = np.reshape(b0, b0.size)
            b1 = np.reshape(b1, b1.size)
            b = np.concatenate((b0,b1))

            na = np.empty(1)
            nb = np.empty(1)

            for i in self.nden[0][0]:
                na = np.append(na, (i + _other.nden[0][0]))
                nb = np.append(nb, (i+ _other.nnum[0][0]))
            na = np.delete(na, 0)

            for j in self.nnum[0][0]:
                nb = np.append(nb, (j + _other.nden[0][0]))
            nb = np.delete(nb, 0)

            return simple(fotf(b,nb,a,na, self.dt))

        else:
            raise ValueError("Cannot handle different delay times")

        # # Preallocate the numerator and denominator of the sum.
        # num = [[[] for j in range(self.inputs)] for i in range(self.outputs)]
        # den = [[[] for j in range(self.inputs)] for i in range(self.outputs)]
        #
        # for i in range(self.outputs):
        #     for j in range(self.inputs):
        #         num[i][j], den[i][j] = _add_siso(self.num[i][j], self.den[i][j],
        #                                          other.num[i][j],
        #                                          other.den[i][j])
        #
        # return FOTransFunc(num, den, dt)

    def __radd__(self, other):
        """Right add two LTI objects (parallel connection)."""
        return self + other

    def __sub__(self, other):
        """Subtract two LTI objects."""
        return self + (-other)

    def __rsub__(self, other):
        """Right subtract two LTI objects."""
        return other + (-self)

    def __mul__(self, other):
        """Multiply two FOTransFunc Object (serial connection)."""
        # Convert the second argument to a transfer function.
        if not isinstance(other, FOTransFunc):
            _other = fotf(other)
        else:
            _other = other

        aa = self.num[0][0]
        othera = _other.num[0][0]
        bb = self.den[0][0]
        otherb = _other.den[0][0]

        aa = np.reshape(aa,(1,aa.shape[0]))
        othera = np.reshape(othera,(1,othera.shape[0]))

        bb = np.reshape(bb,(1,bb.shape[0]))
        otherb = np.reshape(otherb, (1, otherb.shape[0]))

        a = sp.linalg.kron(aa,othera)
        b = sp.linalg.kron(bb,otherb)

        #revert shapes
        a = np.reshape(a, a.size)
        b = np.reshape(b, b.size)

        na = np.empty(1)
        nb = np.empty(1)

        for i in self.nnum[0][0]:
            na = np.append(na, (i + _other.nnum[0][0]))
        na = np.delete(na,0)

        for j in self.nden[0][0]:
            nb = np.append(nb, (j + _other.nden[0][0]))
        nb = np.delete(nb, 0)

        return simple(fotf(a,na,b,nb, self.dt + _other.dt))

    def __rmul__(self, other):
        """Right multiply two FOTransFunc objects (serial connection)."""
        if not isinstance(other, FOTransFunc):
            _other = fotf(other)
        else:
            _other = other

        aa = self.num[0][0]
        othera = _other.num[0][0]
        bb = self.den[0][0]
        otherb = _other.den[0][0]

        aa = np.reshape(aa,(1,aa.shape[0]))
        othera = np.reshape(othera,(1,othera.shape[0]))

        bb = np.reshape(bb,(1,bb.shape[0]))
        otherb = np.reshape(otherb, (1, otherb.shape[0]))

        a = sp.linalg.kron(othera,aa)
        b = sp.linalg.kron(otherb, bb)

        #revert shapes
        a = np.reshape(a, a.size)
        b = np.reshape(b, b.size)

        na = np.empty(1)
        nb = np.empty(1)

        for i in self.nnum[0][0]:
            na = np.append(na, (i + _other.nnum[0][0]))
        na = np.delete(na,0)

        for j in self.nden[0][0]:
            nb = np.append(nb, (j + _other.nden[0][0]))
        nb = np.delete(nb, 0)

        return simple(fotf(a,na,b,nb, self.dt + _other.dt))
        # # Convert the second argument to a transfer function.
        # if isinstance(other, (int, float, complex, np.number)):
        #     other = _convert_to_transfer_function(other, inputs=self.inputs,
        #                                           outputs=self.inputs)
        # else:
        #     other = _convert_to_transfer_function(other)
        #
        # # Check that the input-output sizes are consistent.
        # if other.inputs != self.outputs:
        #     raise ValueError("C = A * B: A has %i column(s) (input(s)), but B has %i "
        #                      "row(s)\n(output(s))." % (other.inputs, self.outputs))
        #
        # inputs = self.inputs
        # outputs = other.outputs
        #
        # # Figure out the sampling time to use
        # if self.dt is None and other.dt is not None:
        #     dt = other.dt  # use dt from second argument
        # elif (other.dt is None and self.dt is not None) \
        #         or (self.dt == other.dt):
        #     dt = self.dt  # use dt from first argument
        # else:
        #     raise ValueError("Systems have different sampling times")
        #
        # # Preallocate the numerator and denominator of the sum.
        # num = [[[0] for j in range(inputs)] for i in range(outputs)]
        # den = [[[1] for j in range(inputs)] for i in range(outputs)]
        #
        # # Temporary storage for the summands needed to find the
        # # (i, j)th element
        # # of the product.
        # num_summand = [[] for k in range(other.inputs)]
        # den_summand = [[] for k in range(other.inputs)]
        #
        # for i in range(outputs):  # Iterate through rows of product.
        #     for j in range(inputs):  # Iterate through columns of product.
        #         for k in range(other.inputs):  # Multiply & add.
        #             num_summand[k] = polymul(other.num[i][k], self.num[k][j])
        #             den_summand[k] = polymul(other.den[i][k], self.den[k][j])
        #             num[i][j], den[i][j] = _add_siso(
        #                 num[i][j], den[i][j],
        #                 num_summand[k], den_summand[k])

        return FOTransFunc(num, den, dt)

    def __truediv__(self, other):
        """Divide two FOTransFunc objects.
        Right division of fractional-order dynamic systems.
        Note: if delays in the systems are present, dividing the two systems may
        result in positive delays thus the overall delay of the system
        will be changed to zero. There with be a warning.

        """


        if not isinstance(other, FOTransFunc):
            other = newfotf(other)

        # Figure out the sampling time to use
        if self.dt is None and other.dt is not None:
            dt = other.dt  # use dt from second argument
        elif other.dt is None and self.dt is not None:
            dt = self.dt  # use dt from first argument
        elif other.dt is not None and self.dt is not None:
            dt = self.dt - other.dt
            if dt < 0:
                dt = 0
                warn("FOTTransFunc.__truedivide__: Resulting FOTransFunc has positive delay: changing to zero")
        else:
            raise ValueError("Systems have different sampling times")

        g = self * other.__invert__()
        g.dt = dt

        return g

    def __rtruediv__(self, other):
        """Right divide two FOTransfunc objects."""
        if not isinstance(other, FOTransFunc):
            other = newfotf(other)
        return other / self

    def __rdiv__(self, other):
        return FOTransFunc.__rtruediv__(self, other)

    def __pow__(self, n):
        """
        Computes the self multiplicaiton of object by n-times

        :param n: int
        :return: self**n

        """
        num,nnum,den,nden,dt = fotfparam(self)
        if not isinstance(n,int):
            raise ValueError("Exponent must be an integer")
        elif den.size * num.size == 1 and nden == 0 and nnum ==1:
            return FOTransFunc(1,0,1,n)  # unity
        else:
            if n >=0:
                y = 1
                for i in range(n):
                    y *=self
            else:
                a = 1
                for i in range(n):
                    a *=self
                    y = a.__invert__()
            if not isinstance(y,FOTransFunc):
                y = fotf(y)
            update = simple(y)
            update.dt = n * self.dt
        return update

    def __eq__(self,other):
        """
        checks is a transfer function is equal, overides the '==' operator
        :param other: A fractional order Transfer function object
        :type other: FOTransFunc
        :return: True or False
        """

        if not isinstance(other,FOTransFunc):
            other = FOTransFunc(other)

        num, nnum,den, nden, dt = fotfparam(self)
        othernum, othernnum, otherden, othernden , otherdt = fotfparam(other)
        if num.all() == othernum.all() and nnum.all() ==othernnum.all() and den.all() ==otherden.all() \
        and nden.all() == othernden.all() and dt == otherdt:
            return True
        else:
            False



    def __invert__(self):
        num, nnum, den,nden, dt = fotfparam(self)
        return fotf(den,nden,num,nnum,-dt)

    def __getitem__(self, key):
        key1, key2 = key

        # pre-process
        if isinstance(key1, int):
            key1 = slice(key1, key1 + 1, 1)
        if isinstance(key2, int):
            key2 = slice(key2, key2 + 1, 1)
        # dim1
        start1, stop1, step1 = key1.start, key1.stop, key1.step
        if step1 is None:
            step1 = 1
        if start1 is None:
            start1 = 0
        if stop1 is None:
            stop1 = len(self.num)
        # dim1
        start2, stop2, step2 = key2.start, key2.stop, key2.step
        if step2 is None:
            step2 = 1
        if start2 is None:
            start2 = 0
        if stop2 is None:
            stop2 = len(self.num[0])

        num = []
        den = []
        for i in range(start1, stop1, step1):
            num_i = []
            den_i = []
            for j in range(start2, stop2, step2):
                num_i.append(self.num[i][j])
                den_i.append(self.den[i][j])
            num.append(num_i)
            den.append(den_i)
        if self.isctime():
            return FOTransFunc(num, den)
        else:
            return FOTransFunc(num, den, self.dt)

    def evalfr(self, omega):
        """Evaluate a transfer function at a single angular frequency.

        self._evalfr(omega) returns the value of the transfer function
        matrix with input value s = i * omega.

        """
        warn("TransferFunction.evalfr(omega) will be deprecated in a "
             "future release of python-control; use evalfr(sys, omega*1j) "
             "instead", PendingDeprecationWarning)
        return self._evalfr(omega)

    def _evalfr(self, omega):
        """Evaluate a transfer function at a single angular frequency."""
        # TODO: implement for discrete time systems
        if isdtime(self, strict=True):
            # Convert the frequency to discrete time
            dt = timebase(self)
            s = exp(1.j * omega * dt)
            if np.any(omega * dt > pi):
                warn("_evalfr: frequency evaluation above Nyquist frequency")
        else:
            s = 1.j * omega

        return self.horner(s)

    def horner(self, s):
        """Evaluate the systems's transfer function for a complex variable

        Returns a matrix of values evaluated at complex variable s.
        """

        # Preallocate the output.
        if getattr(s, '__iter__', False):
            out = empty((self.outputs, self.inputs, len(s)), dtype=complex)
        else:
            out = empty((self.outputs, self.inputs), dtype=complex)

        for i in range(self.outputs):
            for j in range(self.inputs):
                out[i][j] = (polyval(self.num[i][j], s) /
                             polyval(self.den[i][j], s))

        return out

    # Method for generating the frequency response of the system
    def freqresp(self, minExp=None, maxExp=None, numPts=1000):
        """Frequency response of fractional-order transfer functions.
            at which the frequency response must be computed.

            :param minExp: minimum exponent of Frequency to compute
            :type w: int, float
            :param maxExp: maximum exponent of Frequency to compute
            :type maxExp: int, float
            :param numPts: Number of points within interval minExp and maxExp
            :type maxExp: int

            :returns :Magnitude (Db), Phase in Degree, w -   Frequency

        """

        if minExp is None:
            minExp = -5
        if maxExp is None:
            maxExp = 5
        if numPts is None:
            numPts = 1000
        w = np.logspace(minExp, maxExp, numPts)

        _num, _nnum, _den, _nden, _dt = fotfparam(self)
        jj = 1j
        lenW = w.size
        r = np.zeros(lenW, dtype=complex)
        for k in range(lenW):
            bb = np.power((jj * w[k]), _nnum)
            aa = np.power((jj * w[k]), _nden)
            P = _num @ bb
            Q = _den @ aa
            r[k] = P / Q

        # Delay
        if self.dt > 0:
            for k2 in range(lenW):
                r[k2] *= np.exp(-jj * w(k2) * self.dt)

        rangle = unwrap(np.angle(r))
        rangleCalcDeg = np.rad2deg(rangle)
        rmagDb = 20 * np.log10(np.absolute(r))

        # # from control library, used basically to plot bode. But noticed an error so had to code it myself
        # H1 = frd(r,w)
        # rmagDb, rangledeg, w = bode(H1, w, dB=True, Plot=True, deg=True)
        # plt.show()

        plt.figure(dpi=128)
        plt.figure(1)
        plt.subplot(2, 1, 1)
        plt.semilogx(w, rmagDb, 'g-')
        plt.ylabel('Magnitude (Db)')
        plt.title('Bode Plot')
        plt.grid(True, axis='both', which='both')

        plt.subplot(2, 1, 2)
        plt.semilogx(w, rangleCalcDeg, 'g-')
        plt.xlabel('Frequency (rad/s)')
        plt.ylabel('Phase (deg)')
        plt.grid(True, axis='both', which='both')
        plt.show()

    def poles(self):
        """Computes the poles of a Fractional Order Transfer function.

        :returns absZeros, err: Zeros ,error in calculation of zeros"""
        if self.inputs > 1 or self.outputs > 1:
            raise NotImplementedError("FOTransFunc.poles is currently only implemented for SISO systems.")
        else:
            # for now, just give poles of a SISO tf
            comm_factor = MIN_COMM_ORDER ** -1
            b = self.den[0][0]
            nb = self.nden[0][0]
            nb1 = nb * comm_factor
            q = comm_order(self, 'den')
            newnb = np.array(nb1 / (comm_factor * q), dtype=np.int32)
            c = np.zeros(newnb[0] + 1)
            c[newnb] = b
            cslice = c[::-1]
            p = np.roots(cslice)
            if p is not None:
                #resolution Checker
                absZeros = np.array(p * (np.abs(p) > np.finfo(float).resolution))

            err = np.linalg.norm(polyval(cslice, absZeros))
            return absZeros, err

    def zeros(self):
        """Computes the zoles of a Fractional-Order Transfer function.

        :returns absZeros, err: poles ,error in calculation of poles"""
        if self.inputs > 1 or self.outputs > 1:
            raise NotImplementedError("FOTransFunc.zeros is currently only implemented for SISO systems.")
        else:
            # for now, just give zeros of a SISO tf
            comm_factor = MIN_COMM_ORDER ** -1
            a = self.num[0][0]
            na = self.num[0][0]
            na1 = na * comm_factor
            q = comm_order(self, 'den')
            newna = np.array(na1 / (comm_factor * q), dtype=np.int32)
            c = np.zeros(newna[0] + 1)
            c[newna] = a
            cslice = c[::-1]
            p = np.roots(cslice)
            if p is not None:
                # resolution Checker
                absZeros = np.array(p * (np.abs(p) > np.finfo(float).resolution))

            err = np.linalg.norm(polyval(cslice, absZeros))
            return absZeros, err

    def feedback(self, other, sign=-1):
        """Feedback connection of two input/output fractional-order transfer functions.
        M = G.feedback(H) computes a closed-loop model for the feedback loop:
        u --->O---->[ G ]---------> y
              ^               |
              |               |
              ------[ H ]<-----

        """
        if not isinstance(other, FOTransFunc):
            other = newfotf(other)
        if self.dt == other.dt:

            aa = self.den[0][0]
            othera = other.den[0][0]
            bb = self.num[0][0]
            otherb = other.num[0][0]

            # change denominator shape to 2d
            aa = np.reshape(aa, (1, aa.shape[0]))
            othera = np.reshape(othera, (1, othera.shape[0]))
            # change numerator shape to 2d
            bb = np.reshape(bb, (1, bb.shape[0]))
            otherb = np.reshape(otherb, (1, otherb.shape[0]))

            # use Kron product
            b = sp.linalg.kron(bb, othera)
            a0 = sp.linalg.kron(bb, otherb)
            a1 = sp.linalg.kron(aa, othera)

            # revert shapes
            b = np.reshape(b, b.size)
            a0 = np.reshape(a0, a0.size)
            a1 = np.reshape(a1, a1.size)
            a = np.concatenate((a0, a1))

            na = np.empty(1)
            nb = np.empty(1)

            for i in self.nnum[0][0]:
                nb = np.append(nb, (i + other.nden[0][0]))
                na = np.append(na, (i + other.nnum[0][0]))
            nb = np.delete(nb, 0)

            for j in self.nden[0][0]:
                na = np.append(na, (j + other.nden[0][0]))
            na = np.delete(na, 0)

            return simple(fotf(b, nb, a, na, self.dt))

        # For MIMO or SISO systems, the analytic expression is
        #     self / (1 - sign * other * self)
        # But this does not work correctly because the state size will be too
        # large.

    def returnScipySignalLTI(self):
        """Return a list of a list of scipy.signal.lti objects.

        For instance,

        >>> out = tfobject.returnScipySignalLTI()
        >>> out[3][5]

        is a signal.scipy.lti object corresponding to the
        transfer function from the 6th input to the 4th output.

        """

        # TODO: implement for discrete time systems
        if self.dt != 0 and self.dt is not None:
            raise NotImplementedError("Function not \
                    implemented in discrete time")

        # Preallocate the output.
        out = [[[] for j in range(self.inputs)] for i in range(self.outputs)]

        for i in range(self.outputs):
            for j in range(self.inputs):
                out[i][j] = lti(self.num[i][j], self.den[i][j])

        return out

    def sample(self, Ts, method='zoh', alpha=None):
        """Convert a continuous-time system to discrete time

        Creates a discrete-time system from a continuous-time system by
        sampling.  Multiple methods of conversion are supported.

        Parameters
        ----------
        Ts : float
            Sampling period
        method :  {"gbt", "bilinear", "euler", "backward_diff", "zoh", "matched"}
            Which method to use:

               * gbt: generalized bilinear transformation
               * bilinear: Tustin's approximation ("gbt" with alpha=0.5)
               * euler: Euler (or forward differencing) method ("gbt" with alpha=0)
               * backward_diff: Backwards differencing ("gbt" with alpha=1.0)
               * zoh: zero-order hold (default)

        alpha : float within [0, 1]
            The generalized bilinear transformation weighting parameter, which
            should only be specified with method="gbt", and is ignored otherwise

        Returns
        -------
        sysd : StateSpace system
            Discrete time system, with sampling rate Ts

        Notes
        -----
        1. Available only for SISO systems

        2. Uses the command `cont2discrete` from `scipy.signal`

        Examples
        --------
        >>> sys = FOTransFunc(1, [1,1])
        >>> sysd = sys.sample(0.5, method='bilinear')

        """
        if not self.isctime():
            raise ValueError("System must be continuous time system")
        if not self.issiso():
            raise NotImplementedError("MIMO implementation not available")
        if method == "matched":
            return _c2d_matched(self, Ts)
        sys = (self.num[0][0], self.den[0][0])
        numd, dend, dt = cont2discrete(sys, Ts, method, alpha)
        return FOTransFunc(numd[0, :], dend, dt)

    def dcgain(self):
        """Return the zero-frequency (or DC) gain

        For a continous-time transfer function G(s), the DC gain is G(0)
        For a discrete-time transfer function G(z), the DC gain is G(1)

        Returns
        -------
        gain : ndarray
            The zero-frequency gain
        """
        if self.isctime():
            return self._dcgain_cont()
        else:
            return self(1)

    def _dcgain_cont(self):
        """_dcgain_cont() -> DC gain as matrix or scalar

        Special cased evaluation at 0 for continuous-time systems."""
        gain = np.empty((self.outputs, self.inputs), dtype=float)
        for i in range(self.outputs):
            for j in range(self.inputs):
                num = self.num[i][j][-1]
                den = self.den[i][j][-1]
                if den:
                    gain[i][j] = num / den
                else:
                    if num:
                        # numerator nonzero: infinite gain
                        gain[i][j] = np.inf
                    else:
                        # numerator is zero too: give up
                        gain[i][j] = np.nan
        return np.squeeze(gain)

# c2d function contributed by Benjamin White, Oct 2012

def _c2d_matched(sysC, Ts):
    # Pole-zero match method of continuous to discrete time conversion
    szeros, spoles, sgain = tf2zpk(sysC.num[0][0], sysC.den[0][0])
    zzeros = [0] * len(szeros)
    zpoles = [0] * len(spoles)
    pregainnum = [0] * len(szeros)
    pregainden = [0] * len(spoles)
    for idx, s in enumerate(szeros):
        sTs = s * Ts
        z = exp(sTs)
        zzeros[idx] = z
        pregainnum[idx] = 1 - z
    for idx, s in enumerate(spoles):
        sTs = s * Ts
        z = exp(sTs)
        zpoles[idx] = z
        pregainden[idx] = 1 - z
    zgain = np.multiply.reduce(pregainnum) / np.multiply.reduce(pregainden)
    gain = sgain / zgain
    sysDnum, sysDden = zpk2tf(zzeros, zpoles, gain)
    return FOTransFunc(sysDnum, sysDden, Ts)


def poly2str(coeffs, powcoeffs, var='s'):
    """Converts a Fractional Order transfer function polynomial to a string

    Output:
        String

    Arguments:
        Coefficients, Power Coefficients and  string value of Transform 'z' or 's'[Optional].
        's' is default

    """

    if var == 's' or 'z' or None:
        pass
    else:
        raise ValueError("Input should be 's' or 'z' not {}.".format(var))

    thestr = ""
    # Compute the number of coefficients
    N = len(coeffs) - 1

    for k, coef in enumerate(coeffs):
        if k == 0:
            thestr += "{0:.2f}{1}^[{2:.2f}] ".format(coef, var, powcoeffs[k])
        else:
            thestr += "+ {0:.2f}{1}^[{2:.2f}] ".format(coef, var, powcoeffs[k])

    if var == 's' or 'z':
        thestr = thestr.replace('{}^[0.00]'.format(var), '')
        thestr = thestr.replace('1.00{}'.format(var), var)
    thestr = thestr.replace('^[0.00]', '')
    thestr = thestr.replace('+-', '-')
    thestr = thestr.replace('+ -', '- ')
    thestr = thestr.replace('^[1.00]', '')
    thestr = thestr.replace('.00', '')
    thestr = thestr.replace('[', '{')
    thestr = thestr.replace(']', '}')

    return thestr

def str2poly(polystr, bases=None):
    """Converts a string representation of a Fractional order transfer function to Polynomial
    represented by a list

    Args: Fractional order string, base i,e 's' or 'z'.
        default value is 's'
    """
    if isinstance(polystr, str):
        polystr = polystr.replace(" ", "")
        polystr = polystr.replace("{", "")
        polystr = polystr.replace("[", "")
        polystr = polystr.replace("}-", "+-")
        polystr = polystr.replace("]-", "+-")
        polystr = polystr.replace("}", "")
        polystr = polystr.replace("]", "")
        polystr = polystr.replace("*", "")

        if bases is None:
            bases = 's'
        polystr = polystr.split('+')
        for k in range(len(polystr)):
            polystr[k] = polystr[k].split('^')

        row = len(polystr)
        for i in range(row):
            for j in range(len(polystr[i])):
                if bases in polystr[i][j]:
                    polystr[i][j] = polystr[i][j].replace(bases, "")
                    if polystr[i][j] == '':
                        polystr[i][j] = 0
                    else:
                        polystr[i][j] = float(polystr[i][j])
                else:
                    polystr[i][j] = float(polystr[i][j])
                if len(polystr[i]) == 1:
                    polystr[i].append(0)

        outstring = []

        for i in range(2):
            polyb = []
            for j in range(row):
                polyb.append(polystr[j][i])
            outstring.append(polyb)

        return outstring[0], outstring[1]
    else:
        raise ValueError("Input should be of format 'str' with bases type 'z' or 's'")

# TODO:  Solve for FOTF Object
def _add_siso(num1, den1, num2, den2):
    """Return num/den = num1/den1 + num2/den2.

    Each numerator and denominator is a list of polynomial coefficients.

    """

    num = polyadd(polymul(num1, den2), polymul(num2, den1))
    den = polymul(den1, den2)

    return num, den

# TODO:  Check well for compactibility with FOTF Object
def _convert_to_transfer_function(sys, **kw):
    """Convert a system to transfer function form (if needed).

    If sys is already a transfer function, then it is returned.  If sys is a
    state space object, then it is converted to a transfer function and
    returned.  If sys is a scalar, then the number of inputs and outputs can be
    specified manually, as in:

    >>> sys = _convert_to_transfer_function(3.) # Assumes inputs = outputs = 1
    >>> sys = _convert_to_transfer_function(1., inputs=3, outputs=2)

    In the latter example, sys's matrix transfer function is [[1., 1., 1.]
                                                              [1., 1., 1.]].

    If sys is an array-like type, then it is converted to a constant-gain
    transfer function.

    >>> sys = _convert_to_transfer_function([[1., 0.], [2., 3.]])

    In this example, the numerator matrix will be
       [[[1.0], [0.0]], [[2.0], [3.0]]]
    and the denominator matrix [[[1.0], [1.0]], [[1.0], [1.0]]]

    """
    from statesp import StateSpace

    if isinstance(sys, FOTransFunc):
        if len(kw):
            raise TypeError("If sys is a FOTTransferFunction, " +
                            "_convertTo FOTransFunc cannot take keywords.")

        return sys
    elif isinstance(sys, StateSpace):

        if 0 == sys.states:
            # Slycot doesn't like static SS->TF conversion, so handle
            # it first.  Can't join this with the no-Slycot branch,
            # since that doesn't handle general MIMO systems
            num = [[[sys.D[i, j]] for j in range(sys.inputs)] for i in range(sys.outputs)]
            den = [[[1.] for j in range(sys.inputs)] for i in range(sys.outputs)]
        else:
            try:
                from slycot import tb04ad
                if len(kw):
                    raise TypeError(
                        "If sys is a StateSpace, " +
                        "_convertToTransferFunction cannot take keywords.")

                # Use Slycot to make the transformation
                # Make sure to convert system matrices to numpy arrays
                tfout = tb04ad(sys.states, sys.inputs, sys.outputs, array(sys.A),
                               array(sys.B), array(sys.C), array(sys.D), tol1=0.0)

                # Preallocate outputs.
                num = [[[] for j in range(sys.inputs)] for i in range(sys.outputs)]
                den = [[[] for j in range(sys.inputs)] for i in range(sys.outputs)]

                for i in range(sys.outputs):
                    for j in range(sys.inputs):
                        num[i][j] = list(tfout[6][i, j, :])
                        # Each transfer function matrix row
                        # has a common denominator.
                        den[i][j] = list(tfout[5][i, :])

            except ImportError:
                # If slycot is not available, use signal.lti (SISO only)
                if sys.inputs != 1 or sys.outputs != 1:
                    raise TypeError("No support for MIMO without slycot.")

                # Do the conversion using sp.signal.ss2tf
                # Note that this returns a 2D array for the numerator
                num, den = sp.signal.ss2tf(sys.A, sys.B, sys.C, sys.D)
                num = squeeze(num)  # Convert to 1D array
                den = squeeze(den)  # Probably not needed

        return FOTransFunc(num, den, sys.dt)

    elif isinstance(sys, (int, float, complex, np.number)):
        if "inputs" in kw:
            inputs = kw["inputs"]
        else:
            inputs = 1
        if "outputs" in kw:
            outputs = kw["outputs"]
        else:
            outputs = 1

        num = [[[sys] for j in range(inputs)] for i in range(outputs)]
        den = [[[1] for j in range(inputs)] for i in range(outputs)]

        return FOTransFunc(num, den)

    # If this is array-like, try to create a constant feedthrough
    try:
        D = array(sys)
        outputs, inputs = D.shape
        num = [[[D[i, j]] for j in range(inputs)] for i in range(outputs)]
        den = [[[1] for j in range(inputs)] for i in range(outputs)]
        return FOTransFunc(num, den)
    except Exception as e:
        print("Failure to assume argument is matrix-like in"
              " _convertToTransferFunction, result %s" % e)

    raise TypeError("Can't convert given type to TransferFunction system.")

def fotf(*args):
    """fotf(num, nnum, den, nden[, dt])

    Create a transfer function system. May be able to create a MIMO systems in the future.

    The function accepts either 3, 4 or 5 parameters:

    ``fotf(sys)``
        Convert a linear system into transfer function form. Always creates
        a new system, even if sys is already a TransferFunction object.

    ``fotf(num, nnum, den, nden)``
        Create a transfer function system from its numerator and denominator
        polynomial coefficients.

        If `num` and `den` are 1D array_like objects, the function creates a
        SISO system.

        To create a MIMO system, `num` and `den` need to be 2D nested lists
        of array_like objects. (A 3 dimensional data structure in total.)
        (For details see note below.)

    ``fotf(num, nnum, den, nden, dt)``
        Create a discrete time transfer function system; dt can either be a
        positive number indicating the sampling time or 'True' if no
        specific timebase is given.


    Parameters
    ----------
    sys: LTI (StateSpace or TransferFunction) from Control Toolbox
        A linear system
    num: array_like, or list of list of array_like
        Polynomial coefficients of the numerator
    nden: array_like, or List of list of the coefficients
        of the numerators' power of 's'
    den: array_like, or list of list of array_like
        Polynomial coefficients of the denominator
    nden: array_like, or List of list of the coefficients
        of the denominators' power of 's'

    Returns
    -------
    out: :class:`TransferFunction`
        The new fractional transfer function linear system

    Raises
    ------
    ValueError
        if length of `num` and `den` or if length of nnum != length of num
        or if length of 'den' != length of 'nden' have invalid or unequal dimensions
    TypeError
        if `num` or `den` or 'nnum' or 'nden' are of incorrect type

    See Also
    --------
    TransferFunction
    ss
    ss2tf
    tf2ss

    Notes
    -----
    ``num[i][j]`` contains the polynomial coefficients of the numerator
    for the Fractional order  transfer function from the (j+1)st input to the (i+1)st output.
    ``den[i][j]``, nnum[i][j], nden[i][j] works the same way.


    The list ``[2, 3, 4], [0.5, 0.25, 0]`` denotes the polynomial :math:`2s^0.5 + 3s0.2 + 4`.

    Examples
    --------
    >>> # Create a MIMO transfer function object
    >>> # The transfer function from the 2nd input to the 1st output is
    >>> # (s^1.25 + 2) / (9s^2.5 + 8s^0.25 + 7).
    >>> # (3s^1.25 + 4) / (6.5s^2.5 + 5s^0.25 + 4).
    >>> # (5^1.25 + 6) / (3s^2.5 + 2s^0.25 + 1).
    >>> # (7s^1.25 + 8) / (-s^2.5 + -s^0.25 + 3).
    >>> num = [[[1., 2.], [3., 4.]], [[5., 6.], [7., 8.]]]
    >>> nnum = [[[1.25, 0.], [1.25, 0.]], [[1.25, 0.], [1.25, 0.]]]
    >>> den = [[[9., 8., 7.], [6.5, 5., 4.]], [[3., 2., 1.], [-1., -2., -3.]]]
    >>> nden = [[[2.5, 0.25, 0.], [2.5, 0.25, 0.]], [[2.5, 0.25, 0.], [2.5, 0.25, 0.]]]
    >>> dt=[0]
    >>> sys1 = fotf(num, nnum, den, nden)



    >>> # Convert a StateSpace to a TransferFunction object.
    >>> sys_ss = ss("1. -2; 3. -4", "5.; 7", "6. 8", "9.")
    >>> sys2 = fotf(sys1)

    """

    if len(args) == 1 or len(args) == 3 or len(args) == 4 or len(args) == 5:
        return FOTransFunc(*args)

    # elif len(args) == 1:
    #     from statesp import StateSpace
    #     sys = args[0]
    #     if isinstance(sys, StateSpace):
    #         return ss2tf(sys)
    #     elif isinstance(sys, FOTransFunc):
    #         return deepcopy(sys)
    #     else:
    #         raise TypeError("tf(sys): sys must be a StateSpace or TransferFunction object. "
    #                         "It is %s." % type(sys))
    elif len(args) == 2 and isinstance(args[0], (str,float)) and isinstance(args[1], (str,float)):
        return FOTransFunc(str2poly(args[0], 's'), str2poly(args[1], 's'))
    else:
        raise ValueError("fotf: Needs 1 or 3 or 4 pr 5 arguments; received %i." % len(args))

def newfotf(*args):
    """
    Creates a new FOTF object

    Usgae:
        G=newfotf(SI,S2,t) creates a new FOTF objects from provided
         S1 (zero polynomial) and S2 (pole polynomial). t - optional ioDelay parameter [sec].
         Polynomials can also be independently represented by matched float vectors [num, expnum]
         and [den, expden]
    :param zeroPoly: zero polynomial strings or float vectors in form [num, expnum]
    :param polePoly: zero polynomial strings or float vectors in form [den, expden]
    :param t: - optional ioDelay parameter [sec].
    :return: FOTransFunc
    """
    if len(args) < 2:
        raise ValueError("fotf.newfotf:NotEnoughInputArguments")
    elif len(args) == 2:
        dt = 0
        bases = 's'
    elif len(args) == 3:
        delay = args[2]
        if isinstance(delay,(float,int)) and delay is not 0:
            bases = 'z'
            dt = delay
        else:
            bases = 's'
            dt = 0
    else:
        raise ValueError("fotf.newfotf:TooManyInputArguments")

    #ZerosPoly
    if isinstance(args[0],(list,ndarray)):
        _num = args[0]
        len_num = len(_num)/2
        num = _num[0:len_num]
        nnum = _num[len_num+1:-1]
    elif isinstance(args[0],(int,float)):
        num = args[0]
        nnum = 0
    else:
        num,nnum = str2poly(args[0],bases)

    #PolesPoly
    if isinstance(args[1],(list,ndarray)):
        _den = args[1]
        len_num = len(_den)/2
        den = _num[0:len_num]
        nden = _num[len_num+1:-1]
    elif isinstance(args[1],(int,float)):
        den = args[1]
        nden = 0
    else:
        den,nden = str2poly(args[1], bases)

    return FOTransFunc(num, nnum, den, nden, dt)

#TODO: Sort for FOTransfunc
def ss2tf(*args):
    """ss2tf(sys)

    Transform a state space system to a transfer function.

    The function accepts either 1 or 4 parameters:

    ``ss2tf(sys)``
        Convert a linear system into space system form. Always creates a
        new system, even if sys is already a StateSpace object.

    ``ss2tf(A, B, C, D)``
        Create a state space system from the matrices of its state and
        output equations.

        For details see: :func:`ss`

    Parameters
    ----------
    sys: StateSpace
        A linear system
    A: array_like or string
        System matrix
    B: array_like or string
        Control matrix
    C: array_like or string
        Output matrix
    D: array_like or string
        Feedthrough matrix

    Returns
    -------
    out: TransferFunction
        New linear system in transfer function form

    Raises
    ------
    ValueError
        if matrix sizes are not self-consistent, or if an invalid number of
        arguments is passed in
    TypeError
        if `sys` is not a StateSpace object

    See Also
    --------
    tf
    ss
    tf2ss

    Examples
    --------
    >>> A = [[1., -2], [3, -4]]
    >>> B = [[5.], [7]]
    >>> C = [[6., 8]]
    >>> D = [[9.]]
    >>> sys1 = ss2tf(A, B, C, D)

    >>> sys_ss = ss(A, B, C, D)
    >>> sys2 = ss2tf(sys_ss)

    """

    from .statesp import StateSpace
    if len(args) == 4 or len(args) == 4 or len(args) == 5:
        # Assume we were given the A, B, C, D matrix and (optional) dt
        return _convert_to_transfer_function(StateSpace(*args))

    elif len(args) == 1:
        sys = args[0]
        if isinstance(sys, StateSpace):
            return _convert_to_transfer_function(sys)
        else:
            raise TypeError("ss2tf(sys): sys must be a StateSpace object.  It is %s." % type(sys))
    else:
        raise ValueError("Needs 1 or 4 arguments; received %i." % len(args))


# TODO : sETTLE FOR FOTF
def _clean_part(data):
    """
    Return a valid, cleaned up numerator or denominator
    for the TransferFunction class.

    Parameters
    ----------
    data: numerator or denominator of a transfer function.

    Returns
    -------
    data: list of lists of ndarrays, with int converted to float
    """
    valid_types = (int, float, complex, np.number)
    valid_collection = (list, tuple, ndarray)

    if (isinstance(data, valid_types) or
            (isinstance(data, ndarray) and data.ndim == 0)):
        # Data is a scalar (including 0d ndarray)
        data = [[array([data])]]
    elif (isinstance(data, ndarray) and data.ndim == 3 and
          isinstance(data[0, 0, 0], valid_types)):
        data = [[array(data[i, j])
                 for j in range(data.shape[1])]
                for i in range(data.shape[0])]
    elif (isinstance(data, valid_collection) and
          all([isinstance(d, valid_types) for d in data])):
        data = [[array(data)]]
    elif (isinstance(data, (list, tuple)) and
          isinstance(data[0], (list, tuple)) and
          (isinstance(data[0][0], valid_collection) and
           all([isinstance(d, valid_types) for d in data[0][0]]))):
        data = list(data)
        for j in range(len(data)):
            data[j] = list(data[j])
            for k in range(len(data[j])):
                data[j][k] = array(data[j][k])
    else:
        # If the user passed in anything else, then it's unclear what
        # the meaning is.
        raise TypeError("The numerator and denominator inputs must be scalars or vectors "
                        "(for\nSISO), or lists of lists of vectors (for SISO or MIMO).")

    # Check for coefficients that are ints and convert to floats
    for i in range(len(data)):
        for j in range(len(data[i])):
            for k in range(len(data[i][j])):
                if isinstance(data[i][j][k], (int, np.int)):
                    data[i][j][k] = float(data[i][j][k])

    return data

def fotfparam(fotfobject):
    """
    fotfparam get FOTransFunc object parameters

    :param fotfobject: fotfobject
    :return:  FOTransFunc.num, FOTransFunc.nnum, FOTransFunc.den, FOTransFunc.nden,FOTransFunc.den, FOTransFunc.dt
        =>  den, nden - pole polynomial coefficients and exponents,
            num, nnum - zero polynomial coefficients and exponents,
            dt    - ioDelay [sec]
    """
    if isinstance(fotfobject, FOTransFunc):
        return fotfobject.num[0][0], fotfobject.nnum[0][0], fotfobject.den[0][0], fotfobject.nden[0][0], fotfobject.dt
    else:
        raise ValueError("fotf.fotfparam: Input should be of type 'FOTransFunc' not {}".format(type(fotfobject)))

def fix_s(b):
    """
    Extract integer part from a given number
    :param b: List/numpy Array of real numbers
    :return a: List of integer numbers
    """

    if isinstance(b, list):
        a = []  # output array
        for number in b:
            a = a.append(int(number))  # extrct only the integer part
    elif isinstance(b, np.ndarray):
        a = np.array(b, dtype=np.intc)
    else:
        raise ValueError("fotf.fix_s: input should be a type 'list' on 'numpy.array' with real numbers not a {}".format(type(b)))

    return a


def comm_order(G, type=None):
    comm_factor = MIN_COMM_ORDER ** -1
    n = None
    ord = None
    if isinstance(G, FOTransFunc):
        num, nnum, den, nden, dt = fotfparam(G)
        if type is None:
            a, b = nnum * comm_factor, nden * comm_factor
            newa = np.array(a, dtype=np.int32)
            newb = np.array(b, dtype=np.int32)
            c = np.concatenate((newa, newb), axis=None)
            n = np.gcd.reduce(c)
        else:
            if type is 'num':
                a = nnum * comm_factor
                newa = np.array(a, dtype=np.int32)
            elif type is 'den':
                a = nden * comm_factor
                newa = np.array(a, dtype=np.int32)
            else:
                raise ValueError("fotf.comm_order: Second variable must be either 'num' or 'den'.")

            # Check if array has a single element
            if newa.size == 1:
                if newa[0] < 1:
                    newa[0] = 1
                n = newa[0]
            else:
                # Number of elements greater than 1
                n = np.gcd.reduce(newa)
        ord = n / comm_factor
        return ord
    else:
        raise ValueError("fotf.comm_order: paramter should be of type 'FOTransFunc' not {}".format(type(G)))



    return [rmagDb, rangleCalcDeg, w]#, [ mag, phase, w]
    # inverse fourier transform. optimization(for identification- closed loop) and control.
    # open CUA foR Communication in industries
    # step/impulse response, Simulation, Convolution


def lsim(G, u, t, plot = False):
    """Linear simulation of a fractional-order dynamic system.
    :param G: fractional-order transfer function
    :type G: FOTransFunc
    :param u: input signal vector
    :type u: list
    :param t: time sample vector
    :type t: ndarray or list
    :param plot: if you want a graph
    :type plot: bool

    :returns :y - Time Domain Simulation

    """

    num, nnum, den, nden, dt = fotfparam(G)
    if not isinstance(t, np.ndarray):
        t = np.array(t)

    if not isinstance(u, np.ndarray):
        u = np.array(u,dtype=np.float_)

    sizeden = den.size
    detlaT=t[1]-t[0]
    D= np.sum(den/np.power(detlaT,nden))

    sizeT = t.size
    rnT= range(sizeT)
    vec = np.append(nden,nnum)
    D1 = num/np.power(detlaT,nnum)
    y1 = np.zeros(sizeT)
    W = np.ones((sizeT,vec.size))
    for j in rnT[1:]:
        W[j] = W[j-1]*(1-(vec+1)/j)

    for i in rnT[1:]:
        if i is 1:
            A = y1[i-1::-1] * W[i, 0:sizeden]
        else:
            abb = y1[i-1::-1]
            acc = W[1:i+1,0:sizeden]
            A = abb @ acc
        y1[i] = (u[i] - np.sum((A * den)/np.power(detlaT,nden)))/D

    y = np.zeros(sizeT)
    for i in rnT[1:]:
        bbb = W[0:i+1,sizeden:]
        bcc = bbb @ D1
        bdd = y1[i::-1]
        y[i] = bcc @ bdd

    ysize = y.size
    #Account for I/O delay
    if dt is not None and dt > 0:
        ii= np.where(t>dt)
        ii = np.array(ii, dtype=int)
        # There is a possibility that the value of the sampling interval
        # is greater or equal to the delay. In this case we disregard it.
        if ii is not None:
            lz = np.zeros(ii[0][0])
            ystrip= y[:(ysize - lz.size)]
            y = np.concatenate([lz,ystrip]) # check that y is a column vector also

    if plot:
        plt.figure()
        plt.plot(t, y)
        plt.title('Linear system simulation')
        plt.xlabel('Time [s]')
        plt.ylabel('Amplitude')
        plt.grid()
        plt.show()

    else:
        return y

def step_auto_range(G):
    """
    STEP_AUTO_RANGE Automatically locate appropriate step response length
    """

    #No. of points for auto-ranging (if dcgain exists) and vector generation
    AUTORANGE_NUM_PTS_AUTO = 10
    AUTORANGE_NUM_PTS_GEN  = 10000

    #Detection coefficient: error band
    AUTORANGE_DET_COEF = 0.2

    #Add this value to the exponent in case of oscillating processes
    AUTORANGE_OSC_ADD  = 1

    #Default end point in case dcgain gives "0" or "Inf"
    AUTORANGE_DEFAULT_END = 100

    #Test DC gain
    myGain = dcgain(G)

    if np.isinf(myGain) or (np.abs(myGain) < np.finfo(float).resolution):
        t  = np.linspace(0,AUTORANGE_DEFAULT_END,AUTORANGE_NUM_PTS_GEN)
    else:
        #Final time range locator
        t_exp = -4
        t_exp_max = 9

        while t_exp<t_exp_max:
            t = linspace(0, 10**t_exp, AUTORANGE_NUM_PTS_AUTO)
            u = np.ones(len(t))
            y = lsim(G,u,t,plot=False)
            if (abs(y[-1]/np.abs(myGain))>=AUTORANGE_DET_COEF):
                # Check if this is an oscillating process, and add a few
                # exponent values to ensure that it is covered in full
                if max(np.abs(y)) > np.abs(y[-1]):
                    t_exp += AUTORANGE_OSC_ADD
                break
            t_exp += 1
        # Use the obtained value to generate the vector
        t = linspace(0,10**t_exp, AUTORANGE_NUM_PTS_GEN)
    return t



def dcgain(G):
    num, nnum, den, nden, dt = fotfparam(G)
    #evaluate value at s = 0. from observation only the terms that have an exponent of s == 0
    dcnum = np.sum(num * (nnum == 0),dtype=float)
    dcden = np.sum(den * (nden == 0),dtype=float)

    K=dcnum/dcden
    return K

def tfdata(sys):
    """
    Returns numerator and denominator of a fractional transfer function in a
    rational transfer function form with commensurate order q

    Parameters
    ----------
    sys: FOTransFunc (Fractional order TransferFunction)

    Returns
    -------
    (num, den, q): numerator and denominator arrays and  commesurate order
    """
    #tf = _convert_to_transfer_function(sys)
    if isinstance(sys,FOTransFunc):
        #Get commensurate order and max order
        q = comm_order(sys)
        num,nnum,den,nden,dt = fotfparam(sys)

        #Get elements positions
        a1=fix_s(nden/q)
        b1=fix_s(nnum/q)

        #Numerator
        newnum = np.zeros(b1[0]+1,dtype=np.float_)
        for i in range(num.size):
            ii=b1[i]
            newnum[ii] = num[i]
        newnum = np.flip(newnum)

        # Denumerator
        newden = np.zeros(a1[0]+1,dtype=np.float_)
        for j in range(den.size):
            jj=a1[j]
            newden[jj] = den[j]
        newden = np.flip(newden)
    return newnum, newden, q

def simple(G):
    _num,_nnum,_den,_nden, _dt = fotfparam(G)
    num, nnum = polyuniq(_num,_nnum)
    den, nden = polyuniq(_den,_nden)

    nn = min([nnum.min(), nden.min()])
    nnum -= nn
    nden -= nn

    return FOTransFunc(num,nnum,den,nden,_dt)


def polyuniq(num, nnum):
    if not isinstance(num,ndarray):
        _num = np.array(num)
    else:
        _num = deepcopy(num)

    if not isinstance(nnum,ndarray):
        _nnum = np.array(nnum)
    else:
        _nnum = deepcopy(nnum)

    ind = np.lexsort((_num,_nnum))
    if ind.size != 1:
        for i in range(len(ind)):
            indx = ind[i]
            num[i] = _num[indx]
            nnum[i] = _nnum[indx]
        return num[::-1], nnum[::-1]
    else:
        return num, nnum


def impulse(G,tt = None, plot = True):

    if tt is None:
        t = step_auto_range(G)
    elif not isinstance(tt,ndarray):
        t = np.array(tt)
    else:
        t=tt

    #An approximation of impulse response
    y = step(G * fotf('s'), t, output=True)

    if plot:
        plt.figure()
        plt.plot(t,y)
        plt.title('Impulse Response')
        plt.xlabel('Time [s]')
        plt.ylabel('Amplitude')
        plt.grid()
        plt.show()

        #Plot Final value if present
        #Test DCGain

        myGain = dcgain(G)
        if np.isinf(myGain) or (np.abs(myGain) < np.finfo(float).resolution):
            pass
        else:
            plt.figure()
            plt.plot([t[0], t[-1]], [myGain, myGain], ':k')
            plt.title('Dc Gain')
            plt.xlabel('Time [s]')
            plt.ylabel('Amplitude')
            plt.grid()
            plt.show()
    return t,y

def trunc(G,numAcc, nnumAcc):
    """
    Truncate the exponents and coefficients of a fractional-order transfer function.
    :param G: Initial System
    :type G:    FOTransFunc
    :param numAcc: Coefficient Accuracy
    :param nnumAcc: Exponent Accuracy e.g 1e-2 means 2 decimal places
    :return: FOTransfunc
    """

    num,nnum,den,nden,dt = fotfparam(G)
    if numAcc != 0:
        num = fix_s(num/numAcc) * numAcc
        den = fix_s(den/numAcc) * numAcc

    if nnumAcc != 0:
        nnum = fix_s(nnum / numAcc) * nnumAcc
        nden = fix_s(nden / numAcc) * nnumAcc

    return simple(fotf(num,nnum,den,nden,dt))

