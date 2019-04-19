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
from itertools import chain
import re
from lti import LTI

__all__ = ['FOTransFunc', 'fotf', 'ss2tf', 'fotfdata', 'poly2str', 'str2poly']


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
            if len(_args) == 5:
                _dt = _args[4]
        elif len(_args) == 1 and type(_args[0]) == float:
            _num = 1.
            _nnum = 0
            _den = _args[0]
            _nden = 0
            _dt = 0

        elif len(_args) == 1 and type(_args[0]) == 's':
            _num = 1.
            _nnum = 0
            _den = 1.
            _nden = 1.
            _dt = 0

        elif len(_args) == 1 and isinstance((_args[0]), FOTransFunc):
            [_num, _nnum, _den, _nden]= _args[0:4]
            if len(_args)==5:
                _dt = _args[4]

        elif len(_args) >= 2 and isinstance(_args[0], str) and isinstance(_args[1], str):
            _num = _args[0]
            _nnum = _args[1]

            if len(_args) >= 3 and (isinstance(_args[2], float) or isinstance(_args[2], int)):
                _dt = abs(_args[2])

        else:
            raise ValueError("Needs 1, 2 , 3 ,4 or 5 arguments; received {}.".format(len(args)))

        _num = _clean_part(_num)
        _den = _clean_part(_den)
        _nnum = _clean_part(_nnum)
        _nden = _clean_part(_nden)

        inputs = len(_num[0])
        outputs = len(_num)
        npower = len(_nnum[0])

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

        super().__init__(self, inputs, outputs, _dt)
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
        """Add two LTI objects (parallel connection)."""
        from .statesp import StateSpace

        # Convert the second argument to a transfer function.
        if isinstance(other, StateSpace):
            other = _convert_to_transfer_function(other)
        elif not isinstance(other, FOTransFunc):
            other = _convert_to_transfer_function(other, inputs=self.inputs,
                                                  outputs=self.outputs)

        # Check that the input-output sizes are consistent.
        if self.inputs != other.inputs:
            raise ValueError("The first summand has %i input(s), but the second has %i."
                             % (self.inputs, other.inputs))
        if self.outputs != other.outputs:
            raise ValueError("The first summand has %i output(s), but the second has %i."
                             % (self.outputs, other.outputs))

        # Figure out the sampling time to use
        if self.dt is None and other.dt is not None:
            dt = other.dt  # use dt from second argument
        elif (other.dt is None and self.dt is not None) or (timebaseEqual(self, other)):
            dt = self.dt  # use dt from first argument
        else:
            raise ValueError("Systems have different sampling times")

        # Preallocate the numerator and denominator of the sum.
        num = [[[] for j in range(self.inputs)] for i in range(self.outputs)]
        den = [[[] for j in range(self.inputs)] for i in range(self.outputs)]

        for i in range(self.outputs):
            for j in range(self.inputs):
                num[i][j], den[i][j] = _add_siso(self.num[i][j], self.den[i][j],
                                                 other.num[i][j],
                                                 other.den[i][j])

        return FOTransFunc(num, den, dt)

    def __radd__(self, other):
        """Right add two LTI objects (parallel connection)."""
        return self + other

    def __sub__(self, other):
        """Subtract two LTI objects."""
        return self + (-other)

    def __rsub__(self, other):
        """Right subtract two LTI objects."""
        return other + (-self)

    # TODO: investigate on how to do for siso and mimo, will be needed for feedback
    def __mul__(self, other):
        """Multiply two LTI objects (serial connection)."""
        # Convert the second argument to a transfer function.
        if isinstance(other, (int, float, complex, np.number)):
            other = _convert_to_transfer_function(other, inputs=self.inputs,
                                                  outputs=self.inputs)
        else:
            other = _convert_to_transfer_function(other)

        # Check that the input-output sizes are consistent.
        if self.inputs != other.outputs:
            raise ValueError("C = A * B: A has %i column(s) (input(s)), but B has %i "
                             "row(s)\n(output(s))." % (self.inputs, other.outputs))

        inputs = other.inputs
        outputs = self.outputs

        # Figure out the sampling time to use
        if self.dt is None and other.dt is not None:
            dt = other.dt  # use dt from second argument
        elif (other.dt is None and self.dt is not None) or (self.dt == other.dt):
            dt = self.dt  # use dt from first argument
        else:
            raise ValueError("Systems have different sampling times")

        # Preallocate the numerator and denominator of the sum.
        num = [[[0] for j in range(inputs)] for i in range(outputs)]
        den = [[[1] for j in range(inputs)] for i in range(outputs)]

        # Temporary storage for the summands needed to find the (i, j)th element of the product.
        num_summand = [[] for k in range(self.inputs)]
        den_summand = [[] for k in range(self.inputs)]

        # Multiply & add.
        for row in range(outputs):
            for col in range(inputs):
                for k in range(self.inputs):
                    num_summand[k] = polymul(self.num[row][k], other.num[k][col])
                    den_summand[k] = polymul(self.den[row][k], other.den[k][col])
                    num[row][col], den[row][col] = _add_siso(
                        num[row][col], den[row][col],
                        num_summand[k], den_summand[k])

        return FOTransFunc(num, den, dt)

    # TODO: investigate on how to do for siso and mimo, will be needed for feedback
    def __rmul__(self, other):
        """Right multiply two LTI objects (serial connection)."""

        # Convert the second argument to a transfer function.
        if isinstance(other, (int, float, complex, np.number)):
            other = _convert_to_transfer_function(other, inputs=self.inputs,
                                                  outputs=self.inputs)
        else:
            other = _convert_to_transfer_function(other)

        # Check that the input-output sizes are consistent.
        if other.inputs != self.outputs:
            raise ValueError("C = A * B: A has %i column(s) (input(s)), but B has %i "
                             "row(s)\n(output(s))." % (other.inputs, self.outputs))

        inputs = self.inputs
        outputs = other.outputs

        # Figure out the sampling time to use
        if self.dt is None and other.dt is not None:
            dt = other.dt  # use dt from second argument
        elif (other.dt is None and self.dt is not None) \
                or (self.dt == other.dt):
            dt = self.dt  # use dt from first argument
        else:
            raise ValueError("Systems have different sampling times")

        # Preallocate the numerator and denominator of the sum.
        num = [[[0] for j in range(inputs)] for i in range(outputs)]
        den = [[[1] for j in range(inputs)] for i in range(outputs)]

        # Temporary storage for the summands needed to find the
        # (i, j)th element
        # of the product.
        num_summand = [[] for k in range(other.inputs)]
        den_summand = [[] for k in range(other.inputs)]

        for i in range(outputs):  # Iterate through rows of product.
            for j in range(inputs):  # Iterate through columns of product.
                for k in range(other.inputs):  # Multiply & add.
                    num_summand[k] = polymul(other.num[i][k], self.num[k][j])
                    den_summand[k] = polymul(other.den[i][k], self.den[k][j])
                    num[i][j], den[i][j] = _add_siso(
                        num[i][j], den[i][j],
                        num_summand[k], den_summand[k])

        return FOTransFunc(num, den, dt)

    # TODO: Division of MIMO transfer function objects is not written yet.
    def __truediv__(self, other):
        """Divide two LTI objects."""

        if isinstance(other, (int, float, complex, np.number)):
            other = _convert_to_transfer_function(
                other, inputs=self.inputs,
                outputs=self.inputs)
        else:
            other = _convert_to_transfer_function(other)

        if (self.inputs > 1 or self.outputs > 1 or
                other.inputs > 1 or other.outputs > 1):
            raise NotImplementedError(
                "TransferFunction.__truediv__ is currently \
                implemented only for SISO systems.")

        # Figure out the sampling time to use
        if self.dt is None and other.dt is not None:
            dt = other.dt  # use dt from second argument
        elif (other.dt is None and self.dt is not None) or (self.dt == other.dt):
            dt = self.dt  # use dt from first argument
        else:
            raise ValueError("Systems have different sampling times")

        num = polymul(self.num[0][0], other.den[0][0])
        den = polymul(self.den[0][0], other.num[0][0])

        return FOTransFunc(num, den, dt)

    # TODO: Remove when transition to python3 complete
    def __div__(self, other):
        return FOTransFunc.__truediv__(self, other)

    # TODO: Division of MIMO transfer function objects is not written yet.
    def __rtruediv__(self, other):
        """Right divide two LTI objects."""
        if isinstance(other, (int, float, complex, np.number)):
            other = _convert_to_transfer_function(
                other, inputs=self.inputs,
                outputs=self.inputs)
        else:
            other = _convert_to_transfer_function(other)

        if (self.inputs > 1 or self.outputs > 1 or
                other.inputs > 1 or other.outputs > 1):
            raise NotImplementedError(
                "TransferFunction.__rtruediv__ is currently implemented only for SISO systems.")

        return other / self

    # TODO: Remove when transition to python3 complete
    def __rdiv__(self, other):
        return FOTransFunc.__rtruediv__(self, other)

    def __pow__(self, other):
        if not type(other) == int:
            raise ValueError("Exponent must be an integer")
        if other == 0:
            return FOTransFunc([1], [1])  # unity
        if other > 0:
            return self * (self ** (other - 1))
        if other < 0:
            return (FOTransFunc([1], [1]) / self) * (self ** (other + 1))

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
    def freqresp(self, omega):
        """Evaluate a transfer function at a list of angular frequencies.

        mag, phase, omega = self.freqresp(omega)

        reports the value of the magnitude, phase, and angular frequency of the
        transfer function matrix evaluated at s = i * omega, where omega is a
        list of angular frequencies, and is a sorted
        version of the input omega.
        """

        # Preallocate outputs.
        numfreq = len(omega)
        mag = empty((self.outputs, self.inputs, numfreq))
        phase = empty((self.outputs, self.inputs, numfreq))

        # Figure out the frequencies
        omega.sort()
        if isdtime(self, strict=True):
            dt = timebase(self)
            slist = np.array([exp(1.j * w * dt) for w in omega])
            if max(omega) * dt > pi:
                warn("freqresp: frequency evaluation above Nyquist frequency")
        else:
            slist = np.array([1j * w for w in omega])

        # Compute frequency response for each input/output pair
        for i in range(self.outputs):
            for j in range(self.inputs):
                fresp = (polyval(self.num[i][j], slist) /
                         polyval(self.den[i][j], slist))
                mag[i, j, :] = abs(fresp)
                phase[i, j, :] = angle(fresp)

        return mag, phase, omega

    def pole(self):
        """Compute the poles of a transfer function."""
        num, den, denorder = self._common_den()
        rts = []
        for d, o in zip(den, denorder):
            rts.extend(roots(d[:o + 1]))
        return np.array(rts)

    def zero(self):
        """Compute the zeros of a transfer function."""
        if self.inputs > 1 or self.outputs > 1:
            raise NotImplementedError("TransferFunction.zero is currently only implemented "
                                      "for SISO systems.")
        else:
            # for now, just give zeros of a SISO tf
            return roots(self.num[0][0])

    def feedback(self, other=1, sign=-1):
        """Feedback interconnection between two LTI objects."""
        other = _convert_to_transfer_function(other)

        if (self.inputs > 1 or self.outputs > 1 or
                other.inputs > 1 or other.outputs > 1):
            # TODO: MIMO feedback
            raise NotImplementedError("TransferFunction.feedback is currently only implemented "
                                      "for SISO functions.")

        # Figure out the sampling time to use
        if self.dt is None and other.dt is not None:
            dt = other.dt  # use dt from second argument
        elif (other.dt is None and self.dt is not None) or (self.dt == other.dt):
            dt = self.dt  # use dt from first argument
        else:
            raise ValueError("Systems have different sampling times")

        num1 = self.num[0][0]
        den1 = self.den[0][0]
        num2 = other.num[0][0]
        den2 = other.den[0][0]

        num = polymul(num1, den2)
        den = polyadd(polymul(den2, den1), -sign * polymul(num2, num1))

        return FOTransFunc(num, den, dt)

        # For MIMO or SISO systems, the analytic expression is
        #     self / (1 - sign * other * self)
        # But this does not work correctly because the state size will be too
        # large.

    def minreal(self, tol=None):
        """Remove cancelling pole/zero pairs from a transfer function"""
        # based on octave minreal

        # default accuracy
        from sys import float_info
        sqrt_eps = sqrt(float_info.epsilon)

        # pre-allocate arrays
        num = [[[] for j in range(self.inputs)] for i in range(self.outputs)]
        den = [[[] for j in range(self.inputs)] for i in range(self.outputs)]

        for i in range(self.outputs):
            for j in range(self.inputs):

                # split up in zeros, poles and gain
                newzeros = []
                zeros = roots(self.num[i][j])
                poles = roots(self.den[i][j])
                gain = self.num[i][j][0] / self.den[i][j][0]

                # check all zeros
                for z in zeros:
                    t = tol or \
                        1000 * max(float_info.epsilon, abs(z) * sqrt_eps)
                    idx = where(abs(z - poles) < t)[0]
                    if len(idx):
                        # cancel this zero against one of the poles
                        poles = delete(poles, idx[0])
                    else:
                        # keep this zero
                        newzeros.append(z)

                # poly([]) returns a scalar, but we always want a 1d array
                num[i][j] = np.atleast_1d(gain * real(poly(newzeros)))
                den[i][j] = np.atleast_1d(real(poly(poles)))

        # end result
        return FOTransFunc(num, den, self.dt)

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

    def _common_den(self, imag_tol=None):
        """
        Compute MIMO common denominators; return them and adjusted numerators.

        This function computes the denominators per input containing all
        the poles of sys.den, and reports it as the array den.  The
        output numerator array num is modified to use the common
        denominator for this input/column; the coefficient arrays are also
        padded with zeros to be the same size for all num/den.
        num is an sys.outputs by sys.inputs
        by len(d) array.

        Parameters
        ----------
        imag_tol: float
            Threshold for the imaginary part of a root to use in detecting
            complex poles

        Returns
        -------
        num: array
            Multi-dimensional array of numerator coefficients. num[i][j]
            gives the numerator coefficient array for the ith input and jth
            output, also prepared for use in td04ad; matches the denorder
            order; highest coefficient starts on the left.

        den: array
            Multi-dimensional array of coefficients for common denominator
            polynomial, one row per input. The array is prepared for use in
            slycot td04ad, the first element is the highest-order polynomial
            coefficiend of s, matching the order in denorder, if denorder <
            number of columns in den, the den is padded with zeros

        denorder: array of int, orders of den, one per input



        Examples
        --------
        >>> num, den, denorder = sys._common_den()

        """

        # Machine precision for floats.
        eps = finfo(float).eps

        # Decide on the tolerance to use in deciding of a pole is complex
        if (imag_tol is None):
            imag_tol = 1e-8  # TODO: figure out the right number to use

        # A list to keep track of cumulative poles found as we scan
        # self.den[..][..]
        poles = [[] for j in range(self.inputs)]

        # RvP, new implementation 180526, issue #194

        # pre-calculate the poles for all num, den
        # has zeros, poles, gain, list for pole indices not in den,
        # number of poles known at the time analyzed

        # do not calculate minreal. Rory's hint .minreal()
        poleset = []
        for i in range(self.outputs):
            poleset.append([])
            for j in range(self.inputs):
                if abs(self.num[i][j]).max() <= eps:
                    poleset[-1].append([array([], dtype=float),
                                        roots(self.den[i][j]), 0.0, [], 0])
                else:
                    z, p, k = tf2zpk(self.num[i][j], self.den[i][j])
                    poleset[-1].append([z, p, k, [], 0])

        # collect all individual poles
        epsnm = eps * self.inputs * self.outputs
        for j in range(self.inputs):
            for i in range(self.outputs):
                currentpoles = poleset[i][j][1]
                nothave = ones(currentpoles.shape, dtype=bool)
                for ip, p in enumerate(poles[j]):
                    idx, = nonzero(
                        (abs(currentpoles - p) < epsnm) * nothave)
                    if len(idx):
                        nothave[idx[0]] = False
                    else:
                        # remember id of pole not in tf
                        poleset[i][j][3].append(ip)
                for h, c in zip(nothave, currentpoles):
                    if h:
                        poles[j].append(c)
                # remember how many poles now known
                poleset[i][j][4] = len(poles[j])

        # figure out maximum number of poles, for sizing the den
        npmax = max([len(p) for p in poles])
        den = zeros((self.inputs, npmax + 1), dtype=float)
        num = zeros((max(1, self.outputs, self.inputs),
                     max(1, self.outputs, self.inputs), npmax + 1), dtype=float)
        denorder = zeros((self.inputs,), dtype=int)

        for j in range(self.inputs):
            if not len(poles[j]):
                # no poles matching this input; only one or more gains
                den[j, 0] = 1.0
                for i in range(self.outputs):
                    num[i, j, 0] = poleset[i][j][2]
            else:
                # create the denominator matching this input
                # polyfromroots gives coeffs in opposite order from what we use
                # coefficients should be padded on right, ending at np
                np = len(poles[j])
                den[j, np::-1] = polyfromroots(poles[j]).real
                denorder[j] = np

                # now create the numerator, also padded on the right
                for i in range(self.outputs):
                    # start with the current set of zeros for this output
                    nwzeros = list(poleset[i][j][0])
                    # add all poles not found in the original denominator,
                    # and the ones later added from other denominators
                    for ip in chain(poleset[i][j][3],
                                    range(poleset[i][j][4], np)):
                        nwzeros.append(poles[j][ip])

                    numpoly = poleset[i][j][2] * polyfromroots(nwzeros).real
                    # print(numpoly, den[j])
                    # polyfromroots gives coeffs in opposite order => invert
                    # numerator polynomial should be padded on left and right
                    #   ending at np to line up with what td04ad expects...
                    num[i, j, np + 1 - len(numpoly):np + 1] = numpoly[::-1]
                    # print(num[i, j])

        if (abs(den.imag) > epsnm).any():
            print("Warning: The denominator has a nontrivial imaginary part: %f"
                  % abs(den.imag).max())
        den = den.real

        return num, den, denorder

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


# Utility function to convert a transfer function polynomial to a string
# Borrowed from poly1d library


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

@FOTransFunc
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
        #return FOTransFunc(*args)
        return args
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
    else:
        raise ValueError("fotf: Needs 1 or 3 or 4 pr 5 arguments; received %i." % len(args))


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


def fotfdata(sys):
    """
    Return transfer function data objects for a system

    Parameters
    ----------
    sys: LTI (StateSpace, or TransferFunction)
        LTI system whose data will be returned

    Returns
    -------
    (num, den): numerator and denominator arrays
        Transfer function coefficients (SISO only)
    """
    tf = _convert_to_transfer_function(sys)

    return tf.num, tf.den


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


def str2poly(polystr, bases='s'):
    """Converts a a sting representation of a Fractional order transfer function to Polynomial
    represented by a list/array

    Args: Fractional order string, Bases i,e 's' or 'z'.
        default walue is 's'
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

        column = len(polystr[0])
        outstring = []

        for i in range(row):
            polyb = []
            for j in range(column):
                polyb.append(polystr[j][i])
            outstring.append(polyb)

        return outstring
    else:
        raise ValueError("Input should be of format 'str' with bases type 'z' or 's'")

def fotfparam(fotfobject):
    """
    fotfparam get FOTransFunc object parameters

    :param fotfobject: [[1], [den,nden]] or [fotfobject]
    :return:
        [FOTransFunc.den, FOTransFunc.nden] => pole polynomial coefficients and exponents
        [FOTransFunc.num, FOTransFunc.nnum, FOTransFunc.den, FOTransFunc.nden,FOTransFunc.den, FOTransFunc.dt]
        =>  A, NA - pole polynomial coefficients and exponents,
            B, NB - zero polynomial coefficients and exponents,
            T     - ioDelay [sec]
    """
    if isinstance(fotfobject, FOTransFunc) and len(fotfobject) == 2:
        return [fotfobject.den, fotfobject.nden]
    else:
        return [fotfobject.num, fotfobject.nnum, fotfobject.den, fotfobject.nden, fotfobject.dt]
