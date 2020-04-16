'''
Class for creating FOPID controllers
'''

import numpy as np
import math
import warnings
import sys
import time
from unsync import unsync

# Precision
EPSILON = sys.float_info.epsilon

# This implements a discrete oustaloup approximation of a FOPID controller
# in the form of second-order sections of an IIR filter
class fpid_2iir:
    # FOPID parameters
    Kp = Ki = Kd = lam = mu = 0

    # Approximation parameters
    a_wb = 0.0001
    a_wh = 10000
    a_N = 5
    a_Ts = 0.01

    # Flags. Since Python is not a real-time language in any way,
    # this may be overkill. But as always, safety first!
    COMPUTING_FOPID_APPROXIMATION = False
    COMPUTING_CONTROL_LAW = False

    # Storage for all necessary parameters

    # Second order filters
    KIc = None  # Fractional I component IIR filter gain
    I_zsos = None  # Approximation zeros and poles
    I_psos = None
    KDc = None  # Same for D component
    D_zsos = None
    D_psos = None

    # Last call timestamp
    last_timestamp = None
    first_control = True

    # Allowable clock jitter in %
    ALLOWABLE_CLOCK_JITTER = 5

    # Running state information
    s_I = None
    s_D = None
    s_IntMem = None  # Conventional integrator memory
    in_Mem = None  # Previous state

    # Initialize the object
    def __init__(self, params):
        if type(params) is not dict:
            raise Exception("params must be a dict with info about the controller and optionally the approximation")

        # Get the params out
        fpid_params = params['fpid']  # Parameters of the controller

        # Controller
        self.Kp = fpid_params['Kp']
        self.Ki = fpid_params['Ki']
        self.Kd = fpid_params['Kd']
        self.lam = fpid_params['lam']
        self.mu = fpid_params['mu']

        # If not, use default parameters
        if "oust" in params:
            oust_params = params['oust']  # Parameters of approximation

            # Approximation
            self.a_wb = oust_params['wb']
            self.a_wh = oust_params['wh']
            self.a_N = oust_params['N']
            self.a_Ts = oust_params['Ts']

        self.compute_fopid_approximation()

    # Create the approximation
    def compute_fopid_approximation(self):
        self.COMPUTING_FOPID_APPROXIMATION = True

        # Fractional integral and differential component approximations
        self.KIc, self.I_zsos, self.I_psos = self.compute_iir_sos_oustaloup(1 - self.lam)
        self.KDc, self.D_zsos, self.D_psos = self.compute_iir_sos_oustaloup(self.mu)

        # Clear IIR memory
        self.clear_iir_memory()
        self.COMPUTING_FOPID_APPROXIMATION = False

    # Oustaloup approximation
    def compute_iir_sos_oustaloup(self, alpha):
        T = self.a_Ts
        wb = self.a_wb
        wh = self.a_wh
        N = self.a_N

        # Set correct upper frequency bound
        wh = 2 / T if wh > (2 / T) else wh
        omu = wh / wb

        # First we compute the poles and zeros of the approximation
        # and map them onto discrete time
        zz = np.zeros((1, 2 * N + 1))
        zp = np.zeros((1, 2 * N + 1))

        for k in range(-N, N + 1):
            w_kp = wb * math.pow(omu, ((k + N + 0.5 - 0.5 * alpha) / (2 * N + 1)))
            w_k = wb * math.pow(omu, ((k + N + 0.5 + 0.5 * alpha) / (2 * N + 1)))
            zz[0, k + N] = math.exp(-T * w_kp)
            zp[0, k + N] = math.exp(-T * w_k)

        wu = math.sqrt(wb * wh)
        Ks = math.pow(wu, alpha)

        theta = math.cos(wu * T)
        Ku = 1

        # Need to compute correction gain
        for k in range(0, 2 * N + 1):
            nk = (1 - 2 * zz[0, k] * theta + zz[0, k] * zz[0, k])
            dk = (1 - 2 * zp[0, k] * theta + zp[0, k] * zp[0, k])

            if nk > EPSILON and dk > EPSILON:
                Ku *= (nk / dk)

        Ku = math.sqrt(Ku)
        Ks = Ks / Ku

        # Second order section form
        zca = np.zeros((N + 1, 2))
        pca = np.zeros((N + 1, 2))

        # First order section
        zca[0, 0] = -zz[0, 2 * N]
        zca[0, 1] = 0
        pca[0, 0] = -zp[0, 2 * N]
        pca[0, 1] = 0
        for k in range(N - 1, -1, -1):
            zca[N - k, 0] = -zz[0, 2 * k + 1] - zz[0, 2 * k]
            zca[N - k, 1] = zz[0, 2 * k + 1] * zz[0, 2 * k]
            pca[N - k, 0] = -zp[0, 2 * k + 1] - zp[0, 2 * k]
            pca[N - k, 1] = zp[0, 2 * k + 1] * zp[0, 2 * k]

        # Return tuple
        return (Ks, zca, pca)

    # Clear IIR memory
    def clear_iir_memory(self):
        self.s_I = np.zeros((self.a_N + 1, 2))
        self.s_D = np.zeros((self.a_N + 1, 2))
        self.s_IntMem = 0

    # Do IIR filtering

    def do_iir_filtering(self, z, p, kc, mem, inp):
        # Local names for arrays for convenience
        b = z
        a = p
        s = mem
        u_n = kc * inp
        for m in range(0, self.a_N + 1):
            y_n = u_n + s[m, 0]
            s[m, 0] = b[m, 0] * u_n - a[m, 0] * y_n + s[m, 1]
            s[m, 1] = b[m, 1] * u_n - a[m, 1] * y_n
            u_n = y_n

        return u_n

    # Compute the actual control law
    def compute_fopid_control_law(self, err):

        self.COMPUTING_CONTROL_LAW = True
        t_old = self.last_timestamp
        self.last_timestamp = time.time()

        if self.COMPUTING_FOPID_APPROXIMATION:
            warnings.warn("Cannot compute control law while the approximation is being computed.")
            return

        # Because Python does not support hard real-time operation,
        # we need to check whether clock jitter does not factor into control
        # if not self.first_control:
        #     t_diff = self.last_timestamp - t_old
        #     clockJitter = (t_diff / self.a_Ts)*100
        #     if clockJitter > self.ALLOWABLE_CLOCK_JITTER:    #Changed '*100' to '/self.a_Ts'
        #         print("Maximum allowable clock jitter: {0}.\nCurrent Clock Jitter: {1} exceeded.".format(self.ALLOWABLE_CLOCK_JITTER,clockJitter))

        self.first_control = False
        inp = err
        self.in_Mem = err

        foi_out = self.Ki * self.do_iir_filtering(self.I_zsos, self.I_psos, self.KIc, self.s_I, inp)
        fod_out = self.Kd * self.do_iir_filtering(self.D_zsos, self.D_psos, self.KDc, self.s_D, inp)

        self.s_IntMem += self.a_Ts * foi_out
        i_out = self.s_IntMem
        out = self.Kp * inp + i_out + fod_out
        self.COMPUTING_CONTROL_LAW = False

        return out