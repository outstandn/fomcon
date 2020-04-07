import numpy as np
from addict import Dict
import math

# ****************************
# Custom structure definitions using dict()
# ****************************

# IO -----------------------------------------------
T_reading = dict(T_tc = '', T_internal = '', state = '')
# --------------------------------------------------

# Oustaloup's approximation parameters
oustapp_params = dict(wb = 0, wh = 0, N = 0, Ts = 0)

# FO PID controller parameters
fopid = dict(Kp = 1, Ki = 0, Kd = 0, lamda = 1, mu = 1)

# FO FOPDT model parameters
fofopdt_model = dict(K=0, L =0, T =0, alpha =0)

# FO control system tuning parameters
design_specs = dict(wc=0, pm =0, opt_norm =0)

# Solution to a system of linear equations (3)
sle_sol = dict(x1=0, x2=0, x3=0)

#**************
# Constants
# **************
NMAX = 10			    # Max order of Oustaloup approximation
NSEC = NMAX+1
PZMAX = (NMAX*2)+1      # Max number of poles (zeros)
EPS	= 1.0e-10           #Epsilon
NUMMAX = 1.0e10
M_PI  = math.pi           #3.14159265358979323846

# **************
# I/O saturation
# **************
# Input scaling
MAX_SCALED_IN = 1.0
MIN_SCALED_IN = -1.0
INPUT_MARGIN = 0.6      # Input change rate limit before sample memory is cleared

#Output scaling
MAX_SCALED_OUT = 1.0
MIN_SCALED_OUT = -1.0

#Offset
MED_SCALED_IN = 0.5 * (MAX_SCALED_IN + MIN_SCALED_IN)
MED_SCALED_OUT = 0.5 * (MAX_SCALED_OUT + MIN_SCALED_OUT)

#Factor
SCALE_FACTOR_IN = (MAX_SCALED_IN  - MED_SCALED_IN)
SCALE_FACTOR_OUT = (MAX_SCALED_OUT - MED_SCALED_OUT)

#******* Required sample time *******
DT = 0.01  		         #[s]
SAMPLE_RATE	= 1 / DT
#************************************

#************************************
# GLOBAL variables
#************************************

# Coefficient generation temporary static storage
zz = np.zeros((PZMAX,1), dtype=float, order='C')		# Discrete-time zero array
zp= np.zeros((PZMAX,1), dtype=float, order='C')        #Discrete-time pole array

# IIR SO section coefficient storage for I and D components
KIc = 0.0	#IIR filter gain
I_zsos = np.zeros((NSEC,2), dtype=float, order='C')	    #Zero polynomial second-order sections
I_psos = np.zeros((NSEC,2), dtype=float, order='C')	#Pole polynomial second-order sections
KDc = 0.0
D_zsos = np.zeros((NSEC,2), dtype=float, order='C')
D_psos = np.zeros((NSEC,2), dtype=float, order='C')

# States
s_I = np.zeros((NSEC,2), dtype=float, order='C')
s_D = np.zeros((NSEC,2), dtype=float, order='C')
s_IntMem = np.zeros((NSEC,2), dtype=float, order='C')		 # Regular integrator memory
in_Mem  = np.zeros((NSEC,2), dtype=float, order='C')		 # Previous sample

# FOPID controller parameters: defaults to regular PID controller with direct feed-through
the_fopid = Dict(fopid)
# the_fopid.Kp = 1
# the_fopid.Ki = 0
# the_fopid.Kd = 0
# the_fopid.lamda = 1
# the_fopid.mu = 1

# The FOPID for which parameters are sought
und_fopid = Dict(fopid)
# und_fopid.Kp = 1
# und_fopid.Ki = 0
# und_fopid.Kd = 0
# und_fopid.lamda = 1
# und_fopid.mu = 1

# Search for

# Computation flags
flag_FOPID_Ready = False				# FOPID coefficients are being computed
flag_FOPID_Computing_Output = False		# FOPID is still computing the output sample
flag_FOPID_Schedule_Generation = False		# FOPID generation has been scheduled and will take place once the flag_FOPID_Computing_Output will be cleared.

params = Dict(oustapp_params)							# Approximation parameters
# params.wb = 0
# params.wh = 0
# params.N = 0
# params.Ts = 0

the_fofopdt = Dict(fofopdt_model)						# The model
dspecs = Dict(design_specs)							# Design specifications

# Optimization norm
OPT_NORM = 1e-4
OPT_MAX_ITER = 10

# Jacobian and specification function vector are shared
Jac = np.zeros((3,3), dtype =float , order = 'C')
kappa_vec = np.zeros((3,1), dtype =float , order = 'C')

# DEBUG: iterations
numIters = 0
ACTIVATE_TUNING = True

def mainFOFOPIDOPT():
	global flag_FOPID_Schedule_Generation, ACTIVATE_TUNING, flag_FOPID_Computing_Output,the_fofopdt,und_fopid
	# Define the model
	the_fofopdt.K = 66.16
	the_fofopdt.L = 1.93
	the_fofopdt.T = 12.72
	the_fofopdt.alpha = 0.5

	the_fopid.Kp = -0.002934
	the_fopid.Ki = 0.01030
	the_fopid.Kd = 0.05335
	the_fopid.lamda = 0.9
	the_fopid.mu = 0.5

	# Controller parameters to be tuned
	und_fopid.Kp = 1 / the_fofopdt.K
	und_fopid.Ki = 1 / the_fofopdt.K
	und_fopid.Kd = 1 / the_fofopdt.K
	und_fopid.lamda = 0.9
	und_fopid.mu = 0.5

	# Specifications
	dspecs.wc = 0.1
	dspecs.pm = (60 * M_PI) / 180 # Convert to radians
	dspecs.opt_norm = 0.001 # Optimization termination criterion

	# Set approximation parameters and generate a FOPID controller
	params.wb = 0.001
	params.wh = 1000
	params.N = 5
	params.Ts = DT

	# If ACTIVATE_TUNING is on, use do the tuning here
	if ACTIVATE_TUNING:
		print("Controller Tuning Started")
		Do_FOPID_Optimization()
		print(the_fofopdt)
		print(the_fopid)
		print(und_fopid)
		print("Controller Tuninng Finished")

	# Generate the FOPID controller
	Generate_FOPID_Controller()

	while True:
		# If a FOPID generation has been scheduled,
		if flag_FOPID_Schedule_Generation:
		
			# check whether FOPID output sample is still being computed
			if not flag_FOPID_Computing_Output:
			
				# And if not, generate the new controller and clear flag
				Generate_FOPID_Controller()
				flag_FOPID_Schedule_Generation = False


# ********************************
# Coefficient computation function
# ********************************

def Compute_IIR_SOS_Oustaloup(zCoeffArray, pCoeffArray,	Kc, params, alpha):
	global zz,zp,EPS
	# Fetch approximation parameter values from params structure
	T = params.Ts
	wb = params.wb
	wh = params.wh
	N = params.N

	# Check input: upper frequency bound
	if (wh > (2 / T)):
		wh = (2 / T)

	# Calculate discrete-time transfer function zeros and poles
	omu = wh / wb
	# wb = -1 * wb  # The minus sign is important. Was the cause of many debugging issues when compare with matlab
	# global zz = [np.exp((omu ** ((kz + N + 0.5 * (1 - alpha)) / (2 * N + 1))) * wb * -T) for kz in range(-N, N + 1, 1)]  # Zeros
	# global zp = [np.exp((omu ** ((kp + N + 0.5 * (1 + alpha)) / (2 * N + 1))) * wb * -T) for kp in range(-N, N + 1, 1)]  # Poles
	for k in range(-N, N + 1, 1):
		w_kz = (omu ** ((k + N + 0.5 * (1 - alpha)) / (2 * N + 1))) * wb
		w_kp = (omu ** ((k + N + 0.5 * (1 + alpha)) / (2 * N + 1))) * wb

		#Discrete Mapping
		zz[k+N] = np.exp(-T*w_kz)
		zp[k+N] = np.exp(-T*w_kp)


	# K = wh ** alpha # gain
	# ttf = zpk2tf(w_kz, w_kp, K)
	#
	# return tf(ttf[0], ttf[1])
	#
	# ttf = zpk2tf(w_kz, w_kp, K)

	# Compute center frequency and correct gain
	wu = np.sqrt(wb*wh)
	Ks = wu ** alpha

	# Theta
	theta = np.cos(wu*T)

	# Compute absolute value ||H(z)|| at wu rad/s
	nk, dk, Ku = 1,1,1
	for k in range(0, 2*N + 1, 1):
		nk = 1 - (2 * zz[k] * theta) + zz[k]**2		# Zeros #nk = sqrt(1-2*zz[k]*theta+zz[k]*zz[k])
		dk = 1 - (2 * zp[k] * theta) + zz[k]**2		# Poles #dk = sqrt(1-2*zp[k]*theta+zp[k]*zp[k])

		if (nk >EPS) and (dk > EPS):
			Ku=Ku * nk/dk

	Ku = np.sqrt(Ku)

	# Compute the correct gain
	Kc = Ks / Ku

	# Compute second-order section form
	zCoeffArray[0][0] = -zz[2 * N]
	zCoeffArray[0][1] = 0.0
	pCoeffArray[0][0] = -zp[2 * N]
	pCoeffArray[0][1] = 0.0
	for k in range(N - 1, 0, -1):
		zCoeffArray[N - k][0] = -zz[2 * k + 1] - zz[2 * k]
		zCoeffArray[N - k][1] = zz[2 * k + 1] * zz[2 * k]
		pCoeffArray[N - k][0] = -zp[2 * k + 1] - zp[2 * k]
		pCoeffArray[N - k][1] = zp[2 * k + 1] * zp[2 * k]

	# Check stability
	isStable = IIR_SOS_Stability_Test(pCoeffArray)
	if isStable:
		print("Compute_IIR_SOS_Oustaloup computed was 'STABLE'")
	else:
		print("Compute_IIR_SOS_Oustaloup computed was 'UNSTABLE'")


def Generate_FOPID_Controller():
	# Reset FOPID_Ready flag
	global flag_FOPID_Ready
	flag_FOPID_Ready = False

	# Clear the memory banks
	Clear_IIR_Memory()

	# Generate the integrator
	print("Started Generating the integrator")
	Compute_IIR_SOS_Oustaloup(I_zsos, I_psos, KIc, params, 1.0 - the_fopid.lamda)
	print("Finished Generating the integrator")

	# Generate the differentiator
	print("Started Generating the differentiator")
	Compute_IIR_SOS_Oustaloup(D_zsos, D_psos, KDc, params, the_fopid.mu)
	print("Finished Generating the differentiator")

	# All done: Set FOPID_Ready flag
	flag_FOPID_Ready = True


def Clear_IIR_Memory():
	global s_IntMem, s_I, s_D
	for k in range(0,NSEC,1):
		s_I[k][0] = 0
		s_I[k][1] = 0
		s_D[k][0] = 0
		s_D[k][1] = 0
	
	s_IntMem = 0


def IIR_SOS_Stability_Test(pCoeffArray):
	global params
	# Use the "triangle rule" for SOS poles
	for i in range(0, params.N, 1):
		d1 = pCoeffArray[i][0]
		d2 = pCoeffArray[i][1]
		d1_abs = np.abs(d1)
		d2_abs = np.abs(d2)
		d1_abs_m_1 = d1_abs - 1.0
		d2_p_1 = d2 + 1.0

		if not (np.abs(d1) < d2 + 1.0):
			return False 			# Condition (1)
		if not (d2_abs < 1):
			return False           # Condition (2)

	# All coefficients are stable
	return True


# ***********************
# IIR filtering functions
# ***********************
def Do_IIR_Filtering(zCoeffArray, pCoeffArray, Kc, memArray, input):
	global params,s_IntMem
	# Assign local pointers
	b = zCoeffArray
	a = pCoeffArray
	s = memArray

	# Filter signal
	u_n = Kc * input
	for m in range(0,params.N, 1):
		y_n = u_n + s[m][0]
		s[m][0] = b[m][0] * u_n - a[m][0] * y_n + s[m][1]
		s[m][1] = b[m][1] * u_n - a[m][1] * y_n
		u_n = y_n
	

	# Check memory for under/overflow
	# NB! TODO: why do we do BOTH checks here? Should actually be part of DO_FOPID_Control function
	for k in range(0,NSEC, 1):
	
		s_I[k][0] = double_scale_saturation(s_I[k][0], EPS, NUMMAX)
		s_I[k][1] = double_scale_saturation(s_I[k][1], EPS, NUMMAX)
		s_D[k][0] = double_scale_saturation(s_D[k][0], EPS, NUMMAX)
		s_D[k][1] = double_scale_saturation(s_D[k][1], EPS, NUMMAX)
	
	s_IntMem = double_scale_saturation(s_IntMem, EPS, NUMMAX)

	# Assign output
	return u_n

def Do_FOPID_Control(err):
	global flag_FOPID_Computing_Output,in_Mem
	# Begin computing the output sample
	flag_FOPID_Computing_Output = True

	# Get the scaled ADC value
	inn = err

	# Check the input margin clear memory if exceeded
	if np.abs(inn - in_Mem) > INPUT_MARGIN:
		Clear_IIR_Memory()
	
	in_Mem = inn

	# FOPID computation for this sample
	foi_out = the_fopid.Ki*Do_IIR_Filtering(I_zsos, I_psos, KIc, s_I, inn)
	fod_out = the_fopid.Kd*Do_IIR_Filtering(D_zsos, D_psos, KDc, s_D, inn)
	s_IntMem += params.Ts * foi_out
	i_out = s_IntMem
	out = the_fopid.Kp*inn + i_out + fod_out

	# Done computing the output sample
	flag_FOPID_Computing_Output = False

	# Set the scaled value to DAC
	return Set_Scaled_Output(out)

def Set_Scaled_Output(out):

	# Saturate at top values
	scaled_out = out
	if (out > MAX_SCALED_OUT): scaled_out = MAX_SCALED_OUT
	if (out < MIN_SCALED_OUT): scaled_out = MIN_SCALED_OUT

	# Transform to normalized range -1...1
	scaled_out = (scaled_out - MED_SCALED_OUT) / (SCALE_FACTOR_OUT)

	# Return the scaled output
	return scaled_out

# Additional mathematical functions
def sign(x):
	return (x > 0) - (x < 0)


def double_scale_saturation(x, min, max):
	y = x
	if (np.abs(x) < min): y = 0
	if (np.abs(x) > max): y = sign(x) * max
	return y


# Helper trigonometric functions
def _sf(x):
	return np.sin((M_PI * x) / 2)

def _cf(x):
	return np.cos((M_PI * x) / 2)

# *****************************************
# Magnitude and phase response computations
# *****************************************

# Plant magnitude
def magng(w):
	return np.abs(the_fofopdt.K) / np.sqrt(1 + (the_fofopdt.T**2) * w**(2 * the_fofopdt.alpha) + 2 * the_fofopdt.T* (w** the_fofopdt.alpha) * _cf(the_fofopdt.alpha))


# Controller magnitude
def magnfopid(w):
	CR = und_fopid.Kp + (w**-und_fopid.lamda)*und_fopid.Ki*_cf(und_fopid.lamda) + (w**und_fopid.mu) * und_fopid.Kd * _cf(und_fopid.mu)
	CI = -(w**-und_fopid.lamda) * und_fopid.Ki * _sf(und_fopid.lamda) + (w**und_fopid.mu) * und_fopid.Kd * _sf(und_fopid.mu)

	return np.sqrt(CR**2 + CI**2)

# Plant phase
def phg(w):
	return -the_fofopdt.L*w - math.atan(((the_fofopdt.T)*_sf(the_fofopdt.alpha)) / (w**-the_fofopdt.alpha + (the_fofopdt.T * _cf(the_fofopdt.alpha))))

# Controller phase
def phfopid(w):
	CN = pow(w, und_fopid.lamda + und_fopid.mu)*und_fopid.Kd*_sf(und_fopid.mu) - und_fopid.Ki * _sf(und_fopid.lamda)
	CD = und_fopid.Ki*_cf(und_fopid.lamda) + pow(w, und_fopid.lamda) * (pow(w, und_fopid.mu)*und_fopid.Kd*_cf(und_fopid.mu) + und_fopid.Kp)

	return math.atan(CN / CD)


# The PSI functions and their derivatives
# Phase margin
def psi_pm(w):
	return (magng(w) * magnfopid(w)) - 1

def dpsi_pm(w):
	# System magnitude responses at w
	GM = magng(w)
	CM = magnfopid(w)

	# Helper values
	A11 = w**(-1 - (2 * und_fopid.lamda)) * (und_fopid.mu * (w**(2 * (und_fopid.lamda + und_fopid.mu))) * und_fopid.Kd**2 - und_fopid.lamda * und_fopid.Ki *
			( und_fopid.Ki + pow(w, und_fopid.lamda)*und_fopid.Kp*_cf(und_fopid.lamda)) + pow(w, und_fopid.lamda + und_fopid.mu) * und_fopid.Kd *
			( (und_fopid.lamda - und_fopid.mu) * und_fopid.Ki * _sf(und_fopid.lamda + und_fopid.mu)	+ und_fopid.mu * pow(w, und_fopid.lamda) * und_fopid.Kp * _cf(und_fopid.mu)))

	A1d = A11 / CM
	A2d = -((the_fofopdt.T*the_fofopdt.alpha*pow(w, the_fofopdt.alpha - 1)*	(the_fofopdt.T * pow(w, the_fofopdt.alpha) + _cf(the_fofopdt.alpha))) * GM) /\
		  (1 + pow(the_fofopdt.T, 2)*pow(w, 2 * the_fofopdt.alpha)	+ 2 * the_fofopdt.T*pow(w, the_fofopdt.alpha)*_cf(the_fofopdt.alpha))

	return (A1d * GM) + (CM * A2d)

# Gain margin
def psi_gm(w):
	return phg(w) + phfopid(w) + M_PI


def dpsi_gm(w):
	vd = dphg(w)
	B1d = dphfopid(w)
	return B1d + vd

	# Simple Newton's method test
def NRM_simple(w0,f,df):
	# Parameters
	N = 25
	k = 0
	x = w0
	xo = x
	NRM_eps = 0.001
	gamma = 1.5
	f_val = f(x)
	df_val = df(x)

	while (k+1 < N and f_val>NRM_eps):
		f_val = f(x)
		df_val = df(x)
		xo = x
		x = x - f_val / df_val
		if (x < 0):
			x = xo * gamma
	return x


# Derivative of plant phase response 
def dphg(w):
	return -the_fofopdt.L -	(the_fofopdt.alpha*the_fofopdt.T*_sf(the_fofopdt.alpha)) / (w*(2 * the_fofopdt.T*_cf(the_fofopdt.alpha) + pow(w, -the_fofopdt.alpha)
			+ pow(the_fofopdt.T, 2)*pow(w, the_fofopdt.alpha)))

# Derivative of controller phase response
def dphfopid(w):
	B11 = (und_fopid.lamda + und_fopid.mu)*und_fopid.Ki*_cf(und_fopid.lamda + und_fopid.mu - 1)	+ und_fopid.mu * pow(w, und_fopid.lamda) * und_fopid.Kp * _sf(und_fopid.mu)
	B20 = w * (w**( und_fopid.lamda + 2 * und_fopid.mu) * und_fopid.Kd** 2 + (w**-und_fopid.lamda) * (und_fopid.Ki ** 2) + (2 * und_fopid.Kp * und_fopid.Ki * _cf(und_fopid.lamda))
		+ ((w**und_fopid.lamda) * (und_fopid.Kp**2)) - 2 * pow(w, und_fopid.mu) * und_fopid.Kd * (und_fopid.Ki*_sf(und_fopid.lamda + und_fopid.mu - 1) -
																							  (w** und_fopid.lamda) * und_fopid.Kp * _cf(und_fopid.mu)))
	B10 = und_fopid.lamda * und_fopid.Kp * und_fopid.Ki * _sf(und_fopid.lamda) + pow(w, und_fopid.mu) * und_fopid.Kd * B11
	return B10 / B20

# FOPID optimization based on frequency-domain specifications
def Do_FOPID_Optimization():
	global OPT_MAX_ITER,und_fopid,numIters
	for k in range(0, OPT_MAX_ITER,1):
		numIters = k							# Number of used iterations
		b =  np.array([-kappa1(), -kappa2(), -kappa3()],dtype=float, order='C')	# Cost functions
		compute_specs_J()						# Update Jacobian
		dx = compute_cramer3(Jac, b)	# Solve the system of equations
		xn =  [und_fopid.Kp + dx.x1,			# Update the solution
					   und_fopid.Ki + dx.x2,
					   und_fopid.Kd + dx.x3]
		und_fopid.Kp = xn[0]							# Set the new solution
		und_fopid.Ki = xn[1]
		und_fopid.Kd = xn[2]

		# Compute norm and check it
		if (norm2_v3(b) < OPT_NORM):
			return

#Vector norm
def norm2_v3(b):
	return math.sqrt(b[0]**2 + b[1]**2 + b[2]**2)
	

# Critical frequency specification
def kappa1():
	return magnfopid(dspecs.wc) * magng(dspecs.wc) - 1

# Phase margin specification
def kappa2():
	return phfopid(dspecs.wc) + phg(dspecs.wc) + M_PI - dspecs.pm

# Phase response flatness at wc
def kappa3():
	return dphg(dspecs.wc) + dphfopid(dspecs.wc)

# Jacobian computation
def compute_specs_J():
	global und_fopid, dspecs,Jac
	# Enhance readability
	# |
	# Controller parameters
	Kp = und_fopid.Kp
	Ki = und_fopid.Ki
	Kd = und_fopid.Kd
	lam = und_fopid.lamda
	mu = und_fopid.mu
	# |
	# Critical frequency
	w = dspecs.wc
	A12 = pow(w, -2 * lam)*(-pow(w, lam + mu)*_sf(lam + mu - 1)*Kd + Ki + pow(w, lam)*_cf(lam)*Kp)
	A13 = pow(w, -lam + mu)*(pow(w, lam + mu)*Kd - _sf(lam + mu - 1)*Ki + pow(w, lam)*_cf(mu)*Kp)
	A2 = pow(w, 2 * (lam + mu))*pow(Kd, 2) + pow(Ki, 2) + 2 * pow(w, lam)*_cf(lam)*Ki*Kp + pow(w, 2 * lam)*pow(Kp, 2) + \
		 2 * pow(w, lam + mu)*Kd*(-_sf(lam + mu - 1)*Ki + pow(w, lam)*_cf(mu)*Kp)
	A3 = pow(A2, 2)
	A31 = pow(w, lam - 1)*(	mu*pow(w, 3 * (lam + mu))*_sf(mu)*pow(Kd, 3) - pow(w, 2 * (lam + mu))*(2 * mu*_sf(lam) +
						lam * _sf(lam + 2 * mu))*pow(Kd, 2)*Ki + lam * _sf(lam)*Ki*(pow(Ki, 2) - pow(w, 2 * lam)*pow(Kp, 2)) -
						pow(w, lam + mu)*Kd* ((2 * lam*_sf(mu) + mu * _sf(2 * lam + mu))*pow(Ki, 2) + 2 * (lam + mu)*pow(w, lam)*_sf(lam + mu - 1)*Ki*Kp +
						mu * pow(w, 2 * lam)*_sf(mu)*pow(Kp, 2)))

	A32 = pow(w, lam - 1)*((lam + mu)*pow(w, 2 * lam + 3 * mu)*_cf(lam + mu - 1)*pow(Kd, 3) +
		pow(w, 2 * (lam + mu))*(2 * (lam + mu)*_sf(lam) + lam * _sf(lam + 2 * mu))*pow(Kd, 2)*Kp +	lam * _sf(lam)*Kp*(-pow(Ki, 2) + pow(w, 2 * lam)*pow(Kp, 2)) +
		pow(w, mu)*Kd*(-(lam + mu)*_cf(mu + lam - 1)*pow(Ki, 2) - 2 * mu*pow(w, lam)*_sf(mu)*Ki*Kp + pow(w, 2 * lam)*(2 * lam*_cf(lam + mu - 1) + (lam + mu)*_sf(lam - mu))*pow(Kp, 2)))

	A33 = pow(w, lam + mu - 1)*((lam + mu)*_cf(lam + mu - 1)*pow(Ki, 3) - 2 * lam*pow(w, 2 * lam + mu)*_sf(lam)*Kd*Ki*Kp +
			pow(w, lam) * (2 * (lam + mu)*_sf(mu) + mu * _sf(2 * lam + mu))*pow(Ki, 2)*Kp + pow(w, 2 * lam)*(2 * mu*_cf(lam + mu - 1) - (lam + mu)*_sf(lam - mu))*Ki*pow(Kp, 2)
		+ mu * pow(w, 3 * lam)*_sf(mu)*pow(Kp, 3) - pow(w, 2 * (lam + mu))*pow(Kd, 2)*((lam + mu)*_cf(lam + mu - 1)*Ki + mu * pow(w, lam)*_sf(mu)*Kp))

	ACR = Kp + pow(w, -lam)*_cf(lam)*Ki + pow(w, mu)*_cf(mu)*Kd

	magndiv = magng(w) / magnfopid(w)

	Jac[0][0] = magndiv * ACR
	Jac[0][1] = magndiv * A12
	Jac[0][2] = magndiv * A13
	Jac[1][0] = (pow(w, lam)*(-pow(w, lam + mu)*_sf(mu)*Kd + _sf(lam)*Ki)) / A2
	Jac[1][1] = -(pow(w, lam)*(pow(w, mu)*_cf(mu + lam - 1)*Kd + _sf(lam)*Kp)) / A2
	Jac[1][2] = (pow(w, lam + mu)*(_cf(lam + mu - 1)*Ki + pow(w, lam)*_sf(mu)*Kp)) / A2
	Jac[2][0] = A31 / A3
	Jac[2][1] = A32 / A3
	Jac[2][2] = A33 / A3


# 3x3 matrix: determinant computation
def compute_det3(A):
	return A[0][0] * (A[1][1] * A[2][2] - A[1][2] * A[2][1]) - \
		   A[0][1] * (A[1][0] * A[2][2] - A[1][2] * A[2][0]) + \
		   A[0][2] * (A[1][0] * A[2][1] - A[1][1] * A[2][0])


# 3x3 matrix: Cramer rule
def compute_cramer3(A,  b):
	if isinstance(A, np.ndarray):
		if A.shape == (3,3):
			pass
		else:
			ValueError("FOPIDControllerPrototypeMarkII.compute_cramer3: A is ndarray but shape is not equal to (3,3)")
	else:
		ValueError("FOPIDControllerPrototypeMarkII.compute_cramer3: A is not a numpy array")
	if isinstance(b, np.ndarray):
		if b.size == 3:
			pass
		else:
			ValueError("FOPIDControllerPrototypeMarkII.compute_cramer3: b is ndarray but shape is not equal to (3,1)")
	else:
		ValueError("FOPIDControllerPrototypeMarkII.compute_cramer3: b is not a numpy array")

	# EPS
	my_eps = 1e-7

	# The solution
	system_solution = Dict(sle_sol)
	system_solution.exists = False

	# Compute the determinant of A and check it
	detA = compute_det3(A)
	#TODO: if np.abs(detA) == compute_det3(A) replace compute_det3 function
	print(detA)
	print(np.linalg.det(A))
	if detA == np.linalg.det(A):
		print("TODO: if np.abs(detA) == compute_det3(A) replace compute_det3 function")

	if (np.abs(detA) < my_eps):
		return system_solution
	

	# Solve the system
	mat_x1 =  np.array([[b[0],A[0][1],A[0][2]],[b[1],A[1][1],A[1][2]],[b[2],A[2][1],A[2][2]]],dtype=float)
	# mat_x1 = mat_x1.reshape((3,3))
	mat_x2 =  np.array([[A[0][0],b[0],A[0][2]],[A[1][0],b[1],A[1][2]],[A[2][0],b[2],A[2][2]]],dtype=float)
	# mat_x2 = mat_x2.reshape((3,3))
	mat_x3 =  np.array([[A[0][0],A[0][1],b[0]],[A[1][0],A[1][1],b[1]],[A[2][0],A[2][1],b[2]]],dtype=float)
	# mat_x3 = mat_x3.reshape((3, 3))

	system_solution.x1 = compute_det3(mat_x1) / detA
	system_solution.x2 = compute_det3(mat_x2) / detA
	system_solution.x3 = compute_det3(mat_x3) / detA
	system_solution.exists = True
	return system_solution


if __name__ == '__main__':
	mainFOFOPIDOPT()