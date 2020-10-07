import scipy as sp
from scipy import stats
from scipy import signal
from scipy import fft
from scipy import io
from scipy import interpolate
from scipy import optimize
from scipy import linalg
from scipy import integrate

import numpy as np
from matplotlib import pyplot as plt


class SyntheticData_Gen():
    def __init__(self, AR1_coef, AR1_var, AR1_dt, AR1_N, diff_len, dt_sample, meas_noise):
        '''
            Initializes the class with the following arguments:

                AR1_coef:       [float] AR coeffient(only 1, bc AR-1). [0,1)
                AR1_var:        [float] Variance of the initial AR1 series (set as a random process)
                AR1_dt:         [float] Spacing of z array (depth array)
                AR1_N:          [int] Number of points to be created by AR-1 process
                diff_len:       [float] Diffusion length to simulate [m]
                dt_sample:      [float] Discrete sampling interval (Must be more than 3 times smalle than diff_len)
                meas_noise:     [float] Standard deviation of the measurement noise added to the final time series [permil]
        '''

        self.AR1_coef = AR1_coef
        self.AR1_var = AR1_var
        self.AR1_dt = AR1_dt
        self.AR1_N = AR1_N
        self.diff_len = diff_len
        self.dt_sample = dt_sample
        self.meas_noise = meas_noise

        return

    def __call__(self):
        '''
            Performs the steps needed to generate synthetic data. Computes the raw synthetic
            AR1 process, then resamples it to the wanted sample size, and finally computes the
            diffused signal through convolution of signal with gaussian window.
            Returns signal data for all three steps.
        '''

        zAR1, xAR1 = self.synthetic_AR1()

        zAR1_dt, xAR1_dt = self.sample(z = zAR1, x = xAR1, dt = self.dt_sample, meas_noise = self.meas_noise)

        zConv, xConv = self.diffuse(z = zAR1_dt, x = xAR1_dt, diff_len = self.diff_len)

        return zAR1, xAR1, zAR1_dt, xAR1_dt, zConv, xConv

    def synthetic_AR1(self):
        '''
            Create a random autoregressive (first order) process.

            Arguments:
            ----------
                Self

            Returns:
            --------
                z_arr:      [array of floats] Time/depth array
                x_AR:       [array of floats] Array containing values of the AR-1 process
        '''


        # Define the maximum independent variable value
        z_max = self.AR1_N * self.AR1_dt
        # Create an array containing independent variables
        z_arr = np.arange(0, z_max, self.AR1_dt)
        # Generate an array of normally distributed random values to work from
        x = np.random.randn(self.AR1_N) * np.sqrt(self.AR1_var)
        x_AR = np.zeros_like(x)
        # Set the order
        AR_order = 1

        # Compute the autoreressive process
        for j in range(AR_order, len(x_AR) - AR_order - 1):
            # Define the previous values to compute present value from
            prevVals = x[j - AR_order : j]
            # Define the sum coefficients
            sumCoefs = np.zeros(AR_order) + self.AR1_coef
            # Compute the present value on the basis of the preceeding
            x_AR[j] = np.sum(prevVals * sumCoefs) + x[j]

        return z_arr, x_AR

    def smooth(self, z, x, window_len = 11, window = 'hamming'):
        '''
            Function which smooths the data with a specific window type of a given size.
            Based on convolution of a scaled specific window with the signal.
            The signal is adjusted to minimize transient parts in the beginning and ending
            of the array. This is done by introducing reflected copies(of same size as the window)
            of the signal at both ends of the array.

            Arguments:
            ----------
                window_len:     [int] Length of the smoothing window. Should be odd
                window:         [string] Type of window. ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']
                                        Flat will return a moving average smoothing.
            Returns:
            --------
                y:              [array of floats] Array of the smoothed signal


        '''

        #Test input params and raise errors if not acceptable
        if x.ndim != 1:
            raise ValueError("Function 'smooth' only accepts 1D arrays")
        if x.size < window_len:
            raise ValueError("Input array must be larger than window size")
        if window_len < 3:
            return z, x
        if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman', 'gaussian']:
            raise ValueError("'window' must be one of:\n['flat', 'hanning', 'hamming', 'bartlett', 'blackman', 'gaussian']")

        # Generate 'extra' windows at beginning and end of array (elongate array)
        x_mov = np.concatenate([2*x[0] - x[window_len:1:-1], x, 2*x[-1] - x[-1:-window_len:-1]])

        # Choose what window type to smooth with
        if window == 'flat': #moving average
            w = np.ones(window_len, 'd')
        elif window == 'gaussian': #gaussian
            w = eval('sp.signal.gaussian(5*window_len, window_len)')
        else: #any other typesample
            w = eval('np.' + window + '(window_len)')

        # Evaluate convolution btw. signal and window
        y = np.convolve(w/w.sum(), x_mov, mode='same')
        # Get only the wanted values (without the 'extra' windows)
        y = y[window_len - 1 : -window_len + 1]

        return z, y

    def smooth_step(self, x, window_len = 11, window = 'flat'):
        '''
            Similar to 'smooth'-function, except that the window slides a number of points equal
            to its own size at every step.
            Supports only the moving average window ('flat').

            Arguments:
            ----------
                window_len:     [int] Length of the smoothing window. Should be odd
                window:         [string] Type of window. ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']
                                        Flat will return a moving average smoothing.
            Returns:
            --------
                y:              [array of floats] Array of the smoothed signal with size floor(N/window_len)

        '''

        #Test input params and raise errors if not acceptable
        if x.ndim != 1:
            raise ValueError("Function 'smooth' only accepts 1D arrays")
        if x.size < window_len:
            raise ValueError("Input array must be larger than window size")

        # Set sizes of arrays needed
        N = len(x)
        M = window_len
        N_prime = np.int(np.floor(N/M))
        y = np.zeros(N_prime)

        # Calculate mean of window (step)
        for j in range(N_prime):
            y[j] = np.mean(x[j*M : M*(j+1) - 1])

        return y

    def diffuse(self, z, x, diff_len):
        '''
            Diffuse a time/depth series x_ar given a specific diffusion length

            Arguments:
            ----------
                diff_len:       [float] Diffusion length in [m]

            Returns:
                z:              [array of floats] The depth series
                x_AR_conv:      [array of floats] The diffused AR-1 process
            --------


        '''

        #Test input params and raise errors if not acceptable
        if np.size(z) != np.size(x):
            raise ValueError("Z, Y arrays must be of same size.")

        dt_fine = z[1] - z[0]
        N_fine = len(z)

        try:
            win_len = np.int(np.around(diff_len/dt_fine))
            x_AR_conv = self.smooth(z, x, window_len = win_len, window = 'gaussian')[1]
        except:
            print("Error calculating the convolution")

        return z, x_AR_conv

    def sample(self, z, x, dt = 0.005, meas_noise = 0.07):
        '''
            Performs discrete sampling of time series x.

            Arguments:
            ----------
                dt:             [float] The discrete sampling interval
                meas_noise:     [float] Standard deviation of the measurement noise added to the final time series

            Returns:
            --------
                z_dt:           [array of floats] Discrete sampled independent variable
                x_dt:           [array of floats] Discrete sampled dependent variable
        '''

        dt_fine = z[1] - z[0]
        if dt_fine > dt:
            print('The discrete sampling  width is finer than the initial resolution. dt is set to dt_fine')
            dt = dt_fine

        z_dt = np.arange(0, z[-1], dt)
        x_dt = self.smooth_step(x = x, window_len = np.int(dt/dt_fine))
        x_dt = x_dt + np.random.randn(len(x_dt)) * meas_noise

        if len(x_dt) > len(z_dt):
            x_dt = x_dt[:len(z_dt)]
        elif len(x_dt) < len(z_dt):
            z_dt = z_dt[:len(x_dt)]
        else:
            pass


        return z_dt, x_dt
