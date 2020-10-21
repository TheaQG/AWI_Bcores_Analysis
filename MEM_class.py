import numpy as np
import matplotlib.pyplot as plt

class MEM():

    """
        Class to determine the spectral equivalent of a time series through
        Burg's Maximum Entropy Method (MEM), Burg 1967.
    """

    """
        Methods available:

            __init__():
                    Initializes the class given three arguments, t, y, M.

            __call__():
                    Calling function, returns the frequency and spectral
                    signal of the time series, (w, P).

            Coef(self, M):
                    Calculates the maximum entropy coefficients on the basis
                    of the number of poles passed.

            Power(self, N, M):
                    Calculates the power spectrum(spectral signal) on the
                    basis of MEM coefficients, given the number of points
                    wanted and the number of poles passed.

    """
    def __init__(self, t_data, y_data, M):
        self.t_data = t_data
        self.y_data = y_data
        self.M = M

        return

    def __call__(self, t_data, y_data, M, N = None, view = True, print_coefs = False):
        """

                Arguments:
                ----------
                    t_data:         [array of floats] Time series x data.
                    y_data:         [array of floats] Time series y data.
                    M:              [int] Number of poles to perform MEM with
                    N:              [int] Default: None. Number of points to perform MEM with.
                    view:           [bool] Default: True. View plot of PSD?
                    print_coefs:    [bool] Defulat: False. Print MEM coeffients?

                Returns:
                --------
                    power[0]:       [array of floats] Frequencies. Transform of t_data.
                    power[1]:       [array of floats] Power (PSD). Transform of y_data.
        """
        if N == None:
            N = np.size(y_data)

        mem = MEM(t_data, y_data, M)
        power = mem.Power(N, M)

        if view == True:
            plt.ion()
            plt.clf()
            plt.subplot(211)
            line1 = plt.plot(t_data, y_data)
            plt.subplot(212)
            line2 = plt.semilogy(power[0], power[1])

        if print_coefs == True:
            print(MEM.Coef(M, N)[1])

        return power[0], power[1]

    def Coef(self, M):
        """
            On the basis of time seriess y_data and number of poles passed, M,
            computes the coefficients according to the MEM spectral transformation
            method.

                Arguments:
                ----------
                    M:              [int] Number of poles to use.

                Returns:
                --------
                    P:
                    a:
                    ref:
        """

        y = self.y_data
        N = np.size(y)
        y -= np.average(y)

        P0 = np.sum(y**2)/N

        m = 0

        b1 = np.copy(y[0:-1])
        b2 = np.copy(y[1:N])

        nom = np.sum(b1[0:N-1]*b2[0:N-1])
        denom = np.sum(b1[0:N-1]**2 + b2[0:N-1]**2)

        a = np.array(())
        ref = np.array(())

        a = np.append(a,2*nom/denom)
        ref = np.append(ref,2*nom/denom)

        P = np.array(())

        P = np.append(P,P0*(1 - a[m]**2))

        for m in range(1,M):
            aa = np.array(())
            for k in range(0,m):
                aa = np.append(aa,a[k])

            for k in range(0,N - (m+1)):
                b1[k] = b1[k] - aa[m-1]*b2[k]
                b2[k] = b2[k+1] - aa[m-1]*b1[k+1]

            nom = np.sum(b1[0:N - (m+1)] * b2[0:N - (m+1)])
            denom = np.sum(b1[0:N - (m+1)]**2 + b2[0:N - (m+1)]**2)

            a = np.append(a, 2*nom/denom)
            ref = np.append(ref, 2*nom/denom)
            P = np.append(P,P[m-1]*(1 - a[m]**2))

            for kk in range(0,m):
                a[kk] = aa[kk] - a[m]*aa[m-kk-1]

        return (P, a, ref)

    def Power(self, N, M):
        """
            On the basis of number of points wanted and number of poles passed,
            computes the PSD of the time series through Burg's MEM.

                Arguments:
                ----------
                    N:              [int] Number of points wanted.
                    M:              [int] Number of poles passed.

                Returns:
                --------
                    f:              [array of floats] Computed frequencies.
                    Pf:             [array of floats] Computed power (PSD).
        """

        P = self.Coef(M)[0]
        a = self.Coef(M)[1]

        Pm = P[-1]
        dt = self.t_data[1] - self.t_data[0]

        i = 0. + 1.j

        nq = 1/(2*dt)

        f = np.arange(N/2)*2*nq/N

        nomF = Pm*dt

        sumF = 0
        for k in range(1,len(a)+1):
            sumF += a[k-1]*np.exp(-2*np.pi*i*f*k*dt)

        denomF =np.abs(1-sumF)**2

        Pf = nomF/denomF

        return (f, Pf)
