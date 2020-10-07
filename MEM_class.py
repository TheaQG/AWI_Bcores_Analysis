import numpy as np
import matplotlib.pyplot as plt

class MEM():
    def __init__(self, t_data, y_data, M):
        self.t_data = t_data
        self.y_data = y_data
        self.M = M

        return

    def __call__(self, t_data, y_data, M, N = None, view = True, print_coefs = False):

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

        return power

    def Coef(self, M):
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
