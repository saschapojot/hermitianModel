from consts import *

def x(n):
    return 2 * n * omegaF


def y(n):
    return (2 * n + 1) * omegaF


xmValsAll = [x(m) for m in range(0, N)]
ymValsAll = [y(m) for m in range(0, N)]



def zerothVector(firstVec,secondVec,g):
    '''
    the 0th vector for RK2
    :param firstVec: wvfcnt at even sites
    :param secondVec: wvfcnt at odd sites
    :param g: nonlinearity
    :return: K0, M0
    '''

    K0=[]
    M0=[]

    K0.append(
        -1j*g*secondVec[0]**2*np.conj(firstVec[0])
    )

    #n=1,...,N-1
    for n in range(1,N):
        K0.append(
            -1j*g*(secondVec[n]**2*np.conj(firstVec[n])+secondVec[n-1]**2*np.conj(firstVec[n]))
        )
    #n=0,...,N-2
    for n in range(0,N-1):
        M0.append(
            -1j*g*(firstVec[n]**2*np.conj(secondVec[n])+firstVec[n+1]**2*np.conj(secondVec[n]))
        )
    M0.append(
        -1j*g*firstVec[N-1]**2*np.conj(secondVec[N-1])
    )

    return K0,M0


def firstVector(firstVec, secondVec,g,K0,M0):
    '''
     the 1st vector for RK2
    :param firstVec: wvfcnt at even sites
    :param secondVec: wvfcnt at odd sites
    :param g: nonlinearity
    :param K0: K0 coef vector
    :param M0: M0 coef vector
    :return: K1, M1
    '''
    K1=[]
    M1=[]

    K1.append(
        -1j*g*(secondVec[0]+1/2*dt*M0[0])**2\
        *(np.conj(firstVec[0])+1/2*dt*np.conj(K0[0]))
    )
    #n=1,...,N-1
    for n in range(1,N):
        K1.append(
            -1j*g*(
                (secondVec[n]+1/2*dt*M0[n])**2
                *(np.conj(firstVec[n])+1/2*dt*np.conj(K0[n]))
                +(secondVec[n-1]+1/2*dt*M0[n-1])**2
                *(np.conj(firstVec[n])+1/2*dt*np.conj(K0[n]))
            )
        )
    #n=0,..., N-2
    for n in range(0,N-1):
        M1.append(
            -1j*g*(
                (firstVec[n]+1/2*dt*K0[n])**2
                *(np.conj(secondVec[n])+1/2*dt*np.conj(M0[n]))
                +(firstVec[n+1]+1/2*dt*K0[n+1])**2
                *(np.conj(secondVec[n])+1/2*dt*np.conj(M0[n]))
            )
        )

    M1.append(
        -1j*g*(firstVec[N-1]+1/2*dt*K0[N-1])**2
        *(np.conj(secondVec[N-1])+1/2*dt*np.conj(M0[N-1]))
    )

    return K1,M1


def oneStepRk2(firstVec,secondVec,g):
    '''

    :param firstVec: wvfcnt at even sites
    :param secondVec: wvfcnt at odd sites
    :param g: nonlinearity
    :return: out1, out2
    '''
    K0, M0 = zerothVector(firstVec, secondVec, g)
    K1, M1 = firstVector(firstVec, secondVec, g, K0, M0)

    out1 = []
    out2 = []
    for n in range(0, N):
        out1.append(firstVec[n] + 1 / 4 * dt * (K0[n] + K1[n]))
        out2.append(secondVec[n] + 1 / 4 * dt * (M0[n] + M1[n]))
    return out1, out2


def phi1A2(deltaTau, uQ, alphaVec):
    '''

    :param deltaTau: time step
    :param uQ: u value at time (q+1/2)dt
    :param alphaVec: input wvfcnt at even sites
    :return:
    '''
    for m in range(0, N):
        alphaVec[m] *= np.exp(-1j * deltaTau * (xmValsAll[m] + uQ))


def phi1B2(deltaTau, uQ, betaVec):
    '''

    :param deltaTau: time step
    :param uQ: u value at time step (q+1/2)dt
    :param betaVec: input wvfcnt at odd sites
    :return:
    '''
    for m in range(0, N):
        betaVec[m] *= np.exp(-1j * deltaTau * (ymValsAll[m] - uQ))


def phi1C1(deltaTau, vQ, alphaVec, betaVec):
    '''

    :param deltaTau: time step
    :param vQ: v value at time step (q+1/2)dt
    :param alphaVec: input wvfcnt at even sites
    :param betaVec: input wvfcnt at odd sites
    :return:
    '''
    for m in range(0, N):
        alphaVec[m] -= 1j * deltaTau * vQ * betaVec[m]


def phi1D1(deltaTau, vQ, alphaVec, betaVec):
    '''

    :param deltaTau: time step
    :param vQ: v value at time step (q+1/2)dt
    :param alphaVec: input wvfcnt at even sites
    :param betaVec: input wvfcnt at odd sites
    :return:
    '''
    for m in range(0, N):
        betaVec[m] -= 1j * deltaTau * vQ * alphaVec[m]


def phi21(deltaTau, wQ, alphaVec, betaVec):
    '''

    :param deltaTau: time step
    :param wQ: w value at (q+1/2)dt
    :param alphaVec: input wvfcnt at even sites
    :param betaVec: input wvfcnt at odd sites
    :return:
    '''
    for j in range(1, N):
        alphaVec[j] -= 1j * deltaTau * wQ * betaVec[j - 1]


def phi22(deltaTau, wQ, alphaVec, betaVec):
    '''

    :param deltaTau: time step
    :param wQ: w value at time step (q+1/2)dt
    :param alphaVec: input wvfcnt at even sites
    :param betaVec: input wvfcnt at odd sites
    :return:
    '''
    for j in range(0, N - 1):
        betaVec[j] -= 1j * deltaTau * wQ * alphaVec[j + 1]


def zeta1(uQ, vQ, alphaVec, betaVec):
    '''
    composite mapping zeta1
    :param uQ: u value at time (q+1/2)dt
    :param vQ: v value at time (q+1/2)dt
    :param alphaVec: input wvfcnt at even sites
    :param betaVec: input wvfcnt at odd sites
    :return:
    '''
    phi1B2(dt / 2, uQ, betaVec)
    phi1A2(dt / 2, uQ, alphaVec)
    phi1C1(dt / 2, vQ, alphaVec, betaVec)
    phi1D1(dt, vQ, alphaVec, betaVec)
    phi1C1(dt / 2, vQ, alphaVec, betaVec)
    phi1B2(dt / 2, uQ, betaVec)
    phi1A2(dt / 2, uQ, alphaVec)


def zeta2(wQ, alphaVec, betaVec):
    '''
    composite mapping zeta2
    :param wQ: w value at time step (q+1/2)dt
    :param alphaVec: input wvfcnt at even sites
    :param betaVec: input wvfcnt at odd sites
    :return:
    '''
    phi21(dt / 4, wQ, alphaVec, betaVec)
    phi22(dt / 2, wQ, alphaVec, betaVec)
    phi21(dt / 4, wQ, alphaVec, betaVec)


def expmidtF(q, alphaVec, betaVec):
    '''
    linear step in S2
    :param q: time step q
    :param alphaVec: input wvfcnt at even sites
    :param betaVec: input wvfcnt at odd sites
    :return:
    '''

    t = (q + 1 / 2) * dt
    uq = u(t)
    vq = v(t)
    wq = w(t)
    zeta2(wq, alphaVec, betaVec)
    zeta1(uq, vq, alphaVec, betaVec)
    zeta2(wq, alphaVec, betaVec)

def oneStepS2(firstVec, secondVec, q, g):
    '''

    :param firstVec: input wvfcnt at even sites
    :param secondVec: input wvfcnt at odd sites
    :param q: time step q
    :param g: nonlinearity
    :return:
    '''
    alphaVec, betaVec = oneStepRk2(firstVec, secondVec, g)
    expmidtF(q, alphaVec, betaVec)
    r, s = oneStepRk2(alphaVec, betaVec, g)
    return r, s


def meanX(psiQ):
    '''

    :param psiQ: wvfcnt at time step q
    :return: mean position at time q
    '''
    xOut=0
    #norm2=0
    for j in range(0,len(psiQ)):
        xOut+=j*np.abs(psiQ[j])**2
        #norm2+=np.abs(psiQ[j])**2
    return xOut
    #/norm2