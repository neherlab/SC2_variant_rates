import numpy as np
from scipy.stats import nbinom, poisson
from scipy.optimize import minimize
import matplotlib.pyplot as plt
from collections import defaultdict

TINY = 1e-100
LogTINY = np.log(TINY)
birth_rate = 1e-6
logCutoff = 50

def log_propagator(n, m, rho, k, mu):
    '''
    nvec = np.arange(100)
    k,rho = 0.1, 1.1
    pn = propagator(nvec, 5, rho,k)
    avg = np.sum(pn*nvec) # should be rho*5
    '''
    p = k/(k+rho)
    if m:
        return nbinom.logpmf(n,m*k,p)
    else:
        tmp = LogTINY*np.ones_like(n)
        tmp[n==0] = -mu
        tmp[n==1] = np.log(mu)
        return tmp

def log_propagator_mvec(n, m, rho, k, mu):
    p = k/(k+rho)
    res = LogTINY*np.ones_like(m)
    res[m>0] = nbinom.logpmf(n,m[m>0]*k,p)
    if n==0:
        res[m==0] = -mu
    elif n==1:
        res[m==0] = np.log(mu)

    return res

def propagator(n, m, rho, k):
    return np.exp(log_propagator(n,m,rho,k))


def log_sampling_prob(n_samples, n_cases, eps):
    mean_density = eps*n_cases
    return poisson.logpmf(n_samples, mean_density)


def back_propagate(samples_t, logR_tp1, rho, k, eps, mu, nmin=0, nmax=1000, daughters=False, daughter_rate=0):
    nvals = np.arange(nmin, nmax)
    if logR_tp1 is None:
        logR_t = np.array([nvals, log_sampling_prob(samples_t, nvals, eps)]).T
    else:
        peak = logR_tp1[:,1].max()
        logR_t = np.array([(n_t, log_sampling_prob(samples_t,n_t, eps)
                       + peak + np.log(np.sum(np.exp(logR_tp1[:,1] - peak + log_propagator(logR_tp1[:,0], n_t, rho, k, mu)))))
                       for n_t in nvals])
    if daughters:
        logR_t[:,1] += np.log(1-np.exp(-daughter_rate*logR_t[:,0]))

    peak = logR_t[:,1].max()

    return logR_t[logR_t[:,1]>peak - logCutoff]

def fwd_propagate(samples_tm1, logQ_tm1, rho, k, eps, mu, nmin=0, nmax=1000, daughters=False, daughter_rate=0):
    nvals = np.arange(nmin, nmax)
    if logQ_tm1 is None:
        logQ_t = np.array([nvals, LogTINY*np.ones_like(nvals)]).T
        logQ_t[nvals==0,1] = 0
    else:
        peak = logQ_tm1[:,1].max()
        logQ_t = np.array([(n_t, log_sampling_prob(samples_tm1,n_t, eps)
                       + peak + np.log(np.sum(np.exp(logQ_tm1[:,1] - peak + log_propagator_mvec(n_t, logQ_tm1[:,0], rho, k, mu)))))
                       for n_t in nvals])
    if daughters:
        logQ_t[:,1] += np.log(1-np.exp(-daughter_rate*logQ_t[:,0]))

    peak = logQ_t[:,1].max()
    return logQ_t[logQ_t[:,1]>peak - logCutoff]


def back_trace(data):
    sorted_data = data.sort_index(ascending=False)
    logR = {}
    for i, (t, row) in enumerate(sorted_data.iterrows()):

        prev = logR[sorted_data.iloc[i-1].name] if i else None
        logR[t] = back_propagate(row.cases, prev, row.rho, row.k, row.eps, row.mu, nmin=row.cases, nmax=3*(row.cases+3)/eps)

    return logR


def fwd_trace(data):
    sorted_data = data.sort_index(ascending=True)
    logQ = {}
    for i, (t, row) in enumerate(sorted_data.iterrows()):
        prev = logQ[sorted_data.iloc[i-1].name] if i else None
        prev_cases = sorted_data.iloc[i-1].cases if i else None
        logQ[t] = fwd_propagate(prev_cases, prev, row.rho, row.k, row.eps, row.mu, nmin=row.cases, nmax=3*(row.cases+3)/eps)

    return logQ

def negLogLH(data, t0):
    logR = back_trace(data)
    return -logR[t0][0,1]

def cost(x, data, t0, fields):
    print(x)
    for y, n in zip(x,fields):
        data[n] = y**2
    #data['eps'] = x[0]
    return negLogLH(data, t0)


def marginal_distribution(logQ, logR):
    nval_range = int(max(logQ[0,0],logR[0,0])), int(min(logQ[-1,0], logR[-1,0])+1)
    offset = int(min(logQ[0,0],logR[0,0]))
    indR = (logR[:,0]>=nval_range[0]) & (logR[:,0]<nval_range[1])
    indQ = (logQ[:,0]>=nval_range[0]) & (logQ[:,0]<nval_range[1])
    tmplogP = logR[indR,1] + logQ[indQ,1]
    return np.array([logR[indR,0], tmplogP]).T


if __name__=="__main__":
    import pandas as pd
    from scipy.optimize import minimize

    parent_cases =   [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 5, 10, 40, 70, 100]
    daughter_cases = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 2, 8, 3, 6, 20, 30, 50]
    rho = 1.3
    k=.2
    eps=0.05
    data = pd.DataFrame({i:{'cases':parent_cases[i], 'rho':rho, 'eps':eps, 'k':k, 'mu':birth_rate} for i in range(len(parent_cases))}).T
    sol = minimize(cost, [rho], (data, 0, ['rho']))
    rho = sol['x'][0]

    data_daughter = pd.DataFrame({i:{'cases':daughter_cases[i], 'rho':rho, 'eps':eps, 'k':k, 'mu':TINY} for i in range(len(parent_cases))}).T

    logR_daughter = back_trace(data_daughter)
    logR = back_trace(data)

    daughter_seeding = {t:lR[lR[:,0]==1,1][0] for t,lR in logR_daughter.items() if lR[0,0]<=1}
    logR_wdaughter = {}
    for seeding_time, lR_daughter in daughter_seeding.items():
        lR_wd = {}
        for t,lR in logR.items():
            if t>=seeding_time:
                lR_wd[t] = lR
            else:
                row = data.loc[t]
                lR_wd[t] = back_propagate(row.cases, lR_wd[t+1], row.rho, row.k, row.eps, row.mu,
                                          nmin=row.cases, nmax=3*(row.cases+3)/eps,
                                          daughters=(t+1==seeding_time), daughter_rate=0.02)
        logR_wdaughter[seeding_time] = lR_wd

    print(logR[0][0,1])
    for t in logR_wdaughter:
        print(t,logR_wdaughter[t][0][0,1]+daughter_seeding[t], daughter_seeding[t])



    plt.figure()
    for t in logR:
        plt.plot(logR[t][:,0]+0.5, logR[t][:,1] - logR[t][:,1].max(), label=f'cases={data.at[t,"cases"]}')

    plt.xscale('log')
    plt.legend()

    logQ = fwd_trace(data)
    plt.figure()
    for t in logQ:
        plt.plot(logQ[t][:,0]+0.5, logQ[t][:,1] - logQ[t][:,1].max(), label=f'cases={data.at[t,"cases"]}')

    plt.xscale('log')
    plt.legend()



    time_points = range(1,len(data))
    logP = {}
    for i in time_points:
        t = data.iloc[i].name
        logP[t] = marginal_distribution(logQ[t], logR[t])
    plt.xscale('log')
    plt.legend()


    daughter_seeding = {t:lR[lR[:,0]==1,1][0] for t,lR in logR_daughter.items() if lR[0,0]<=1}
    logP_wdaughter = {}
    peak_val = -np.inf
    for i in time_points:
        tmp = {}
        tm1 = data.iloc[i-1].name
        t = data.iloc[i].name
        for seeding_time, lR_daughter in daughter_seeding.items():
            tmp[seeding_time] = marginal_distribution(logQ[tm1], logR_wdaughter[seeding_time][t], data.loc[t])
            tmp[seeding_time][:,1] += lR_daughter
            peak_val = max(peak_val, tmp[seeding_time][:,1].max())
            # P = np.exp(tmp[seeding_time][:,1] - tmp[seeding_time][:,1].max())
            # P/=P.sum()
            # plt.plot(tmp[seeding_time][:,0]+0.5, P, c=f"C{i%10}", lw=1)

        logP_wdaughter[t] = tmp

    marginal_dis = {}
    for t, logPseeding in logP_wdaughter.items():
        tmpP = defaultdict(float)
        for seeding_time in logPseeding:
            for n,v in zip(logPseeding[seeding_time][:,0], np.exp(logPseeding[seeding_time][:,1] - peak_val)):
                tmpP[n] += v
        marginal_dis[t] = np.array([(n,np.log(tmpP[n])+peak_val) for n in sorted(tmpP.keys())])

    plt.figure()
    for t in time_points:
        P = np.exp(marginal_dis[t][:,1] - marginal_dis[t][:,1].max())
        P/=P.sum()
        plt.plot(marginal_dis[t][:,0]+0.5, P, c=f"C{t%10}", lw=2, label=f'with c={parent_cases[t]}', ls='--')
        P = np.exp(logP[t][:,1])
        P/=P.sum()
        plt.plot(logP[t][:,0]+0.5, P, c=f"C{t%10}", lw=2, label=f'with c={daughter_cases[t]}')
    plt.xscale('log')
    plt.yscale('log')
    plt.legend()




