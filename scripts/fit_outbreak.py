import numpy as np
from scipy.stats import nbinom, poisson
from scipy.optimize import minimize
import matplotlib.pyplot as plt
from collections import defaultdict

TINY = 1e-100
LogTINY = np.log(TINY)
birth_rate = 1e-6
logCutoff = 50

def log_propagator(n, m, rho, k):
    '''
    nvec = np.arange(100)
    k,rho = 0.1, 1.1
    pn = propagator(nvec, 5, rho,k)
    avg = np.sum(pn*nvec) # should be rho*5
    '''
    p = k/(k+rho)
    return nbinom.logpmf(n,m*k,p)

def propagator(n, m, rho, k):
    return np.exp(log_propagator(n,m,rho,k))


def log_sampling_prob(n_samples, eps, n_cases):
    mean_density = eps*n_cases
    return poisson.logpmf(n_samples, mean_density)


def back_propagate_R(sampling_t, logR_tp1, rho, k, nmin=1, nmax=1000):
    nvals = np.arange(max(1,int(nmin)), int(nmax))
    if logR_tp1 is None:
        logR_t = np.array([nvals, np.sum([log_sampling_prob(s_t[0], s_t[1], nvals) for s_t in sampling_t], axis=0)])
    else:
        peak = logR_tp1[1].max()
        sampling =  np.sum([log_sampling_prob(s_t[0], s_t[1], nvals) for s_t in sampling_t], axis=0)
        tmp_R = np.array([np.log(np.sum(np.exp(logR_tp1[1] - peak + log_propagator(logR_tp1[0], n_t, rho, k))))
                       for n_t in nvals])
        logR_t = np.array([nvals, tmp_R + sampling + peak])

    peak = logR_t[1].max()

    return logR_t[:, logR_t[1]>peak - logCutoff]

def fwd_propagate_Q(sampling_tm1, logQ_tm1, rho, k, nmin=1, nmax=1000):
    if logQ_tm1 is None:
        logQ_t = np.array([[1], [0]])
    else:
        nvals = np.arange(max(1,int(nmin)), int(nmax))
        peak = logQ_tm1[1].max()
        sampling =  np.sum([log_sampling_prob(s_t[0], s_t[1], nvals) for s_t in sampling_tm1], axis=0)
        tmpQ = np.array([np.log(np.sum(np.exp(sampling_n \
                                       + logQ_tm1[1] - peak + log_propagator(n_tm1, logQ_tm1[0], rho, k))))
                           for sampling_n, n_tm1 in zip(sampling, nvals)])
        logQ_t = np.array([nvals, tmpQ + peak])

    peak = logQ_t[1].max()
    return logQ_t[:, logQ_t[1]>peak - logCutoff]


def back_trace(data):
    sorted_data = data.sort_index(ascending=False)
    logR = {}
    for i, (t, row) in enumerate(sorted_data.iterrows()):
        prev = logR[sorted_data.iloc[i-1].name] if i else None
        logR[t] = back_propagate_R([(row.cases, row.eps)], prev, row.rho, row.k, nmin=row.cases, nmax=3*(row.cases+3)/eps)

    return logR


def fwd_trace(data, tau):
    sorted_data = data.sort_index(ascending=True)
    sorted_data = sorted_data.loc[sorted_data.index>=tau]
    logQ = {}
    for i, (t, row) in enumerate(sorted_data.iterrows()):
        prev = logQ[sorted_data.iloc[i-1].name] if i else None
        prev_cases = sorted_data.iloc[i-1].cases if i else None
        logQ[t] = fwd_propagate_Q([(prev_cases, row.eps)], prev, row.rho, row.k, nmin=row.cases, nmax=3*(row.cases+3)/eps)

    return logQ

def negLogLH(data):
    logR = back_trace(data)
    peak = np.max([x[1,0] for x in logR.values() if x[0,0]==1])
    res = peak + np.log(np.sum([np.exp(x[1,0] - peak) for x in logR.values() if x[0,0]==1]))
    return -res

def cost(x, data, fields):
    print(x)
    for y, n in zip(x,fields):
        data[n] = y**2
    #data['eps'] = x[0]
    return negLogLH(data)


def marginal_distribution(logQ, logR):
    nmin, nmax = int(max(logQ[0,0],logR[0,0])), int(min(logQ[-1,0], logR[-1,0]))
    indR = (logR[:,0]>=nmin) & (logR[:,0]<nmax)
    indQ = (logQ[:,0]>=nmin) & (logQ[:,0]<nmax)
    tmplogP = logR[1,indR] + logQ[1, indQ]
    return np.array([logR[0,indR], tmplogP])


def optimize_single_trajectory(bins, cases, sampling):
    bin_width = bins[1]-bins[0]
    max_seeding = np.min([ti for ti, x in enumerate(cases) if x>0])
    if max_seeding<10:
        cases = [0]*(10-max_seeding) + list(cases)
        min_bin = bins[0]
        bins = [min_bin - bin_width*i for i in range(10-max_seeding, 0,-1)] + list(bins)

    rho = 1.3
    k=.2
    data = pd.DataFrame({i:{'cases':parent_cases[i], 'rho':rho, 'eps':sampling,
                         'k':k, 'mu':birth_rate} for i in range(len(cases))}).T
    sol = minimize(cost, [rho], (data, ['rho']))
    logR = back_trace(data)
    return bins, [np.exp(logR[t][1,0]) if logR[t][0,0]==1 else -np.inf for t in bins]



if __name__=="__main__":
    import pandas as pd
    from scipy.optimize import minimize

    parent_cases =   [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 5, 10, 40, 70, 100]
    rho = 1.3
    k=.2
    eps=0.02
    data = pd.DataFrame({i:{'cases':parent_cases[i], 'rho':rho, 'eps':eps, 'k':k, 'mu':birth_rate} for i in range(len(parent_cases))}).T
    sol = minimize(cost, [rho], (data, ['rho']))
    rho = sol['x'][0]

    logR = back_trace(data)

    daughter_cases = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 6, 20, 30, 50]
    data_daughter = pd.DataFrame({i:{'cases':daughter_cases[i], 'rho':rho, 'eps':eps, 'k':k, 'mu':TINY} for i in range(len(parent_cases))}).T

    max_seeding_daughter = np.min([ti for ti, x in enumerate(daughter_cases) if x>0])
    if daughter_cases[max_seeding_daughter]>1:
        max_seeding_daughter-=1
    max_seeding_parent = min(max_seeding_daughter, np.min([ti for ti, x in enumerate(parent_cases) if x>0]))
    if parent_cases[max_seeding_parent]>1:
        max_seeding_parent-=1
    logQ = {}
    for tau in range(0,max_seeding_parent+1):
        logQ[tau] = fwd_trace(data, tau)

    logR_daughter = back_trace(data_daughter)

    logP = np.ones((max_seeding_parent+1, max_seeding_daughter+1))*np.nan
    mutation_rate = 0.1
    for tau_parent in logQ:
        for tau_daughter in range(tau_parent,max_seeding_daughter+1):
            nmin = max(logR[tau_daughter][0,0], logQ[tau_parent][tau_daughter][0,0])
            nmax = min(logR[tau_daughter][0,-1], logQ[tau_parent][tau_daughter][0,-1])
            nvals = np.arange(nmin, nmax+1, dtype=int)
            Qind = (logQ[tau_parent][tau_daughter][0]>=nmin)&(logQ[tau_parent][tau_daughter][0]<=nmax)
            Rind = (logR[tau_daughter][0]>=nmin)&(logR[tau_daughter][0]<=nmax)
            res = logR_daughter[tau_daughter][1,0] if logR_daughter[tau_daughter][0,0]==1 else np.nan
            peak = (logQ[tau_parent][tau_daughter][1,Qind] + logR[tau_daughter][1,Rind]).max()
            res += peak + np.log(np.sum([np.exp(r + q - peak)*(1-np.exp(-n*mutation_rate))
                          for n,r,q in zip(logR[tau_daughter][0,Rind], logR[tau_daughter][1,Rind], logQ[tau_parent][tau_daughter][1,Qind])]))

            logP[tau_parent,tau_daughter] = res

    ML = logP[~np.isnan(logP)].max()
    P = np.exp(logP-ML)
    P[np.isnan(P)] = 0

    plt.figure()
    plt.plot(P.sum(axis=1), label='parent')
    plt.plot(P.sum(axis=0), label='daugher')
    plt.xlabel('time')
    plt.legend()
