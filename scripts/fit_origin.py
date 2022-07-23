import argparse,json
import numpy as np
from datetime import datetime
import matplotlib as mpl
from fit_outbreak import optimize_single_trajectory
mpl.rcParams['axes.formatter.useoffset'] = False

import matplotlib.pyplot as plt
import seaborn as sns


if __name__=="__main__":
    parser = argparse.ArgumentParser(
        description="remove time info",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument('--counts', type=str, required=True, help="input data")
    parser.add_argument('--output-plot', type=str, help="figure file")
    args = parser.parse_args()

    with open(args.counts) as fh:
        counts = json.load(fh)

    obs_dates = [datetime.fromordinal(x) for x in counts['bins']]

    fig = plt.figure()
    cutoff = 6
    eps=0.02
    print("fit all samples")
    fit = optimize_single_trajectory(counts['bins'][:cutoff], counts['all_samples'][:cutoff], sampling=eps)
    dates = [datetime.fromordinal(x) for x in fit['timepoints']]
    peak = max(fit['logLH'])
    plt.plot(dates, np.exp(fit['logLH'] - peak))

    print("fit founder")
    fit2 = optimize_single_trajectory(counts['bins'][:cutoff], counts['mutations']['0'][:cutoff], sampling=eps)
    dates = [datetime.fromordinal(x) for x in fit2['timepoints']]
    peak = max(fit2['logLH'])
    plt.plot(dates, np.exp(fit2['logLH'] - peak))

    print("fit one mutation class")
    fit3 = optimize_single_trajectory(counts['bins'][:cutoff], counts['mutations']['1'][:cutoff], sampling=eps)
    dates = [datetime.fromordinal(x) for x in fit3['timepoints']]
    peak = max(fit3['logLH'])
    plt.plot(dates, np.exp(fit3['logLH'] - peak))


    fig.autofmt_xdate()
