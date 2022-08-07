import argparse,glob
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

def get_freq(d, cov_cutoff=1000):
    cov = d.Good_depth
    freq = np.array([d[f"Count_{n}"]/cov for n in 'ACGT']).T
    freq[cov<cov_cutoff] = np.nan
    return freq

def get_diversity(freq):
    return 1-np.sum(freq**2, axis=1)

if __name__=="__main__":
    parser = argparse.ArgumentParser(
        description="plot within host evolution",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument('--allele-counts-dir', type=str, required=True, help="input data")
    parser.add_argument('--output-intra-host', type=str, help="fitness figure")
    args = parser.parse_args()

    fnames = glob.glob(args.allele_counts_dir + "/*")
    freqs = []
    div = []
    for fname in fnames:
        d = pd.read_csv(fname, sep='\t')
        freqs.append(get_freq(d))
        div.append(get_diversity(freqs[-1]))

    masked_diversity = np.ma.array(div, mask=np.isnan(div))


    median_div = np.ma.median(masked_diversity, axis=0)

    smoothed_diversity = np.exp(np.ma.convolve(np.log(median_div.filled(fill_value=1e-5)+1e-5), np.ones(30)/30, mode='same'))

    minor_variants = [[] for i in range(4)]
    for f in freqs:
        cons = np.ma.argmax(f, axis=1)
        valid = ~np.isnan(np.ma.max(f, axis=1))
        for i in range(4):
            ind = (cons==i)&valid
            minor_variants[i].append(np.median(f[ind,:], axis=0))

    M = np.array([np.mean(minor_variants[i], axis=0) for i in range(4)])
    np.fill_diagonal(M,0)

