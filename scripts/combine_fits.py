import argparse
from turtle import fillcolor
import numpy as np
import pandas as pd
from scipy.stats import linregress
import matplotlib as mpl
mpl.rcParams['axes.formatter.useoffset'] = False

import matplotlib.pyplot as plt
import seaborn as sns


if __name__=="__main__":
    parser = argparse.ArgumentParser(
        description="remove time info",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument('--rate-table', type=str, required=True, help="input data")
    parser.add_argument('--output-plot', type=str, help="plot file")
    parser.add_argument('--output-plot-rates', type=str, help="plot file")
    args = parser.parse_args()

    rates = pd.read_csv(args.rate_table, sep='\t', index_col='clade')
    inter_clade_rates = {}
    fig, axs = plt.subplots(2,2, figsize=(12,12))
    for ax, mut_type, ax_label in zip(axs.flatten(),['nuc', 'aa', 'syn'],['total divergence', 'aa divergence', 'syn divergence']):
        ax.set_ylabel(ax_label)
        inter_clade = []
        ci=0
        for clade, row in rates.iterrows():
            dt = np.linspace(0, 0.75,2)
            t = dt + row[f'{mut_type}_origin']
            slope = row[f'{mut_type}_rate']
            ax.plot(t, dt*slope + row[f'{mut_type}_div'], label=clade, c=f"C{ci%10}")
            ax.scatter([[row[f'{mut_type}_origin']]], [[row[f'{mut_type}_div']]], c=f"C{ci%10}")
            if row[f'{mut_type}_origin']>2019.7 and row[f'{mut_type}_origin']<2022.7:
                inter_clade.append([row[f'{mut_type}_origin'], row[f'{mut_type}_div']])
            ci += 1
        inter_clade = np.array(inter_clade)
        reg = linregress(inter_clade[:,0], inter_clade[:,1])
        inter_clade_rates[mut_type] = reg.slope

        x = np.linspace(2019.8, 2022.5, 21)
        y = reg.slope*x + reg.intercept
        std_dev = np.sqrt(np.maximum(0,reg.slope*x + reg.intercept))
        ax.plot(x, y, c='k', lw=3, alpha=0.5)
        ax.fill_between(x, y+std_dev, np.maximum(0, y-std_dev), fc='k', lw=3, alpha=0.1, ec=None)
        ax.set_xlim(x.min(), x.max())
        ax.set_ylim(0)
        if mut_type == 'aa': ax.legend(ncol=2)
    for i, (mut_type, rate) in enumerate(inter_clade_rates.items()):
        axs[-1,-1].plot([i-0.4, i+0.4], [rate, rate], lw=3, c='k', alpha=0.5)
        clade_rates = rates[f"{mut_type}_rate"]
        axs[-1,-1].scatter(i - 0.35 + np.random.random(size=len(clade_rates))*0.7,
                           clade_rates, c=[f"C{ci%10}" for ci in range(len(clade_rates))])
    axs[-1,-1].set_ylabel("substitutions per year")
    axs[-1,-1].set_xticks([0,1,2], ['nuc', 'aa', 'syn'])
    axs[-1,-1].set_ylim(0)

    if args.output_plot:
        plt.savefig(args.output_plot)
    else:
        plt.show()

    plt.figure()
    plt.plot(rates["nuc_origin"], rates["aa_rate"], 'o', label='amino acid rate')
    plt.plot(rates["nuc_origin"], rates["syn_rate"], 'o', label='synonymous rate')
    plt.plot(rates["nuc_origin"], np.ones_like(rates['nuc_origin'])*inter_clade_rates["aa"], label='inter-clade amino acid rate', c=f"C{0}", lw=3)
    plt.plot(rates["nuc_origin"], np.ones_like(rates['nuc_origin'])*inter_clade_rates["syn"], label='inter-clade synonymous rate', c=f"C{1}", lw=3)
    plt.ylabel('rate estimate [1/y]')
    plt.legend()

    if args.output_plot_rates:
        plt.savefig(args.output_plot_rates)
    else:
        plt.show()
