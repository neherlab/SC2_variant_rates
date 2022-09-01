import argparse
import numpy as np
import pandas as pd
from scipy.stats import linregress
import matplotlib as mpl
mpl.rcParams['axes.formatter.useoffset'] = False
from root_to_tip import add_panel_label

import matplotlib.pyplot as plt
import seaborn as sns


if __name__=="__main__":
    parser = argparse.ArgumentParser(
        description="remove time info",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    fs=14
    parser.add_argument('--rate-table', type=str, required=True, help="input data")
    parser.add_argument('--output-plot', type=str, help="plot file")
    parser.add_argument('--output-plot-rates', type=str, help="plot file")
    parser.add_argument('--output-plot-rates-genes', type=str, help="plot file")
    args = parser.parse_args()

    rates = pd.read_csv(args.rate_table, sep='\t', index_col='clade')
    inter_clade_rates = {}
    fig, axs = plt.subplots(2,2, figsize=(12,12))
    xticks = [2020,2020.5,2021,2021.5,2022,2022.5]
    for ax, mut_type, ax_label, panel in zip(axs.flatten(),['nuc', 'aa', 'syn'],
                        ['total divergence', 'amino acid divergence', 'synonymous divergence'], ['A', 'B', 'C']):
        ax.set_ylabel(ax_label, fontsize=fs)
        add_panel_label(ax, panel, fs=fs*1.8)
        inter_clade = []
        ci=0
        for clade, row in rates.iterrows():
            dt = np.linspace(0, 0.75,2)
            t = dt + row[f'{mut_type}_origin']
            slope = row[f'{mut_type}_rate']
            ls = '-' if ci<10 else '--'
            m = 'o' if ci<10 else 's'
            ax.plot(t, dt*slope + row[f'{mut_type}_div'], label=clade, c=f"C{ci%10}", ls=ls)
            ax.scatter([[row[f'{mut_type}_origin']]], [[row[f'{mut_type}_div']]], c=f"C{ci%10}", marker=m)
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
        ax.text( 0.4, 0.06,f"overall rate: {reg.slope:1.1f}/year", fontsize=fs, transform=ax.transAxes)
        if mut_type == 'aa': ax.legend(ncol=2)
        ax.set_xticks(xticks, [str(x) for x in xticks], fontsize=fs)

    for i, (mut_type, rate) in enumerate(inter_clade_rates.items()):
        axs[-1,-1].fill_between([i-0.5, i+0.5], [35,35], facecolor='k', alpha=0.075*(2+i%2))
        axs[-1,-1].plot([i-0.4, i+0.4], [rate, rate], lw=3, c='k', alpha=0.5)
        clade_rates = rates[f"{mut_type}_rate"]
        x_offset = i + np.linspace(-.35, 0.35,len(clade_rates))
        axs[-1,-1].scatter(x_offset[:10],
                           clade_rates[:10], c=[f"C{ci%10}" for ci in range(10)],
                           marker='o')
        axs[-1,-1].scatter(x_offset[10:],
                           clade_rates[10:], c=[f"C{ci%10}" for ci in range(len(clade_rates) - 10)],
                           marker='s')
    axs[-1,-1].set_ylabel("substitutions per year", fontsize=fs)
    axs[-1,-1].set_xticks([0,1,2], ['total', 'amino acid', 'synonymous'], rotation=20,
                          ha='center', fontsize=fs)
    axs[-1,-1].set_ylim(0)
    add_panel_label(axs[-1,-1], 'D', fs=fs*1.8)

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
    plt.ylim(0)
    plt.legend()

    if args.output_plot_rates:
        plt.savefig(args.output_plot_rates)
    else:
        plt.show()

    plt.figure()
    plt.plot(rates["aa_rate"], 'o-', label='Overall amino acid rate')
    plt.plot(rates["spike_rate"], 's-', label='spike protein')
    plt.plot(rates["orf1ab_rate"], 'd-', label='ORF1ab')
    plt.plot(rates["enm_rate"], 'd-', label='E,M,N')
    plt.plot(rates["aa_rate"] - rates["spike_rate"] - rates["orf1ab_rate"]- rates["enm_rate"],
             'v-', label='other ORFs')
    plt.ylabel('rate estimate [subs/y]')
    plt.legend()
    plt.xticks(range(len(rates)), rates.index, rotation=60, ha='right')
    plt.tight_layout()

    if args.output_plot_rates_genes:
        plt.savefig(args.output_plot_rates_genes)
    else:
        plt.show()
