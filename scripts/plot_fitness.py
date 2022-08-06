import argparse,json
import pandas as pd
import numpy as np
from collections import defaultdict
import matplotlib.pyplot as plt
from scipy.stats import scoreatpercentile

if __name__=="__main__":
    parser = argparse.ArgumentParser(
        description="remove time info",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument('--fitness', type=str, required=True, help="input data")
    parser.add_argument('--mutations', type=str, required=True, help="input data")
    parser.add_argument('--output-fitness', type=str, help="fitness figure")
    parser.add_argument('--output-fitness-by-gene', type=str, help="fitness figure")
    parser.add_argument('--output-mutations', type=str, help="fitness figure")
    args = parser.parse_args()

    mutations = pd.read_csv(args.mutations, sep='\t', index_col=0)
    fitness = pd.read_csv(args.fitness, sep='\t', index_col=0)

    ## figure with mutation distributions
    plt.figure()
    sorted_muts = mutations.sort_values('scaled_rate', ascending=False)
    mut_sum = np.sum(sorted_muts.scaled_rate)
    muts = [(t,v/mut_sum) for t,v in sorted_muts.scaled_rate.items()]
    plt.bar(np.arange(len(muts)), height=[m[1] for m in muts])
    plt.xticks(np.arange(len(muts)), [m[0] for m in muts], rotation=30)
    plt.ylabel('fraction')
    plt.savefig(args.output_mutations)

    ## Figure with 1st, 2nd, 3rd positions
    plt.figure()
    for p in range(1,4):
        ind = fitness.pos_in_codon==p
        plt.plot(sorted(fitness.tolerance[ind]), np.linspace(0,1,ind.sum()),
                        label = f"codon pos={p}")

    ind = fitness.pos_in_codon.isna()
    plt.plot(sorted(fitness.tolerance[ind]), np.linspace(0,1,ind.sum()),
                    label = f"non-codon")

    ind = fitness.pos_in_codon==3
    syn_cutoff = scoreatpercentile(fitness.tolerance[ind],10)
    plt.plot([syn_cutoff, syn_cutoff], [0,1], c='k', alpha=0.3)
    for i in range(3):
        ind = fitness.pos_in_codon==i
        print(np.mean(fitness.tolerance[ind]<syn_cutoff))

    plt.xscale('log')
    plt.ylabel("fraction less")
    plt.xlabel("scaled number of lineages with mutations")
    plt.legend()
    plt.savefig(args.output_fitness)

    # # figure with pdfs instead of cdfs
    # # plt.figure()
    # measure = total_events_rescaled
    # bins = np.logspace(0, np.ceil(np.log10(measure.max())), 101)
    # bc = np.sqrt(bins[:-1]*bins[1:])
    # bins[0] = 0
    # rate_estimate = {}
    # for p in range(4):
    #     ind = codon_pos==p
    #     y,x = np.histogram(measure[ind], bins=bins)
    #     rate_estimate[p] = {"mean": np.sum(bc*y/y.sum()),
    #                         "geo-mean": np.exp(np.sum(np.log(bc)*y/y.sum())),
    #                         "median": np.median(measure[ind])}
    #     # plt.plot(bc,y/y.sum(), label = 'non-coding' if p==0 else f"codon pos={p}")

    # two panel figure with 1/2nd and 3rd position mutations
    fig, axs = plt.subplots(2,1, figsize = (6,10), sharex=True)
    pos = np.arange(len(fitness))
    genes = [x for x in fitness.gene.unique() if not pd.isna(x)]
    for i,gene in enumerate(genes):
        c = f"C{i}"
        ls = '--' if i>9 else '-'
        ind = (fitness.pos_in_codon==3) & (fitness.gene==gene)
        axs[1].plot(sorted(fitness.tolerance[ind]), np.linspace(0,1,ind.sum()), ls=ls, c=c,label = f'{gene}')

        ind = (fitness.pos_in_codon>0) & (fitness.pos_in_codon<3) & (fitness.gene==gene)
        axs[0].plot(sorted(fitness.tolerance[ind]), np.linspace(0,1,ind.sum()),
                        ls=ls, c=c)

    ind = fitness.pos_in_codon.isna()
    axs[1].plot(sorted(fitness.tolerance[ind]), np.linspace(0,1,ind.sum()),
                    label = 'non-coding', c='k')

    plt.xscale('log')
    plt.xlim(0.01, 10)
    for ax in axs:
        ax.grid()

    axs[1].legend(ncol=2)
    axs[1].set_title("3rd codon positions or non-coding")
    axs[0].set_title("1st and 2nd codon positions")
    axs[0].set_ylabel("fraction below")
    axs[1].set_ylabel("fraction below")
    axs[1].set_xlabel("scaled number of lineages with mutations")
    plt.savefig(args.output_fitness_by_gene)

    # fitness_cost = np.array(fitness_cost)
    # plt.figure()
    # pos = np.arange(len(ref_array))
    # ws = 10
    # w = np.ones(ws)/ws
    # for i,gene in enumerate(gene_position):
    #     if gene=='ORF9b': continue

    #     gene_range = gene_position[gene]
    #     plt.figure()
    #     for cp in range(3):
    #         #ind = (codon_pos==(cp+1)) & (pos>=gene_range.start) & (pos<gene_range.end) & (~np.isnan(fitness_cost).all(axis=1))
    #         #plt.hist(np.min(np.log(fitness_cost[ind]), axis=1))
    #         ind = (codon_pos==(cp+1)) & (pos>=gene_range.start) & (pos<gene_range.end)
    #         gene_pos = (np.convolve(pos[ind],w, mode='valid') - gene_range.start)/3
    #         plt.plot(gene_pos,
    #                  np.convolve(np.log10(measure[ind]+.03), w, mode='valid'),
    #                             c=f'C{cp}', label=f'codon pos {cp+1}')
    #         plt.plot(gene_pos, np.zeros_like(gene_pos), c='k', alpha=0.3, lw=2)
    #         plt.plot((pos[ind]-gene_range.start)/3,
    #                  np.log10((np.log(number_of_pango_muts[ind]+1)/100+0.01)), 'o', c=f'C{cp}')
    #         plt.ylabel('log10 scaled mutations')
    #     plt.ylim(-1.5, 1)
    #     plt.title(gene)
