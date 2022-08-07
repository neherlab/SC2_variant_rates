import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

if __name__=="__main__":
    parser = argparse.ArgumentParser(
        description="remove time info",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument('--fitness', type=str, required=True, help="input data")
    parser.add_argument('--output-fitness-landscape', type=str, help="fitness figure")
    args = parser.parse_args()

    fitness = pd.read_csv(args.fitness, sep='\t', index_col=0)
    genes = [x for x in fitness.gene.unique() if not pd.isna(x)]
    pos = np.arange(len(fitness))

    fig = plt.figure(figsize=(15,12))
    gene_groups = [ ['ORF3a', 'E', 'M', 'ORF6','ORF7a', 'ORF7b', 'ORF8'],['S', 'N'],['ORF1b'],['ORF1a'] ]
    gene_length = {k:(fitness.gene==k).sum() for k in genes}

    fs = 14
    ws = 20
    w = np.ones(ws)/ws
    ws_fine = 7
    w_fine = np.ones(ws_fine)/ws_fine
    lower_fitness_cutoff = 0.03
    n_rows = len(gene_groups)
    h_spread = 0.003
    h_margin = 0.05
    v_spread = 0.02
    v_margin = 0.05
    height = (1 - v_spread*n_rows - v_margin)/n_rows
    for row,group in enumerate(gene_groups):
        axis_bottom = row*(v_spread + height) + v_margin
        available_width = 1 - len(group)*h_spread - h_margin
        left = h_margin
        total_genome = np.sum([gene_length[g] for g in group])
        for gi,gene in enumerate(group):
            width = gene_length[gene]/total_genome*available_width
            ax = fig.add_axes((left,axis_bottom, width, height))
            left += width + h_spread
            for cp in range(3):
                ind = (fitness.pos_in_codon==(cp+1)) & (fitness.gene==gene)
                gene_pos = (fitness.codon[ind]*3 + cp - 1)/3
                gene_pos_smooth = np.convolve(gene_pos, w, mode='valid')
                ax.plot(gene_pos_smooth,
                        np.convolve(np.log10(fitness.tolerance[ind]+lower_fitness_cutoff), w, mode='valid'),
                                    c=f'C{cp}', label=f'codon pos {cp+1}')

                gene_pos_smooth = np.convolve(gene_pos, w_fine, mode='valid')
                ax.plot(gene_pos_smooth,
                        np.convolve(np.log10(fitness.tolerance[ind]+lower_fitness_cutoff), w_fine, mode='valid'),
                                    c=f'C{cp}', alpha = 0.5)

                ax.plot(gene_pos,
                    np.log10(fitness.lineage_fraction_with_changes[ind]+.001)/3 - 1.2, 'o', c=f'C{cp}')

            ax.plot(gene_pos, np.zeros_like(gene_pos), c='k', alpha=0.3, lw=2)
            ax.set_ylim(-2.0, 1)
            ax.set_xlim(0, len(gene_pos))

            if gi==0:
                ax.set_ylabel('log10 scaled tolerance')
                if row==n_rows-1:
                    ax.legend(ncol=3, fontsize=fs)
            else:
                ax.set_yticklabels([])

            ax.text(0.15,0.75, gene, fontsize=fs*1.1)

    plt.savefig(args.output_fitness_landscape)
