import argparse,json
import pandas as pd
import numpy as np
import matplotlib as mpl
from Bio import SeqIO,Seq
mpl.rcParams['axes.formatter.useoffset'] = False
from collections import defaultdict
from scipy.stats import scoreatpercentile

import matplotlib.pyplot as plt
import seaborn as sns

columns = ['seqName', 'Nextclade_pango', "privateNucMutations.unlabeledSubstitutions", 'qc.overallStatus']

def get_sequence(mods, root_seq):
    seq = root_seq.copy()
    for mut in mods['nuc']:
        a,pos,d = mut[0], int(mut[1:-1])-1, mut[-1]
        if a!=seq[pos]:
            print(seq[pos], mut)
        seq[pos]=d
    return seq

def translate(codon):
    try:
        return Seq.translate("".join(codon))
    except:
        return 'X'

if __name__=="__main__":
    parser = argparse.ArgumentParser(
        description="remove time info",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument('--metadata', type=str, required=True, help="input data")
    parser.add_argument('--reference', type=str, required=True, help="input data")
    parser.add_argument('--pango-gts', type=str, required=True, help="input data")
    parser.add_argument('--output-fitness', type=str, required=True, help="fitness figure")
    parser.add_argument('--output-fitness-by-gene', type=str, required=True, help="fitness figure")
    parser.add_argument('--output-mutations', type=str, required=True, help="fitness figure")
    args = parser.parse_args()

    ## load reference
    ref = SeqIO.read(args.reference, 'genbank')
    ref_array = np.array(ref.seq)
    base_content = {x:ref.seq.count(x) for x in 'ACGT'}

    ## make a map of codon positions, ignore orf9b (overlaps N)
    codon_pos = np.zeros(len(ref))
    map_to_gene = {}
    gene_position = {}
    gene_length = {}
    for feat in ref.features:
        if feat.type=='CDS':
            gene_name = feat.qualifiers['gene'][0]
            gene_position[gene_name] = feat.location
            gene_length[gene_name] = (feat.location.end - feat.location.start)//3
            if gene_name!='ORF9b':
                for gpos, pos in enumerate(feat.location):
                    codon_pos[pos] = (gpos%3)+1
                    map_to_gene[pos]= gene_name


    # load pango genotypes
    with open(args.pango_gts) as fh:
        pango_gts = json.load(fh)

    # load and filter metadata
    d = pd.read_csv(args.metadata, sep='\t', usecols=columns).fillna('')
    d = d.loc[d["qc.overallStatus"]=='good',:]

    # glob all rare mutations by pango lineage
    mutation_counter = defaultdict(lambda: defaultdict(int))
    # pango_counter = defaultdict(int)
    for r, row in d.iterrows():
        pango = row.Nextclade_pango
        # pango_counter[pango] += 1
        if pango[0]=='X':
            continue

        muts = row["privateNucMutations.unlabeledSubstitutions"]
        if muts:
            for m in muts.split(','):
                mutation_counter[pango][m] += 1


    # for each pango lineage, count mutations by position, as well as synoymous and non-synonyomous by codon
    position_counter = defaultdict(lambda: defaultdict(list))
    syn =    defaultdict(lambda: defaultdict(lambda: defaultdict(list)))
    nonsyn = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))

    for pango, pmuts in mutation_counter.items():
        if pango not in pango_gts:
            print("missing pango", pango)
            continue
        seq = get_sequence(pango_gts[pango], ref_array)
        # translate each codon in pango reference sequence
        for mut,c in pmuts.items():
            a, pos, d = mut[0], int(mut[1:-1])-1, mut[-1]
            if a=='-':
                continue
            position_counter[pos][(a,d)].append(c)
            if pos in map_to_gene:
                gene = map_to_gene[pos]
                pos_in_gene = pos - gene_position[gene].start
                frame = pos_in_gene%3
                codon = seq[pos-frame:pos-frame+3]
                codon_number = pos_in_gene//3
                alt = codon.copy()
                alt[frame] = d
                aa = translate(codon)
                alt_aa = translate(alt)
                if aa==alt_aa and aa!='X':
                    syn[gene][codon_number][a].append(c)
                elif 'X' not in [aa,alt_aa]:
                    nonsyn[gene][codon_number][a].append(c)

    # determine the total number of events and counts at each position
    def total_len_filtered_lists(l, cutoff):
        return sum([len([y for y in x if y>cutoff]) for x in l])

    cutoff=1
    total_events = np.array([total_len_filtered_lists(position_counter[pos].values(), cutoff)
                            if pos in position_counter else 0
                            for pos in range(len(ref))])

    total_counts = np.array([sum([sum(x) for x in position_counter[pos].values()])
                             if pos in position_counter else 0
                             for pos in range(len(ref))])

    # sort all mutations by type of mutation
    events_by_transition = defaultdict(int)
    for m in position_counter.values():
        for t in m:
            events_by_transition[t]+=len(m[t])
    # determine the mutations rate for each type, as well as the rate out of a nuc
    transition_rate = dict()
    away_rate = defaultdict(float)
    for a,d in events_by_transition:
        if a in base_content:
            transition_rate[(a,d)] = events_by_transition[(a,d)]/base_content[a]
            away_rate[a] += transition_rate[(a,d)]

    # average mutation rate used to scale rates
    avg_rate = np.sum([base_content[n]*away_rate[n]/len(ref) for n in 'ACGT'])
    total_events_rescaled = np.array([c/away_rate[nuc]*avg_rate for c,nuc in zip(total_events, ref.seq)])

    # calculate synonymous and non-synonymous distributions
    total_events_by_type = {}
    for mut_counts, label in [(syn, 'syn'), (nonsyn, 'nonsyn')]:
        total_events_by_type[label] = {}
        for gene in mut_counts:
            total_events_by_type[label][gene] = [np.sum([len([x for x in mut_counts[gene][pos][a] if x>cutoff])/away_rate[a]*avg_rate
                                                     for a in mut_counts[gene][pos]])
                                                   for pos in range(gene_length[gene])]

    ## Figure with 1st, 2nd, 3rd positions
    measure = total_events_rescaled
    plt.figure()
    for p in range(4):
        ind = codon_pos==p
        plt.plot(sorted(measure[ind]), np.linspace(0,1,ind.sum()),
                        label = 'non-coding' if p==0 else f"codon pos={p}")

    ind = codon_pos==3
    syn_cutoff = scoreatpercentile(measure[ind],10)
    plt.plot([syn_cutoff, syn_cutoff], [0,1], c='k', alpha=0.3)
    for i in range(3):
        ind = codon_pos==i
        print(np.mean(measure[ind]<syn_cutoff))

    plt.xscale('log')
    plt.ylabel("fraction less")
    plt.xlabel("scaled number of lineages with mutations")
    plt.legend()
    plt.savefig(args.output_fitness)

    # figure with pdfs instead of cdfs
    # measure = total_events_rescaled
    # plt.figure()
    # bins = np.logspace(0, np.ceil(np.log10(measure.max())), 31)
    # bc = np.sqrt(bins[:-1]*bins[1:])
    # bins[0] = 0
    # for p in range(4):
    #     ind = codon_pos==p
    #     y,x = np.histogram(measure[ind], bins=bins)
    #     plt.plot(bc,y/y.sum(), label = 'non-coding' if p==0 else f"codon pos={p}")

    # plt.xscale('log')
    # plt.ylabel("distribution")
    # plt.xlabel("scaled number of lineages with mutations")
    # plt.legend()

    ## figure with mutation distributions
    plt.figure()
    muts = sorted(transition_rate.items(), key=lambda x:x[1], reverse=True)
    mut_sum = np.sum([x[1] for x in muts])
    muts = [(t,v/mut_sum) for t,v in muts]
    plt.bar(np.arange(len(muts)), height=[m[1] for m in muts])
    plt.xticks(np.arange(len(muts)), [f"{m[0][0]}->{m[0][1]}" for m in muts], rotation=30)
    plt.ylabel('fraction')
    plt.savefig(args.output_mutations)

    # two panel figure with 1/2nd and 3rd position mutations
    measure = total_events_rescaled
    fig, axs = plt.subplots(2,1, figsize = (6,10), sharex=True)
    pos = np.arange(len(ref_array))
    for i,gene in enumerate(gene_position):
        c = f"C{i}"
        ls = '--' if i>9 else '-'
        gene_range = gene_position[gene]
        ind = (codon_pos==3) & (pos>=gene_range.start) & (pos<gene_range.end)
        axs[1].plot(sorted(measure[ind]), np.linspace(0,1,ind.sum()), ls=ls, c=c,label = f'{gene}')

        ind =  (codon_pos<3) & (pos>=gene_range.start) & (pos<gene_range.end)
        axs[0].plot(sorted(measure[ind]), np.linspace(0,1,ind.sum()),
                        ls=ls, c=c)

    ind = codon_pos==0
    axs[1].plot(sorted(measure[ind]), np.linspace(0,1,ind.sum()),
                    label = 'non-coding', c='k')

    axs[1].legend(ncol=2)
    axs[1].set_title("3rd codon positions or non-coding")
    axs[0].set_title("1st and 2nd codon positions")
    plt.xscale('log')
    axs[0].set_ylabel("fraction below")
    axs[1].set_ylabel("fraction below")
    axs[1].set_xlabel("scaled number of lineages with mutations")
    plt.savefig(args.output_fitness_by_gene)

    # ## strange pattern in E.
    # measure = total_events_rescaled
    # plt.figure()
    # gene='E'
    # gene_range = gene_position[gene]
    # for i in [1,2,3]:
    #     ind =(codon_pos==i) & (pos>=gene_range.start) & (pos<gene_range.end)
    #     plt.plot(pos[ind],measure[ind], 'o', c=f"C{i}")
