import argparse,json
import pandas as pd
import numpy as np
import matplotlib as mpl
from Bio import SeqIO
mpl.rcParams['axes.formatter.useoffset'] = False
from collections import defaultdict

import matplotlib.pyplot as plt
import seaborn as sns

columns = ['seqName', 'Nextclade_pango', "privateNucMutations.unlabeledSubstitutions", 'qc.overallStatus']

if __name__=="__main__":
    parser = argparse.ArgumentParser(
        description="remove time info",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument('--metadata', type=str, required=True, help="input data")
    parser.add_argument('--reference', type=str, required=True, help="input data")
    parser.add_argument('--output-fitness', type=str, required=True, help="fitness figure")
    parser.add_argument('--output-mutations', type=str, required=True, help="fitness figure")
    args = parser.parse_args()

    ref = SeqIO.read(args.reference, 'genbank')

    d = pd.read_csv(args.metadata, sep='\t', usecols=columns).fillna('')
    d = d.loc[d["qc.overallStatus"]=='good',:]

    mutation_counter = defaultdict(lambda: defaultdict(int))

    for r, row in d.iterrows():
        pango = row.Nextclade_pango
        if pango[0]=='X':
            continue

        muts = row["privateNucMutations.unlabeledSubstitutions"]
        if muts:
            for m in muts.split(','):
                mutation_counter[pango][m] += 1


    position_counter = defaultdict(lambda: defaultdict(list))

    for pmuts in mutation_counter.values():
        for mut,c in pmuts.items():
            a, pos, d = mut[0], int(mut[1:-1])-1, mut[-1]
            position_counter[pos][(a,d)].append(c)


    cutoff=1
    total_events = np.array([sum([len([y for y in x if y>cutoff]) for x in position_counter[pos].values()])
                        if pos in position_counter else 0 for pos in range(len(ref))])

    total_counts = np.array([sum([sum(x) for x in position_counter[pos].values()])
                        if pos in position_counter else 0 for pos in range(len(ref))])

    events_by_transition = defaultdict(int)
    for m in position_counter.values():
        for t in m:
            events_by_transition[t]+=len(m[t])

    base_content = {x:ref.seq.count(x) for x in 'ACGT'}
    transition_rate = dict()
    away_rate = defaultdict(float)
    for a,d in events_by_transition:
        if a in base_content:
            transition_rate[(a,d)] = events_by_transition[(a,d)]/base_content[a]
            away_rate[a] += transition_rate[(a,d)]

    avg_rate = np.sum([base_content[n]*away_rate[n]/len(ref) for n in 'ACGT'])
    total_events_rescaled = np.array([c/away_rate[nuc]*avg_rate for c,nuc in zip(total_events, ref.seq)])

    codon_pos = np.zeros(len(ref))
    for feat in ref.features:
        if feat.type=='CDS' and feat.qualifiers['gene'][0]!='ORF9b':
            for gpos, pos in enumerate(feat.location):
                codon_pos[pos] = (gpos%3)+1


    measure = total_events_rescaled
    plt.figure()
    for p in range(4):
        ind = codon_pos==p
        plt.plot(sorted(measure[ind]), np.linspace(0,1,ind.sum()),
                        label = 'non-coding' if p==0 else f"codon pos={p}")

    plt.xscale('log')
    plt.ylabel("fraction less")
    plt.xlabel("scaled dnumber of lineages with mutations")
    plt.legend()
    plt.savefig(args.output_fitness)


    plt.figure()
    muts = sorted(transition_rate.items(), key=lambda x:x[1], reverse=True)
    mut_sum = np.sum([x[1] for x in muts])
    muts = [(t,v/mut_sum) for t,v in muts]
    plt.bar(np.arange(len(muts)), height=[m[1] for m in muts])
    plt.xticks(np.arange(len(muts)), [f"{m[0][0]}->{m[0][1]}" for m in muts], rotation=30)
    plt.ylabel('fraction')
    plt.savefig(args.output_mutations)