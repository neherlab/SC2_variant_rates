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

columns = ['seqName', 'Nextclade_pango', "privateNucMutations.unlabeledSubstitutions", 'qc.overallStatus', 'privateNucMutations.reversionSubstitutions']

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
    parser.add_argument('--rare-cutoff', type=int, default=2, help="minimal number of occurrences to count mutations")
    parser.add_argument('--output-fitness', type=str, required=True, help="fitness table")
    parser.add_argument('--output-events', type=str, required=True, help="events table")
    parser.add_argument('--output-mutations', type=str, required=True, help="fitness table")
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
    pango_counter = defaultdict(int)
    for r, row in d.iterrows():
        pango = row.Nextclade_pango
        if pango[0]=='X':
            continue

        pango_counter[pango] += 1
        muts = row["privateNucMutations.unlabeledSubstitutions"].split(',')
        rev_muts = row["privateNucMutations.reversionSubstitutions"].split(',')
        if len(muts) + len(rev_muts):
            for m in muts + rev_muts:
                if m:
                    mutation_counter[pango][m] += 1

    lineage_size_cutoff = 100
    big_lineages = [pango for pango in pango_counter if pango_counter[pango]>lineage_size_cutoff]
    nlin = len(big_lineages)
    tmp_number_of_pango_muts = defaultdict(int)
    for pango, muts in pango_gts.items():
        if pango_counter[pango]>lineage_size_cutoff:
            for m in muts['nuc']:
                tmp_number_of_pango_muts[int(m[1:-1])-1] += 1
    number_of_pango_muts = np.array([tmp_number_of_pango_muts.get(pos,0)/nlin for pos in range(len(ref))])


    # for each pango lineage, count mutations by position, as well as synoymous and non-synonyomous by codon
    position_counter = defaultdict(lambda: defaultdict(list))
    position_freq = defaultdict(lambda: defaultdict(list))
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
            if pango_counter[pango]>lineage_size_cutoff and c>1:
                position_freq[pos][(a,d)].append(c/pango_counter[pango])

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
        return sum([len([y for y in x if y>=cutoff]) for x in l])

    sorted_transitions = []
    for a in 'ACGT':
        for d in 'ACGT':
            if a!=d:
                sorted_transitions.append((a,d))

    all_events = []
    for pos in range(len(ref)):
        tmp = []
        for a, d in sorted_transitions:
            tmp.append(total_len_filtered_lists([position_counter[pos][(a,d)]], args.rare_cutoff))
        all_events.append(tmp)

    total_events = np.array([total_len_filtered_lists(position_counter[pos].values(), args.rare_cutoff)
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
    total_rate = np.sum(list(events_by_transition.values()))
    for a,d in events_by_transition:
        if a in base_content:
            transition_rate[(a,d)] = events_by_transition[(a,d)]/base_content[a]/total_rate*len(ref)
            away_rate[a] += transition_rate[(a,d)]

    # average mutation rate used to scale rates
    total_events_rescaled = np.array([c/away_rate[nuc] for c,nuc in zip(total_events, ref.seq)])
    total_events_rescaled /= np.median(total_events_rescaled)

    with open(args.output_events, 'w') as fh:
        fh.write('\t'.join(['position', 'ref_state', 'lineage_fraction_with_changes', 'gene', 'codon', 'pos_in_codon']
                            + [f"{a}->{d}" for a,d in sorted_transitions]) + '\n')
        for pos, nuc in enumerate(ref_array):
            gene = map_to_gene.get(pos,"")
            data = "\t".join([f"{all_events[pos][i]}" for i in range(len(sorted_transitions))])
            if gene:
                gene_pos = pos - gene_position[gene].start
                codon = gene_pos//3 + 1
                cp = gene_pos%3 + 1
                fh.write(f'{pos+1}\t{nuc}\t{number_of_pango_muts[pos]:1.3f}\t{gene}\t{codon}\t{cp}\t' + data + "\n")
            else:
                fh.write(f'{pos+1}\t{nuc}\t{number_of_pango_muts[pos]:1.3f}\t\t\t\t' + data + "\n")

    with open(args.output_fitness, 'w') as fh:
        fh.write('\t'.join(['position', 'ref_state', 'lineage_fraction_with_changes', 'gene', 'codon', 'pos_in_codon', 'total_count', 'total_events', 'tolerance']) + '\n')
        for pos, nuc in enumerate(ref_array):
            gene = map_to_gene.get(pos,"")
            if gene:
                gene_pos = pos - gene_position[gene].start
                codon = gene_pos//3 + 1
                cp = gene_pos%3 + 1
                fh.write(f'{pos+1}\t{nuc}\t{number_of_pango_muts[pos]:1.3f}\t{gene}\t{codon}\t{cp}\t{total_counts[pos]}\t{total_events[pos]}\t{total_events_rescaled[pos]:1.3f}\n')
            else:
                fh.write(f'{pos+1}\t{nuc}\t{number_of_pango_muts[pos]:1.3f}\t\t\t\t{total_counts[pos]}\t{total_events[pos]}\t{total_events_rescaled[pos]:1.3f}\n')

    with open(args.output_mutations, 'w') as fh:
        fh.write('\t'.join(['mutation', 'raw_counts', 'origin_sites', 'scaled_rate']) + '\n')
        for a in 'ACTG':
            for d in 'ACGT':
                if a==d: continue
                key = (a,d)
                fh.write(f"{a}->{d}\t{events_by_transition[key]}\t{base_content[a]}\t{transition_rate[key]:1.3f}\n")


    # # calculate synonymous and non-synonymous distributions
    # total_events_by_type = {}
    # for mut_counts, label in [(syn, 'syn'), (nonsyn, 'nonsyn')]:
    #     total_events_by_type[label] = {}
    #     for gene in mut_counts:
    #         total_events_by_type[label][gene] = [np.sum([len([x for x in mut_counts[gene][pos][a] if x>cutoff])/away_rate[a]
    #                                                  for a in mut_counts[gene][pos]])
    #                                                for pos in range(gene_length[gene])]


    # def calc_fitness_cost(freqs, transition):
    #     mu = 0.0004/50
    #     mut_rate = transition_rate[transition]*mu
    #     avg_freq = np.sum(freqs)/nlin
    #     return mut_rate/(avg_freq + 1e-3/nlin)

    # fitness_cost = []
    # for pos,nuc in enumerate(ref_array):
    #     total_muts = np.sum([len(v) for v in position_freq[pos].values()])
    #     reference_muts = np.sum([len(v) for (a,d), v in position_freq[pos].items() if a==nuc])
    #     if reference_muts>0.9*total_muts:
    #         fitness_cost.append([np.inf if d==nuc
    #                              else calc_fitness_cost(position_freq[pos][(nuc,d)],(nuc,d))
    #                             for d in 'ACGT'])
    #     else:
    #         fitness_cost.append([np.nan, np.nan, np.nan, np.nan])

