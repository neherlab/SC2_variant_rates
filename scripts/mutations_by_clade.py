import argparse,json
from importlib.metadata import metadata
import pandas as pd
import numpy as np
import matplotlib as mpl
from Bio import SeqIO,Seq
mpl.rcParams['axes.formatter.useoffset'] = False
from collections import defaultdict
from scipy.stats import scoreatpercentile
import seaborn as sns
import pandas as pd

import matplotlib.pyplot as plt
import seaborn as sns

columns = ['seqName', 'clade', "privateNucMutations.unlabeledSubstitutions", 'qc.overallStatus', 'privateNucMutations.reversionSubstitutions']

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

    # load and filter metadata
    metadata = pd.read_csv(args.metadata, sep='\t', usecols=columns).fillna('')
    metadata = metadata.loc[metadata["qc.overallStatus"]=='good',:]

    # glob all rare mutations by pango lineage
    mutation_counter = {clade: defaultdict(lambda: defaultdict(int)) for clade in metadata.clade.unique()}
    clade_counter = defaultdict(int)
    for r, row in metadata.iterrows():
        clade = row.clade

        clade_counter[clade] += 1
        muts = row["privateNucMutations.unlabeledSubstitutions"].split(',')
        rev_muts = row["privateNucMutations.reversionSubstitutions"].split(',')
        if len(muts) + len(rev_muts):
            for m in muts + rev_muts:
                if m:
                    a,pos, d = m[0], int(m[1:-1]), m[-1]
                    if codon_pos[pos-1] == 2:
                        mutation_counter[clade][pos][(a,d)] += 1

    mutation_counter_filtered = defaultdict(lambda: defaultdict(int))
    for clade in mutation_counter:
        for pos, entry in sorted(mutation_counter[clade].items(), key=lambda x:sum(list(x[1].values())))[:-20]:
            for k,v in entry.items():
                mutation_counter_filtered[clade][k]+=v


    total = {k: sum(list(v.values())) for k,v in mutation_counter_filtered.items()}
    spectra = {}
    for clade in mutation_counter_filtered:
        if total[clade]<10000: continue
        spectra[clade]={}
        for a in 'ACGT':
            for d in 'ACGT':
                if a!=d:
                    spectra[clade][(a,d)] = mutation_counter_filtered[clade][(a,d)]/total[clade]

    spectra_matrix = np.array([list(spectra[c].values()) for c in spectra])

    mutation_similarity = np.zeros((len(spectra), len(spectra)))
    for i1, c1 in enumerate(spectra):
        for i2, c2 in enumerate(spectra):
            mutation_similarity[i1,i2] = np.sum(np.abs([spectra[c1][k] - spectra[c2][k] for k in spectra[c1]]))

    plt.figure()
    plt.matshow(mutation_similarity)
    plt.yticks(range(len(spectra)), list(spectra.keys()))


    plt.figure()
    for i1, c1 in enumerate(spectra):
        plt.plot(list(spectra[c1].values()), label=c1, ls='-' if i1//10 else '--')
    plt.xticks(range(len(spectra[c1])), list(spectra[c1].keys()), rotation=90)
    plt.legend()


    plt.figure()
    sns.clustermap(pd.DataFrame(mutation_similarity, columns=spectra.keys(), index=spectra.keys()))


    plt.figure()
    sns.clustermap(pd.DataFrame(spectra_matrix, columns=spectra[c1].keys(), index=spectra.keys()))
