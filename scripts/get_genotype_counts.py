import pandas as pd
import argparse,json
from datetime import datetime
import numpy as np
from root_to_tip import filter_and_transform, get_clade_gts
from collections import defaultdict

if __name__=="__main__":
    parser = argparse.ArgumentParser(
        description="remove time info",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument('--metadata', nargs='+', type=str, required=True, help="input data")
    parser.add_argument('--clade-gts', type=str, required=True, help="input data")
    parser.add_argument('--clade', type=str, required=True, help="input data")
    parser.add_argument('--sub-clades', type=str, required=True, help="input data")
    parser.add_argument('--min-date', type=float, help="input data")
    parser.add_argument('--max-date', type=float, help="input data")
    parser.add_argument('--query', type=str, help="filters")
    parser.add_argument('--bin-size', type=float, help="input data")
    parser.add_argument('--output-json', type=str, help="rate file")
    args = parser.parse_args()

    clade_gt = get_clade_gts(args.clade_gts, args.sub_clades)

    d = pd.concat([pd.read_csv(x, sep='\t').fillna('') for x in args.metadata])
    filtered_data = filter_and_transform(d, clade_gt, min_date=args.min_date, max_date=args.max_date,
                                         query = args.query,
                                         completeness=0, swap_root=args.clade=='19B+')
    filtered_data["day"] = filtered_data.datetime.apply(lambda x:x.toordinal())

    print("clade", args.clade, "done filtering")
    bins = np.arange(np.min(filtered_data.day),np.max(filtered_data.day), args.bin_size)
    all_sequences = np.histogram(filtered_data.day, bins=bins)[0]
    cumulative_sum = np.zeros_like(bins[:-1])

    mutation_number = {}
    nmax = 5
    for n in range(nmax):
        ind = filtered_data["divergence"]==n
        mutation_number[n] = np.histogram(filtered_data.loc[ind,"day"], bins=bins)[0]
        cumulative_sum += mutation_number[n]
    mutation_number[f'{nmax}+'] = all_sequences - cumulative_sum
    print("clade", args.clade, "done mutation_number")

    mutation_counts = defaultdict(int)
    cutoff = args.min_date + (args.max_date - args.min_date)/2
    window = cutoff - 14/365, cutoff + 14/365

    ind_early = filtered_data.numdate<cutoff
    ind_window = (filtered_data.numdate>=window[0]) & (filtered_data.numdate<window[1])
    n_early = ind_early.sum()
    n_window = ind_window.sum()

    # Find common mutations
    for muts in filtered_data.loc[ind_early, "intra_substitutions"]:
        for m in muts:
            mutation_counts[m] +=1

    relevant_muts = [x[0] for x in sorted(list(mutation_counts.items()), key=lambda k:k[1])[-10:]]
    mutations = {}
    nmax = 5
    for mut in relevant_muts:
        ind = filtered_data["intra_substitutions"].apply(lambda x: mut in x)
        mutations[mut] = np.histogram(filtered_data.loc[ind,"day"], bins=bins)[0]

    # Get mutation _spectrum
    mutation_spectrum = defaultdict(float)
    for muts in filtered_data.loc[ind_window, "intra_substitutions"]:
        for m in muts:
            mutation_spectrum[m] += 1.0/n_window

    relevant_muts = [x[0] for x in sorted(list(mutation_counts.items()), key=lambda k:k[1])[-10:]]
    mutations = {}
    nmax = 5
    for mut in relevant_muts:
        ind = filtered_data["intra_substitutions"].apply(lambda x: mut in x)
        mutations[mut] = np.histogram(filtered_data.loc[ind,"day"], bins=bins)[0]

    print("clade", args.clade, "done mutations")

    intra_geno = filtered_data.loc[ind_early,"intra_substitutions_str"].value_counts()
    genotypes = {}
    for x,i in intra_geno.items():
        if i<10 or i<n_early/10000:
            break
        nmuts = len(x.split(',')) if x else 0
        ind = filtered_data["intra_substitutions_str"]==x
        if nmuts==0 or (nmuts==1 and ind.sum()>0.02*n_early) or (nmuts<4 and ind.sum()>0.05*n_early):
            genotypes[x] = np.histogram(filtered_data.loc[ind,"day"], bins=bins)[0]

    print("clade", args.clade, "done genotypes")
    with open(args.output_json, 'w') as fh:
        json.dump({'bins':[int(x) for x in bins],
                   'all_samples':[int(x) for x in all_sequences],
                   'mutation_number': {k:[int(x) for x in v] for k,v in mutation_number.items()},
                   'mutation_spectrum': {k: float(v) for k,v in mutation_spectrum.items()},
                   'mutations': {k:[int(x) for x in v] for k,v in mutations.items()},
                   'genotypes': {k:[int(x) for x in v] for k,v in genotypes.items()}}, fh)
