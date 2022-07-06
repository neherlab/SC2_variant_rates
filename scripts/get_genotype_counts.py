import pandas as pd
import argparse,json
from datetime import datetime
import numpy as np
from root_to_tip import filter_and_transform


if __name__=="__main__":
    parser = argparse.ArgumentParser(
        description="remove time info",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument('--metadata', nargs='+', type=str, required=True, help="input data")
    parser.add_argument('--clade-gts', type=str, required=True, help="input data")
    parser.add_argument('--clade', type=str, required=True, help="input data")
    parser.add_argument('--min-date', type=float, help="input data")
    parser.add_argument('--max-date', type=float, help="input data")
    parser.add_argument('--bin-size', type=float, help="input data")
    parser.add_argument('--output-json', type=str, help="rate file")
    args = parser.parse_args()

    with open(args.clade_gts) as fh:
        clade_gt = json.load(fh)[args.clade]

    d = pd.concat([pd.read_csv(x, sep='\t').fillna('') for x in args.metadata])
    filtered_data = filter_and_transform(d, clade_gt, min_date=args.min_date, max_date=args.max_date, completeness=0)
    filtered_data["day"] = filtered_data.datetime.apply(lambda x:x.toordinal())
    print("clade", args.clade, "done filtering")
    bins = np.arange(np.min(filtered_data.day),np.max(filtered_data.day), args.bin_size)
    all_sequences = np.histogram(filtered_data.day, bins=bins)[0]
    cumulative_sum = np.zeros_like(bins[:-1])

    mutations = {}
    nmax = 5
    for n in range(nmax):
        ind = filtered_data["divergence"]==n
        mutations[n] = np.histogram(filtered_data.loc[ind,"day"], bins=bins)[0]
        cumulative_sum += mutations[n]
    mutations[f'{nmax}+'] = all_sequences - cumulative_sum

    print("clade", args.clade, "done mutations")
    intra_geno = filtered_data["intra_substitutions_str"].value_counts()
    genotypes = {}
    for x,i in intra_geno.items():
        if i<10 or i<len(filtered_data)/10000:
            break
        nmuts = len(x.split(',')) if x else 0
        ind = filtered_data["intra_substitutions_str"]==x
        if nmuts==0 or (nmuts==1 and ind.sum()>0.003*len(filtered_data)) or (nmuts<3 and ind.sum()>0.01*len(filtered_data)):
            genotypes[x] = np.histogram(filtered_data.loc[ind,"day"], bins=bins)[0]

    print("clade", args.clade, "done genotypes")
    with open(args.output_json, 'w') as fh:
        json.dump({'bins':[int(x) for x in bins],
                   'all_samples':[int(x) for x in all_sequences],
                   'mutations': {k:[int(x) for x in v] for k,v in mutations.items()},
                   'genotypes': {k:[int(x) for x in v] for k,v in genotypes.items()}}, fh)
