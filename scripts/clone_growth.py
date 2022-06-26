import pandas as pd
import argparse,json
from collections import defaultdict
from datetime import datetime
import numpy as np
import matplotlib as mpl
mpl.rcParams['axes.formatter.useoffset'] = False

import matplotlib.pyplot as plt
import seaborn as sns
from root_to_tip import filter_and_transform

if __name__=="__main__":
    parser = argparse.ArgumentParser(
        description="remove time info",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument('--metadata', type=str, required=True, help="input data")
    parser.add_argument('--clade-gts', type=str, required=True, help="input data")
    parser.add_argument('--clade', type=str, required=True, help="input data")
    parser.add_argument('--min-date', type=float, help="input data")
    parser.add_argument('--output-plot', type=str, help="plot file")
    #parser.add_argument('--output-json', type=str, help="rate file")
    args = parser.parse_args()

    with open(args.clade_gts) as fh:
        clade_gt = json.load(fh)[args.clade]

    d = pd.read_csv(args.metadata, sep='\t').fillna('')
    filtered_data = filter_and_transform(d, clade_gt, min_date=args.min_date, max_date=args.min_date + 0.3, completeness=0)

    intra_subs_dis = filtered_data["divergence"].value_counts().sort_index()
    intra_aaSubs_dis = filtered_data["aaDivergence"].value_counts().sort_index()
    intra_geno = filtered_data["intra_substitutions_str"].value_counts()


    ls = ['-', '--', '-.', ':']
    plt.figure()
    ls_counter = defaultdict(int)
    for x,i in intra_geno.items():
        nmuts = len(x.split(',')) if x else 0
        if (nmuts==0 or i>len(filtered_data)/1000) and ls_counter[nmuts]<5:
            ind = filtered_data["intra_substitutions_str"]==x
            if nmuts>5:
                continue
            plt.plot(sorted(filtered_data.loc[ind,"numdate"]), np.arange(i)+1, c=f'C{nmuts}', lw=2 if nmuts else 3,
                            ls=ls[ls_counter[nmuts]%len(ls)], label=f"gt={x if x else 'founder'}"[:30])
            ls_counter[nmuts]+=1
    plt.ylim(0.5,1000)
    plt.yscale('log')
    plt.legend()

    if args.output_plot:
        plt.savefig(args.output_plot)
    else:
        plt.show()
