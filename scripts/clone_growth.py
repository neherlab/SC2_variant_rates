import pandas as pd
import argparse,json
from collections import defaultdict
from datetime import datetime
import numpy as np
import matplotlib as mpl
mpl.rcParams['axes.formatter.useoffset'] = False

import matplotlib.pyplot as plt
import seaborn as sns
from root_to_tip import filter_and_transform, get_clade_gts

if __name__=="__main__":
    parser = argparse.ArgumentParser(
        description="remove time info",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument('--metadata', type=str, nargs='+', required=True, help="input data")
    parser.add_argument('--clade-gts', type=str, required=True, help="input data")
    parser.add_argument('--clade', type=str, required=True, help="input data")
    parser.add_argument('--sub-clades', type=str, required=True, help="input data")
    parser.add_argument('--min-date', type=float, help="input data")
    parser.add_argument('--output-plot', type=str, help="plot file")
    #parser.add_argument('--output-json', type=str, help="rate file")
    args = parser.parse_args()

    clade_gt = get_clade_gts(args.clade_gts, args.sub_clades)

    d = pd.concat([pd.read_csv(x, sep='\t').fillna('') for x in args.metadata])
    filtered_data = filter_and_transform(d, clade_gt, min_date=args.min_date, max_date=args.min_date + 0.3, completeness=0)
    #filtered_data=filtered_data.loc[filtered_data.country!='China']

    intra_subs_dis = filtered_data["divergence"].value_counts().sort_index()
    intra_aaSubs_dis = filtered_data["aaDivergence"].value_counts().sort_index()
    intra_geno = filtered_data["intra_substitutions_str"].value_counts()


    ls = ['-', '--', '-.', ':']
    plt.figure()
    dates = sorted(filtered_data.loc[:,"numdate"])
    plt.plot(dates, (np.arange(len(dates))+1), c='k', lw=3, alpha=0.3, label="all")
    ls_counter = defaultdict(int)
    for x,i in intra_geno.items():
        nmuts = len(x.split(',')) if x else 0
        ind = filtered_data["intra_substitutions_str"]==x
        dates = sorted(filtered_data.loc[ind,"numdate"])
        factor = 3 if 'C8782T' in x else 1
        if dates[0]<args.min_date+0.2:
        #if (nmuts<2 or i>len(filtered_data)/100) and ls_counter[nmuts]<10:
            if nmuts>5:
                continue
            plt.plot(dates, factor*(np.arange(i)+1), c=f'C{nmuts}', lw=2 if nmuts else 3,
                            ls=ls[ls_counter[nmuts]%len(ls)], marker='o' if len(dates)<3 else '',
                            label=f"gt={x if x else 'founder'}"[:30])
            ls_counter[nmuts]+=1
    plt.ylim(0.5,1000)
    plt.xlim(args.min_date,args.min_date+0.3)
    plt.yscale('log')
    plt.legend()

    if args.output_plot:
        plt.savefig(args.output_plot)
    else:
        plt.show()
