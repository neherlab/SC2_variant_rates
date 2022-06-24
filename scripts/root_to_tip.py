import pandas as pd
import argparse
from treetime.utils import numeric_date
from datetime import datetime
from scipy.stats import linregress

if __name__=="__main__":
    parser = argparse.ArgumentParser(
        description="remove time info",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument('--metadata', type=str, required=True, help="input data")
    args = parser.parse_args()

    d = pd.read_csv(args.metadata, sep='\t').fillna('')
    d = d.loc[d.date.apply(lambda x:len(x)==10 and 'X' not in x)]
    d.loc[:,'numdate'] = d.date.apply(lambda x: numeric_date(datetime.strptime(x, '%Y-%m-%d')))
    d["divergence"] = d.substitutions.apply(lambda x: len(x.split(',')))
    d["aaDivergence"] = d.aaSubstitutions.apply(lambda x: len(x.split(',')))

    regression = linregress(d.numdate, d.divergence)

    plt.scatter(d.numdate, d.divergence)
