import pandas as pd
import argparse,json
from treetime.utils import numeric_date
from datetime import datetime
import numpy as np
from scipy.stats import linregress, scoreatpercentile
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns
mpl.rcParams['axes.formatter.useoffset'] = False

date_reference = datetime(2020,1,1).toordinal()
def date_to_week_since2020(d):
    return (d.toordinal() - date_reference)//7


def filter_and_transform(d, clade_gt, min_date=None, max_date=None, completeness=None):
    # filter for incomplete data
    d = d.loc[d.date.apply(lambda x:len(x)==10 and 'X' not in x)]
    d['datetime'] = d.date.apply(lambda x: datetime.strptime(x, '%Y-%m-%d'))
    d['numdate'] = d.datetime.apply(lambda x: numeric_date(x))
    d['CW'] = d.datetime.apply(date_to_week_since2020)
    # filter date range
    if min_date:
        d = d.loc[d.numdate>args.min_date]
    if max_date:
        d = d.loc[d.numdate<args.max_date]

    # look for clade defining substitutions
    d["clade_substitutions"] = d.substitutions.apply(lambda x:     [y for y in x.split(',') if y in clade_gt['nuc']] if x else [])
    # assign number of substitutions missing in the clade definition
    d["missing_subs"] = d.clade_substitutions.apply(lambda x: len(clade_gt['nuc'])-len(x))

    # define "with-in clade substitutions"
    d["intra_substitutions"] = d.substitutions.apply(lambda x:     [y for y in x.split(',') if y not in clade_gt['nuc']] if x else [])
    # make a hashable string representation
    d["intra_substitutions_str"] = d.intra_substitutions.apply(lambda x: ','.join(x))
    # define "with-in clade substitutions"
    d["intra_aaSubstitutions"] = d.aaSubstitutions.apply(lambda x: [y for y in x.split(',') if y not in clade_gt['aa'] and 'ORF9' not in y] if x else [])

    # within clade divergence
    d["divergence"] =   d.intra_substitutions.apply(lambda x:   len(x))
    d["aaDivergence"] = d.intra_aaSubstitutions.apply(lambda x: len(x))

    # filter
    if completeness is not None:
        return d.loc[d.missing_subs<=completeness]

    return d

def regression_by_week(d, field, min_count=5):
    val = d.loc[:,[field,'CW']].groupby('CW').mean()
    std = d.loc[:,[field,'CW']].groupby('CW').std()
    count = d.loc[:,[field,'CW']].groupby('CW').count()
    ind = count[field]>min_count
    reg = linregress(val.index[ind], val.loc[ind, field])
    slope = reg.slope * 365 / 7
    intercept = reg.intercept - 2020*slope
    return {"slope":slope, "intercept":intercept}

if __name__=="__main__":
    parser = argparse.ArgumentParser(
        description="remove time info",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument('--metadata', type=str, required=True, help="input data")
    parser.add_argument('--clade-gts', type=str, required=True, help="input data")
    parser.add_argument('--clade', type=str, required=True, help="input data")
    parser.add_argument('--min-date', type=float, help="input data")
    parser.add_argument('--max-date', type=float, help="input data")
    parser.add_argument('--output-plot', type=str, help="plot file")
    args = parser.parse_args()

    with open(args.clade_gts) as fh:
        clade_gt = json.load(fh)[args.clade]

    d = pd.read_csv(args.metadata, sep='\t').fillna('')
    filtered_data = filter_and_transform(d, clade_gt, min_date=args.min_date, max_date=args.max_date, completeness=0)

    regression = linregress(filtered_data.numdate, filtered_data.divergence)
    filtered_data["residuals"] = filtered_data.apply(lambda x: x.divergence - (regression.intercept + regression.slope*x.numdate), axis=1)
    iqd = scoreatpercentile(filtered_data.residuals, 75) -  scoreatpercentile(filtered_data.residuals, 25)
    filtered_data["outlier"] = filtered_data.residuals.apply(lambda x: np.abs(x)>5*iqd)
    ind = filtered_data.outlier==False
    # regression_clean = linregress(filtered_data.numdate[ind], filtered_data.divergence[ind])
    # regression_clean_aa = linregress(filtered_data.numdate[ind], filtered_data.aaDivergence[ind])
    regression_clean = regression_by_week(filtered_data.loc[ind], "divergence")
    regression_clean_aa = regression_by_week(filtered_data.loc[ind], "aaDivergence")


    intra_subs_dis = filtered_data["divergence"].value_counts().sort_index()
    intra_aaSubs_dis = filtered_data["aaDivergence"].value_counts().sort_index()
    intra_geno = filtered_data["intra_substitutions_str"].value_counts()

    fig, axs = plt.subplots(1,2, figsize=(12,6))
    ymax = 20
    sns.histplot(x=filtered_data.numdate, y=np.minimum(ymax*1.5, filtered_data.divergence), bins=(20,ymax+1), ax=axs[0])
    x = np.array(axs[0].get_xlim())
    axs[0].plot(x, regression_clean["intercept"] + regression_clean["slope"]*x, lw=4, label=f"slope = {regression_clean['slope']:1.1f} subs/year")
    axs[0].legend(loc=2)
    axs[0].set_ylim(0,ymax)


    ymax = 20
    sns.histplot(x=filtered_data.numdate, y=np.minimum(ymax*1.5, filtered_data.aaDivergence), bins=(20,ymax+1), ax=axs[1])
    x = np.array(axs[0].get_xlim())
    axs[1].plot(x, regression_clean_aa["intercept"] + regression_clean_aa["slope"]*x, lw=4, label=f"slope = {regression_clean_aa['slope']:1.1f} subs/year")
    axs[1].legend(loc=2)
    axs[1].set_ylim(0,ymax)

    if args.output_plot:
        plt.savefig(args.output_plot)
    else:
        plt.show()
