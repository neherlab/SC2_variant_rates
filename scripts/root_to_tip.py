import pandas as pd
import argparse,json
from treetime.utils import numeric_date, datestring_from_numeric
from datetime import datetime
import numpy as np
from scipy.stats import linregress, scoreatpercentile
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns
from collections import defaultdict
mpl.rcParams['axes.formatter.useoffset'] = False

date_reference = datetime(2020,1,1).toordinal()
def date_to_week_since2020(d):
    return (d.toordinal() - date_reference)//7

def week_since2020_to_date(d):
    return datetime.fromordinal(int(d*7 + date_reference))

def week_since2020_to_numdate(d):
    return numeric_date(week_since2020_to_date(d))

def filter_and_transform(d, clade_gt, min_date=None, max_date=None, query=None, completeness=None, swap_root=False, max_group=None):
    # filter for incomplete data
    d = d.loc[d.date.apply(lambda x:len(x)==10 and 'X' not in x)]
    d = d.loc[d.QC_overall_status=='good']
    if query:
        d = d.query(query)
    d['datetime'] = d.date.apply(lambda x: datetime.strptime(x, '%Y-%m-%d'))
    d['numdate'] = d.datetime.apply(lambda x: numeric_date(x))
    d['CW'] = d.datetime.apply(date_to_week_since2020)
    # filter date range
    if min_date:
        d = d.loc[d.numdate>min_date]
    if max_date:
        d = d.loc[d.numdate<max_date]

    # look for clade defining substitutions
    d["clade_substitutions"] = d.substitutions.apply(lambda x:     [y for y in x.split(',') if y in clade_gt['nuc']] if x else [])
    # assign number of substitutions missing in the clade definition
    d["missing_subs"] = d.clade_substitutions.apply(lambda x: len(clade_gt['nuc'])-len(x))

    # define "with-in clade substitutions"
    d["intra_substitutions"] = d.substitutions.apply(lambda x:     [y for y in x.split(',')
                                                     if all([y not in clade_gt['nuc'], int(y[1:-1])>150,  int(y[1:-1])<29753])] if x else [])
    # define "with-in clade substitutions"
    d["intra_aaSubstitutions"] = d.aaSubstitutions.apply(lambda x: [y for y in x.split(',') if y not in clade_gt['aa'] and 'ORF9' not in y] if x else [])
    d["intra_SpikeSubstitutions"] = d.aaSubstitutions.apply(lambda x: [y for y in x.split(',') if y not in clade_gt['aa'] and 'ORF9' not in y and y[0]=='S'] if x else [])

    if swap_root:
        muts = [("C8782T","T8782C"), ("T28144C","C28144T")]
        def swap(mutations, pair):
            return [y for y in mutations if y!=pair[0]] if pair[0] in mutations else mutations + [pair[1]]
        for m in muts:
            d["intra_substitutions"] = d.intra_substitutions.apply(lambda x: swap(x,m))
        aa_mut = ("ORF8:L84S","ORF8:S84L")
        d["intra_aaSubstitutions"] = d.intra_aaSubstitutions.apply(lambda x: swap(x,aa_mut))


    # make a hashable string representation
    d["intra_substitutions_str"] = d.intra_substitutions.apply(lambda x: ','.join(x))
    # within clade divergence
    d["divergence"] =   d.intra_substitutions.apply(lambda x:   len(x))
    d["aaDivergence"] = d.intra_aaSubstitutions.apply(lambda x: len(x))
    d["spikeDivergence"] = d.intra_SpikeSubstitutions.apply(lambda x: len(x))
    d["synDivergence"] = d["divergence"] - d["aaDivergence"]


    # filter
    if completeness is not None:
        d = d.loc[d.missing_subs<=completeness]

    if max_group:
        return d.groupby(['CW', 'country']).sample(max_group, replace=True).drop_duplicates(subset='strain')

    return d

def weighted_regression(x,y,w):
    '''
    This function determine slope and intercept by minimizing
    sum_i w_i*(y_i - f(x_i))^2 with f(x) = slope*x + intercept

    sum_i w_i*(y_i - f(x_i)) = 0 => sum_i w_i y_i = intercept * sum_i w_i + slope sum_i w_i x_i
    sum_i w_i*(y_i - f(x_i)) x_i = 0 => sum_i w_i y_i x_i = intercept * sum_i w_i x_i + slope sum_i w_i x_i^2
    '''
    wa=np.array(w)
    xa=np.array(x)
    ya=np.array(y)
    wx = np.sum(xa*wa)
    wy = np.sum(ya*wa)
    wxy = np.sum(xa*ya*wa)
    wxx = np.sum(xa**2*wa)
    wsum = np.sum(wa)
    wmean = np.mean(wa)

    slope = (wy*wx - wxy*wsum)/(wx**2 - wxx*wsum)
    intercept = (wy - slope*wx)/wsum

    # not correct
    # hessianinv = np.linalg.inv(np.array([[wxx, wx], [wx, wsum]])/wsum)
    # stderrs = wmean*np.sqrt(hessianinv.diagonal())

    return {"slope": slope, "intercept":intercept} #, 'slope_err':stderrs[0], 'intercept_err':stderrs[1]}


def regression_by_week(d, field, min_count=5):
    val = d.loc[:,[field,'CW']].groupby('CW').mean()
    std = d.loc[:,[field,'CW']].groupby('CW').std()
    count = d.loc[:,[field,'CW']].groupby('CW').count()
    ind = count[field]>min_count
    #reg = linregress(val.index[ind], val.loc[ind, field])
    # slope = reg.slope * 365 / 7
    # intercept = reg.intercept - 2020*slope
    reg = weighted_regression(val.index[ind], val.loc[ind, field], np.array((count[ind]-min_count)**0.25).squeeze())
    slope = reg["slope"] * 365 / 7
    intercept = reg["intercept"] - 2020*slope

    return {"slope":slope, "intercept":intercept,
            "origin": -intercept/slope,
            "date":[week_since2020_to_numdate(x) for x in val.index[ind]],
            "mean":[x for x in val.loc[ind, field]],
            "stderr":[x for x in std.loc[ind, field]],
            "count":[x for x in count.loc[ind, field]]}

def make_date_ticks(ax):
    ax.set_xlabel('')
    ax.set_xticklabels([datestring_from_numeric(x) for x in ax.get_xticks()], rotation=30, horizontalalignment='right')

def get_clade_gts(all_gts, subclade_str):
    with open(all_gts) as fh:
        clade_gts = json.load(fh)

    subclades = subclade_str.split(',')
    if len(subclades)>1:
        clade_gt = {'nuc':{}, 'aa':{}}
        clade_gt['nuc'] = set.intersection(*[set(clade_gts[x]['nuc']) for x in subclades])
        clade_gt['aa'] = set.intersection(*[set(clade_gts[x]['aa']) for x in subclades])
    else:
        clade_gt = clade_gts[subclade_str]

    return clade_gt

if __name__=="__main__":
    parser = argparse.ArgumentParser(
        description="remove time info",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument('--metadata',  type=str, required=True, help="input data")
    parser.add_argument('--clade-gts', type=str, required=True, help="input data")
    parser.add_argument('--clade', type=str, required=True, help="input data")
    parser.add_argument('--sub-clades', type=str, required=True, help="input data")
    parser.add_argument('--min-date', type=float, help="input data")
    parser.add_argument('--max-date', type=float, help="input data")
    parser.add_argument('--max-group', type=int, default=100, help="input data")
    parser.add_argument('--query', type=str, help="filters")
    parser.add_argument('--output-plot', type=str, help="plot file")
    parser.add_argument('--output-json', type=str, help="rate file")
    args = parser.parse_args()

    clade_gt = get_clade_gts(args.clade_gts, args.sub_clades)

    d = pd.read_csv(args.metadata, sep='\t').fillna('')
    filtered_data = filter_and_transform(d, clade_gt, min_date=args.min_date, max_date=args.max_date,
                                         query = args.query, max_group=args.max_group,
                                         completeness=0, swap_root=args.clade=='19B+')

    regression = linregress(filtered_data.numdate, filtered_data.divergence)
    filtered_data["residuals"] = filtered_data.apply(lambda x: x.divergence - (regression.intercept + regression.slope*x.numdate), axis=1)
 #   iqd = scoreatpercentile(filtered_data.residuals, 75) -  scoreatpercentile(filtered_data.residuals, 25)
 #   filtered_data["outlier"] = filtered_data.residuals.apply(lambda x: np.abs(x)>5*iqd)

    tolerance = lambda t: 3 + 2*np.sqrt(np.maximum(0,(regression.intercept + regression.slope*t)))
    filtered_data["outlier"] = filtered_data.apply(lambda x: np.abs(x.residuals)>tolerance(x.numdate), axis=1)
    ind = filtered_data.outlier==False
    # regression_clean = linregress(filtered_data.numdate[ind], filtered_data.divergence[ind])
    # regression_clean_aa = linregress(filtered_data.numdate[ind], filtered_data.aaDivergence[ind])
    regression_clean = regression_by_week(filtered_data.loc[ind], "divergence")
    regression_clean_aa = regression_by_week(filtered_data.loc[ind], "aaDivergence")
    regression_clean_syn = regression_by_week(filtered_data.loc[ind], "synDivergence")
    regression_clean_spike = regression_by_week(filtered_data.loc[ind], "spikeDivergence")

    fig, axs = plt.subplots(1,3, figsize=(18,6), sharex=True, sharey=True)
    ymax = 20
    bins = bins=(20,np.arange(-0.5,ymax+0.5))
    sns.histplot(x=filtered_data.numdate, y=np.minimum(ymax*1.5, filtered_data.divergence), bins=bins, ax=axs[0])
    x = np.linspace(*axs[0].get_xlim(),101)
    axs[0].set_title('all differences')
    axs[0].plot(x, regression.intercept + regression.slope*x + tolerance(x), lw=4)
    axs[0].plot(x, regression_clean["intercept"] + regression_clean["slope"]*x, lw=4, label=f"slope = {regression_clean['slope']:1.1f} subs/year")
    axs[0].errorbar(regression_clean["date"], regression_clean["mean"], regression_clean["stderr"])

    axs[1].set_title('amino acid differences')
    sns.histplot(x=filtered_data.numdate[ind], y=np.minimum(ymax*1.5, filtered_data.aaDivergence[ind]), bins=bins, ax=axs[1])
    axs[1].plot(x, regression_clean_aa["intercept"] + regression_clean_aa["slope"]*x, lw=4, label=f"slope = {regression_clean_aa['slope']:1.1f} subs/year")
    axs[1].errorbar(regression_clean_aa["date"], regression_clean_aa["mean"], regression_clean_aa["stderr"])

    axs[2].set_title('synonymous differences')
    sns.histplot(x=filtered_data.numdate[ind], y=np.minimum(ymax*1.5, filtered_data.synDivergence[ind]), bins=bins, ax=axs[2])
    axs[2].plot(x, regression_clean_syn["intercept"] + regression_clean_syn["slope"]*x, lw=4, label=f"slope = {regression_clean_syn['slope']:1.1f} subs/year")
    axs[2].errorbar(regression_clean_syn["date"], regression_clean_syn["mean"], regression_clean_syn["stderr"])

    for ax in axs:
        make_date_ticks(ax)
        ax.set_yticks(np.arange(0,ymax,3))
        ax.legend(loc=2)
        ax.set_ylim(-0.5,ymax-0.5)

    if args.output_plot:
        plt.savefig(args.output_plot)
    else:
        plt.show()

    rate_data = {'clade':args.clade, 'nuc':regression_clean, 'aa':regression_clean_aa,
                 'syn':regression_clean_syn, 'spike':regression_clean_spike}
    if args.output_json:
        with open(args.output_json, 'w') as fh:
            json.dump(rate_data, fh)