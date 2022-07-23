import pandas as pd
import numpy as np
import argparse

columns = ['strain', 'virus', 'gisaid_epi_isl', 'genbank_accession',
'date', 'region', 'country', 'length', 'host', 'Nextstrain_clade', 'pango_lineage',
'Nextclade_pango', 'missing_data', 'divergence', 'nonACGTN',
'rare_mutations', 'reversion_mutations', 'potential_contaminants',
'QC_overall_score', 'QC_overall_status', 'frame_shifts', 'deletions',
'insertions', 'substitutions', 'aaSubstitutions', 'clock_deviation']


if __name__=="__main__":
    parser = argparse.ArgumentParser(
        description="remove time info",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument('--metadata', type=str, required=True, help="input data")
    parser.add_argument('--variant-labels', nargs='+', type=str, required=True, help="input data")
    parser.add_argument('--variants', nargs='+', type=str, required=True, help="input data")
    args = parser.parse_args()

    files = [f"subsets/{v}.tsv" for v in args.variant_labels]

    d = pd.read_csv(args.metadata, sep='\t', usecols=columns).fillna('')
    d["Nextstrain_clade"] = d.Nextstrain_clade.apply(lambda x:x.split()[0] if x else '')

    for v, name in zip(args.variants, files):
        if len(v.split(','))>1:
            ind = np.any([d.Nextstrain_clade==y for y in v.split(',')], axis=0)
        else:
            ind = d.Nextstrain_clade==v

        subset = d.loc[ind]
        print(name, v, len(subset))
        subset.to_csv(name, sep='\t')

