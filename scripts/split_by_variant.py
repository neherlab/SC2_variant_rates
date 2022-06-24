import pandas as pd
import sys

columns = ['strain', 'virus', 'gisaid_epi_isl', 'genbank_accession',
'sra_accession', 'date', 'region', 'country', 'division', 'location',
'length', 'host', 'age', 'sex', 'Nextstrain_clade', 'pango_lineage',
'date_submitted', 'sampling_strategy',
'Nextclade_pango', 'missing_data', 'divergence', 'nonACGTN',
'rare_mutations', 'reversion_mutations', 'potential_contaminants',
'QC_missing_data', 'QC_mixed_sites', 'QC_rare_mutations',
'QC_snp_clusters', 'QC_frame_shifts', 'QC_stop_codons',
'QC_overall_score', 'QC_overall_status', 'frame_shifts', 'deletions',
'insertions', 'substitutions', 'aaSubstitutions', 'clock_deviation']



metadata = sys.argv[1]
variants = sys.argv[2:]

files = [f"subsets/{v}.tsv" for v in variants]

d = pd.read_csv(metadata, sep='\t', usecols=columns).fillna('')
d["Nextstrain_clade"] = d.Nextstrain_clade.apply(lambda x:x.split()[0] if x else '')

for v, name in zip(variants, files):
    subset = d.loc[d.Nextstrain_clade==v]
    print(name, v, len(subset))
    subset.to_csv(name, sep='\t')

