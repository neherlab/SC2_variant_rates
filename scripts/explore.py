import pandas as pd
from root_to_tip import filter_and_transform, get_clade_gts
from treetime.utils import numeric_date, datestring_from_numeric
from datetime import datetime
# alpha
clade_gt = get_clade_gts("data/clade_gts.json", "20I")

alpha = pd.read_csv('subsets/20I.tsv', sep='\t')
alpha_filtered = filter_and_transform(alpha, clade_gt, min_date=numeric_date(datetime(2020,8,1)))

print(alpha_filtered.loc[alpha_filtered.date<"2020-10-01",["strain", "date", "intra_substitutions"]])
