variants = ['19A', '19B', '20A', '20B', '20C', '20E', '20H', '20I', '20J', '21I', '21J', '21K', '21L']
date_ranges = {
'19A': (2019.9, 2020.3),
'19B': (2019.9, 2020.3),
'20A': (2020.1, 2020.7),
'20B': (2020.1, 2020.7),
'20C': (2020.1, 2020.7),
'20E': (2020.5, 2020.9),
'20H': (2020.7, 2021.2),
'20I': (2020.7, 2021.2),
'20J': (2020.7, 2021.2),
'21I': (2021.2, 2021.8),
'21J': (2021.2, 2021.8),
'21K': (2021.8, 2022.2),
'21L': (2021.8, 2022.2),
}

rule get_data:
    output:
        "data/metadata.tsv.gz"
    shell:
        """
        curl https://data.nextstrain.org/files/ncov/open/metadata.tsv.gz -o {output}
        """

rule get_clade_data:
    output:
        tree = "data/clade_tree.json",
        root = "data/root-sequence.json"
    shell:
        """
        curl https://data.nextstrain.org/nextclade_sars-cov-2.json | gunzip > {output.tree}
        #curl https://data.nextstrain.org/nextclade_sars-cov-2_root-sequence.json -o {output.root}
        curl "https://nextstrain.org/charon/getDataset?prefix=ncov/open/global/all-time&type=root-sequence" -o {output.root}
        """


rule split_by_variant:
    input:
        metadata = "data/metadata.tsv.gz"
    output:
        files = [f"subsets/{v}.tsv" for v in variants]
    params:
        variants = variants
    shell:
        """
        python3 scripts/split_by_variant.py {input.metadata} {params.variants}
        """


rule clade_genotypes:
    input:
        tree = "data/clade_tree.json",
        root = "data/root-sequence.json"
    output:
        gt = "data/clade_gts.json"
    shell:
        """
        python3 scripts/get_genotypes.py --tree {input.tree} --root {input.root} --output {output.gt}
        """

rule root_to_tip:
    input:
        gt = "data/clade_gts.json",
        metadata = "subsets/{v}.tsv"
    output:
        figure = "figures/{v}_rtt.png",
        json = "rates/{v}_rate.json"
    params:
        clade = lambda w: w.v,
        mindate = lambda w: date_ranges[w.v][0],
        maxdate = lambda w: date_ranges[w.v][1],
    shell:
        """
        python3 scripts/root_to_tip.py --metadata {input.metadata} --clade {params.clade} \
                                       --clade-gts data/clade_gts.json \
                                       --min-date {params.mindate} \
                                       --max-date {params.maxdate} \
                                       --output-plot {output.figure} \
                                       --output-json {output.json}
        """

rule all_rtt:
    input:
        expand("figures/{v}_rtt.png", v=variants)

rule rate_table:
    input:
        rate_files = expand("rates/{v}_rate.json", v=variants),
        gt = "data/clade_gts.json"
    output:
        rate_table = "rates.tsv"
    run:
        import json
        import pandas as pd

        with open(input.gt) as fin:
            clade_gts = json.load(fin)

        data = []
        for fname in input.rate_files:
            with open(fname) as fin:
                d = json.load(fin)
            aa_div = len([x for x in clade_gts[d['clade']]['aa'] if 'ORF9' not in x])
            nuc_div = len(clade_gts[d['clade']]['nuc'])
            data.append({'clade':d['clade'], 'nuc_rate': d['nuc']['slope'], 'nuc_origin':d['nuc']['origin'],
                         'aa_rate': d['aa']['slope'], 'aa_origin':d['aa']['origin'], 'syn_rate': d['nuc']['slope']-d['aa']['slope'],
                         'nuc_div': nuc_div, 'aa_div':aa_div, 'syn_div':nuc_div-aa_div})

        df = pd.DataFrame(data)
        df.to_csv(output.rate_table, sep='\t')
