variants = {'19A':'19A',
            '19B':'19B',
            '20A':'20A',
            '20B':'20B',
            '20C':'20C',
            '20A+':'"20A,20B,20C,20D"',
            '20E':'20E',
            '20H':'20H',
            '20I':'20I',
            '20J':'20J',
            '21D':'21D',
            '21G':'21G',
            '21H':'21H',
            '21I':'21I',
            '21J':'21J',
            '21K':'21K',
            '21L':'21L',
            '22A':'22A',
            '22B':'22B',
#            '22D':'22D',
#            '21L+':'"21L,22C,22D"'
}

offset = 0.5
date_ranges = {
'19A':  (2019.9, 2019.9 + offset),
'19B':  (2019.9, 2019.9 + offset),
'20A':  (2020.1, 2020.1 + offset),
'20B':  (2020.1, 2020.1 + offset),
'20C':  (2020.1, 2020.1 + offset),
'20A+': (2020.1, 2020.1 + offset),
'20E':  (2020.5, 2020.5 + offset),
'20H':  (2020.7, 2020.7 + offset),
'20I':  (2020.7, 2020.7 + offset),
'20J':  (2020.7, 2020.7 + offset),
'21G':  (2021.0, 2021.2 + offset),
'21H':  (2021.0, 2021.2 + offset),
'21I':  (2021.2, 2021.2 + offset),
'21J':  (2021.2, 2021.2 + offset),
'21K':  (2021.8, 2021.8 + offset),
'21L':  (2021.8, 2021.8 + offset),
'22A':  (2022.0, 2022.0 + offset),
'22B':  (2022.0, 2022.0 + offset),
'22D':  (2022.2, 2022.2 + offset),
'21L+': (2021.8, 2021.8 + offset),
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
        curl https://data.nextstrain.org/nextclade_sars-cov-2_root-sequence.json | gunzip > {output.root}
        """


rule split_by_variant:
    input:
        metadata = "data/metadata.tsv.gz"
    output:
        files = [f"subsets/{v}.tsv" for v in variants]
    params:
        variants = [x for x in variants],
        variant_labels = [variants[x] for x in variants]
    shell:
        """
        python3 scripts/split_by_variant.py --metadata {input.metadata} --variants {params.variants} --variant-labels {params.variant_labels}
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
        clades = lambda w: variants[w.v]
    shell:
        """
        python3 scripts/root_to_tip.py --metadata {input.metadata} --clade {params.clade} --sub-clades {params.clades} \
                                       --clade-gts data/clade_gts.json \
                                       --min-date {params.mindate} \
                                       --max-date {params.maxdate} \
                                       --output-plot {output.figure} \
                                       --output-json {output.json}
        """

rule genotype_counts:
    input:
        gt = "data/clade_gts.json",
        metadata = "subsets/{v}.tsv"
    output:
        json = "genotypes/{v}_counts.json"
    params:
        clade = lambda w: w.v,
        clades = lambda w: variants[w.v],
        mindate = lambda w: date_ranges[w.v][0],
        maxdate = lambda w: date_ranges[w.v][1],
        bin_size = 5
    shell:
        """
        python3 scripts/get_genotype_counts.py --metadata {input.metadata} --clade {params.clade} --sub-clades {params.clades} \
                                       --clade-gts data/clade_gts.json \
                                       --min-date {params.mindate} \
                                       --max-date {params.maxdate} \
                                       --bin-size {params.bin_size} \
                                       --output-json {output.json}
        """

rule genotype_count_figures:
    input:
        json = "genotypes/{v}_counts.json"
    output:
        fig = "figures/{v}_counts.png"
    shell:
        """
        python3 scripts/plot_genotype_counts.py --counts {input.json} --output-plot {output.fig}
        """


rule clone_growth:
    input:
        gt = "data/clade_gts.json",
        metadata = "subsets/{v}.tsv"
    output:
        figure = "figures/{v}_clones.png"
    params:
        clade = lambda w: w.v,
        clades = lambda w: variants[w.v],
        mindate = lambda w: date_ranges[w.v][0]
    shell:
        """
        python3 scripts/clone_growth.py --metadata {input.metadata} --clade {params.clade} --sub-clades {params.clades} \
                                       --clade-gts data/clade_gts.json \
                                       --min-date {params.mindate} \
                                       --output-plot {output.figure}
        """


rule all_rtt:
    input:
        expand("figures/{v}_rtt.png", v=variants.keys())

rule all_clones:
    input:
        expand("figures/{v}_clones.png", v=variants.keys())

rule genotypes:
    input:
        expand("figures/{v}_counts.png", v=variants.keys())

rule rate_table:
    input:
        rate_files = expand("rates/{v}_rate.json", v=variants.keys()),
        gt = "data/clade_gts.json"
    output:
        rate_table = "rates.tsv",
        rate_table_tex = "rates.tex"
    run:
        import json
        import pandas as pd
        from treetime.utils import datestring_from_numeric
        with open(input.gt) as fin:
            clade_gts = json.load(fin)

        data = []
        for fname in input.rate_files:
            with open(fname) as fin:
                d = json.load(fin)
            base_clade = d['clade'][:3]
            aa_div = len([x for x in clade_gts[base_clade]['aa'] if 'ORF9' not in x])
            nuc_div = len(clade_gts[base_clade[:3]]['nuc'])
            data.append({'clade':d['clade'], 'nuc_rate': d['nuc']['slope'],
                         'nuc_origin': d['nuc']['origin'],'nuc_origin_date': datestring_from_numeric(d['nuc']['origin']),
                         'aa_rate': d['aa']['slope'],
                         'aa_origin':d['aa']['origin'], 'aa_origin_date':datestring_from_numeric(d['aa']['origin']),
                         'syn_rate': d['syn']['slope'],
                         'syn_origin':d['syn']['origin'],'syn_origin_date':datestring_from_numeric(d['syn']['origin']),
                         'nuc_div': nuc_div, 'aa_div':aa_div, 'syn_div':nuc_div-aa_div})

        df = pd.DataFrame(data)
        df.to_csv(output.rate_table, sep='\t')
        df.to_latex(output.rate_table_tex, float_format="%.2f", index=False)


rule rate_summary:
    input:
        rate_table = "rates.tsv"
    output:
        figure = "figures/rate_summary.pdf",
        figure_rates = "figures/rate_summary.pdf"
    shell:
        """
        python3 scripts/combine_fits.py --rate-table {input.rate_table}\
                                       --output-plot {output.figure} \
                                       --output-plot-rates {output.figure_rates}
        """
