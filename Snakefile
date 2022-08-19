variants = {'19A':'19A',
            '19B':'19B',
            '19B+':'19A,19B',
            '19B++':'19A,19B,20A,20B,20C,20D',
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
'19A':  2019.9,
'19B':  2019.9,
'19B+': 2019.9,
'19B++': 2019.9,
'20A':  2020.1,
'20B':  2020.1,
'20C':  2020.1,
'20A+': 2020.1,
'20E':  2020.47,
'20H':  2020.63,
'20I':  2020.715,
'20J':  2020.91,
'21D':  2020.9,
'21G':  2021.0,
'21H':  2021.0,
'21I':  2021.22,
'21J':  2021.21,
'21K':  2021.85,
'21L':  2021.85,
'22A':  2022.163,
'22B':  2022.163,
'22D':  2022.2,
'21L+': 2021.8,
}

filter_queries = {
    '20A':  "region=='Europe' | region=='North America'",
    '20B':  "region=='Europe' | region=='North America'",
    '20C':  "region=='Europe' | region=='North America'",
    '20A+': "region=='Europe' | region=='North America'",
    '20I': "region=='Europe'"
}

rule get_data:
    output:
        "data/metadata.tsv.gz"
    shell:
        """
        curl https://data.nextstrain.org/files/ncov/open/metadata.tsv.gz -o {output}
        """

rule get_nextclade:
    output:
        "data/nextclade.tsv.gz"
    shell:
        """
        curl https://data.nextstrain.org/files/ncov/open/nextclade.tsv.gz -o {output}
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
        variant_labels = [x for x in variants],
        variants = [variants[x] for x in variants]
    shell:
        """
        python3 scripts/split_by_variant.py --metadata {input.metadata} --variants {params.variants} --variant-labels {params.variant_labels}
        """


rule pango_genotypes:
    input:
        tree = "data/clade_tree.json",
        root = "data/root-sequence.json"
    output:
        gt = "data/pango_gts.json"
    shell:
        """
        python3 scripts/get_genotypes_pango.py --tree {input.tree} --root {input.root} --output {output.gt}
        """

rule root_to_tip:
    input:
        gt = "data/clade_gts.json",
        metadata = "subsets/{v}.tsv"
    output:
        figure = "figures/{v}_rtt.pdf",
        json = "rates/{v}_rate.json"
    params:
        clade = lambda w: w.v,
        mindate = lambda w: date_ranges[w.v],
        maxdate = lambda w: date_ranges[w.v] + offset,
        clades = lambda w: variants[w.v],
        filter_query = lambda w: ('--query ' + f'"{filter_queries[w.v]}"') if w.v in filter_queries else ''
    shell:
        """
        python3 scripts/root_to_tip.py --metadata {input.metadata} --clade {params.clade} --sub-clades {params.clades} \
                                       --clade-gts data/clade_gts.json \
                                       --min-date {params.mindate} \
                                       --max-date {params.maxdate} \
                                       {params.filter_query} \
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
        mindate = lambda w: date_ranges[w.v],
        maxdate = lambda w: date_ranges[w.v] + offset,
        bin_size = 5,
        filter_query = lambda w: ('--query ' + f'"{filter_queries[w.v]}"') if w.v in filter_queries else ''
    shell:
        """
        python3 scripts/get_genotype_counts.py --metadata {input.metadata} --clade {params.clade} --sub-clades {params.clades} \
                                       --clade-gts data/clade_gts.json \
                                       --min-date {params.mindate} \
                                       --max-date {params.maxdate} \
                                       {params.filter_query} \
                                       --bin-size {params.bin_size} \
                                       --output-json {output.json}
        """

rule genotype_count_figures:
    input:
        json = "genotypes/{v}_counts.json"
    output:
        fig = "figures/{v}_counts.pdf"
    shell:
        """
        python3 scripts/plot_genotype_counts.py --counts {input.json} --output-plot {output.fig}
        """


rule clone_growth:
    input:
        gt = "data/clade_gts.json",
        metadata = "subsets/{v}.tsv"
    output:
        figure = "figures/{v}_clones.pdf"
    params:
        clade = lambda w: w.v,
        clades = lambda w: variants[w.v],
        mindate = lambda w: date_ranges[w.v]
    shell:
        """
        python3 scripts/clone_growth.py --metadata {input.metadata} --clade {params.clade} --sub-clades {params.clades} \
                                       --clade-gts data/clade_gts.json \
                                       --min-date {params.mindate} \
                                       --output-plot {output.figure}
        """

rule af:
    input:
        count_files = expand("genotypes/{v}_counts.json", v=variants.keys()),
    output:
        af_fig = "figures/mutation_frequencies.pdf",
        rates = "data/rates_poisson.tsv",
    shell:
        """
        python3 scripts/plot_af.py --counts {input.count_files} --output-plot {output.af_fig} --output-rates {output.rates}
        """

rule all_rtt:
    input:
        expand("figures/{v}_rtt.pdf", v=variants.keys())

rule all_clones:
    input:
        expand("figures/{v}_clones.pdf", v=variants.keys())

rule genotypes:
    input:
        expand("figures/{v}_counts.pdf", v=variants.keys())

rule rate_table:
    input:
        rate_files = expand("rates/{v}_rate.json", v=variants.keys()),
        gt = "data/clade_gts.json"
    output:
        rate_table = "data/rates.tsv",
        rate_table_tex = "manuscript/rates.tex"
    run:
        import json
        import pandas as pd
        from treetime.utils import datestring_from_numeric
        with open(input.gt) as fin:
            clade_gts = json.load(fin)

        offset_nuc = lambda clade: -2 if clade=='19B' else 2
        offset_nonsyn = lambda clade: -1 if clade=='19B' else 1

        data = []
        for fname in input.rate_files:
            with open(fname) as fin:
                d = json.load(fin)
            base_clade = d['clade'][:3]
            aa_div = len([x for x in clade_gts[base_clade]['aa'] if 'ORF9' not in x]) + offset_nonsyn(base_clade)
            nuc_div = len(clade_gts[base_clade[:3]]['nuc']) + offset_nuc(base_clade)
            data.append({'clade':d['clade'], 'nuc_rate': d['nuc']['slope'],
                         'nuc_origin': d['nuc']['origin'], 'nuc_origin_date': datestring_from_numeric(d['nuc']['origin']),
                         'aa_rate': d['aa']['slope'],
                         'aa_origin':d['aa']['origin'], 'aa_origin_date':datestring_from_numeric(d['aa']['origin']),
                         'syn_rate': d['syn']['slope'],
                         'syn_origin':d['syn']['origin'],'syn_origin_date':datestring_from_numeric(d['syn']['origin']),
                         'spike_rate': d['spike']['slope'],
                         'orf1ab_rate': d['orf1']['slope'],
                         'nuc_div': nuc_div, 'aa_div':aa_div, 'syn_div':nuc_div-aa_div})

        df = pd.DataFrame(data)
        df.to_csv(output.rate_table, sep='\t')
        df.to_latex(output.rate_table_tex, float_format="%.2f", index=False)


rule rate_summary:
    input:
        rate_table = "data/rates.tsv"
    output:
        figure = "figures/rate_summary.pdf",
        figure_rates = "figures/rate_progression.pdf",
        figure_rates_genes = "figures/rate_progression_by_gene.pdf"
    shell:
        """
        python3 scripts/combine_fits.py --rate-table {input.rate_table}\
                                       --output-plot {output.figure} \
                                       --output-plot-rates {output.figure_rates} \
                                       --output-plot-rates-genes {output.figure_rates_genes}
        """


rule fitness_costs:
    input:
        ref = 'data/reference.gb',
        nextclade = 'data/nextclade.tsv.gz',
        pango_gts = 'data/pango_gts.json'
    output:
        fitness_costs = "data/fitness.tsv",
        mutation_rates = "data/mutation_rates.tsv",
        all_events = "data/mutation_events.tsv"
    shell:
        """
        python3 scripts/count_mutations.py --metadata {input.nextclade}\
                 --reference {input.ref} \
                 --pango-gts {input.pango_gts} \
                 --output-fitness {output.fitness_costs} \
                 --output-events {output.all_events} \
                 --output-mutations {output.mutation_rates}
        """

rule fitness_figures:
    input:
        fitness_costs = "data/fitness.tsv",
        mutation_rates = "data/mutation_rates.tsv"
    output:
        fitness_figure = "figures/fitness_cost.pdf",
        fitness_figure_by_gene = "figures/fitness_cost_by_gene.pdf",
        mutation_figure = "figures/mutation_distribution.pdf"
    shell:
        """
        python3 scripts/plot_fitness.py --fitness {input.fitness_costs}\
                 --mutations {input.mutation_rates} \
                 --output-fitness {output.fitness_figure} \
                 --output-fitness-by-gene {output.fitness_figure_by_gene} \
                 --output-mutations {output.mutation_figure}
        """

rule fitness_landscape:
    input:
        fitness_costs = "data/fitness.tsv",
    output:
        fitness_landscape = "figures/fitness_landscape.pdf",
    shell:
        """
        python3 scripts/plot_fitness_landscape.py --fitness {input.fitness_costs}\
                 --output-fitness-landscape {output.fitness_landscape}
        """
