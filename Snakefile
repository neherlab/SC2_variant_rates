variants = ['19A', '19B', '20A', '20B', '20C', '20E', '20H', '20I', '20J']


rule get_data:
    output:
        "metadata.tsv.gz"
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
        metadata = "metadata.tsv.gz"
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
