import pandas as pd
from collections import defaultdict
import argparse,json

def genotype_struct():
    return {'nuc':{}, 'aa':defaultdict(dict)}

def assign_genotypes(n):
    gt = n["genotype"]
    if "children" in n:
        for c in n["children"]:
            cgt = genotype_struct()

            cgt["nuc"].update({k:v for k,v in gt["nuc"].items()})
            if "nuc" in c["branch_attrs"]["mutations"]:
                for mut in c["branch_attrs"]["mutations"]["nuc"]:
                    a,pos,d = mut[0], int(mut[1:-1]), mut[-1]
                    cgt["nuc"][pos] = d


            for gene in set(c["branch_attrs"]["mutations"].keys()).union(set(gt["aa"].keys())):
                if gene=='nuc':
                    continue
                cgt["aa"][gene].update({k:v for k,v in gt["aa"][gene].items()})
                if gene in c["branch_attrs"]["mutations"]:
                    for mut in c["branch_attrs"]["mutations"][gene]:
                        a,pos,d = mut[0], int(mut[1:-1]), mut[-1]
                        cgt['aa'][gene][pos] = d
            c["genotype"] = cgt
            assign_genotypes(c)

def get_pango_genotypes(n, pango_gts):
    name = n["name"]
    if "children" in n and len(n["children"]):
        for c in n["children"]:
            get_pango_genotypes(c, pango_gts)
    elif (name in ['A', 'B'] or '.' in name) and not ('/' in name):
        pango_gts[n["name"]] = n["genotype"]


if __name__=="__main__":
    parser = argparse.ArgumentParser(
        description="get variant genotypes",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument('--tree', type=str, required=True, help="input data")
    parser.add_argument('--root', type=str, required=True, help="input data")
    parser.add_argument('--output', type=str, required=True, help="output data")
    args = parser.parse_args()

    with open(args.tree) as fh:
        tree = json.load(fh)['tree']

    with open(args.root) as fh:
        root_sequence = json.load(fh)

    tree["genotype"] = genotype_struct()
    assign_genotypes(tree)
    pango_gts = {}
    get_pango_genotypes(tree, pango_gts)

    subs = {}
    for clade in pango_gts:
        tmp_nuc = []
        if 'nuc' in pango_gts[clade]:
            for pos, d in pango_gts[clade]['nuc'].items():
                a=root_sequence['nuc'][pos-1]
                if d not in ['N', '-'] and d!=a:
                    tmp_nuc.append(f"{a}{pos}{d}")

        tmp_aa = []
        for gene in pango_gts[clade]['aa']:
            for pos, d in pango_gts[clade]['aa'][gene].items():
                a=root_sequence[gene][pos-1]
                if d not in ['X', '-'] and d!=a:
                    tmp_aa.append(f"{gene}:{a}{pos}{d}")
        subs[clade] = {'nuc':tmp_nuc, 'aa':tmp_aa}

    with open(args.output, 'w') as fh:
        json.dump(subs, fh)
