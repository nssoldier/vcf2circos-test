import json
import subprocess
import re
import pandas as pd
from os.path import join as osj
import os
import gzip
from tqdm import tqdm
import vcf2circos

# from natsort import natsort_keygen
from typing import Generator
from functools import wraps
import time
import numpy as np
import pyfiglet
from pprint import pprint

# Globals
# variants_color = {
#    "INS": "red",
#    "INV": "purple",
#    "DEL": "orange",
#    "DUP": "blue",
#    "CNV": "brown",
#    "BND": "blue",
#    "SNV": "dimgray",
#    "INDEL": "dimgray",
#    "OTHER": "dimgray",
# }


# variants_color = {
#    "INS": "#A52A2A",
#    "INV": "#9B30FF",
#    "DEL": "#FF8000",
#    "DUP": "#0000FF",
#    "CNV": "#A52A2A",
#    "BND": "#0000FF",
#    "SNV": "#808080",
#    "INDEL": "#808080",
#    "OTHER": "#808080",
# }
def get_swap_dict(d):
    """
    https://note.nkmk.me/en/python-dict-swap-key-value/
    """
    return {v: k for k, v in d.items()}


def chr_valid():
    chr_valid = ["chr" + str(i) for i in range(1, 23)]
    chr_valid.extend(["chrM", "chrX", "chrY"])
    return chr_valid


def launch():
    """
    https://ascii.co.uk/art/dna
    """
    print(
        """
-._    _.--'"`'--._    _.--'"`'--._    _.--'"`'--._    _   
    '-:`.'|`|"':-.  '-:`.'|`|"':-.  '-:`.'|`|"':-.  '.` : '.   
  '.  '.  | |  | |'.  '.  | |  | |'.  '.  | |  | |'.  '.:   '.  '.
  : '.  '.| |  | |  '.  '.| |  | |  '.  '.| |  | |  '.  '.  : '.  `.
  '   '.  `.:_ | :_.' '.  `.:_ | :_.' '.  `.:_ | :_.' '.  `.'   `.
         `-..,..-'       `-..,..-'       `-..,..-'       `         `
    """
    )
    print(pyfiglet.figlet_format("         vcf2circos", font="big"))
    print("Author: " + vcf2circos.__author__)
    print("Version: " + vcf2circos.__version__)
    print("Last update: " + vcf2circos.__date__)
    print("\n")


class Colorpal:
    def __init__(self, n=None):
        self.color = self.random_rgb()
        if n is not None:
            self.colorpal = [self.random_rgb() for c in range(0, n)]

    # def __call__(self):
    #    if len(self.colorpal) == 0:
    #        print("Specify length of color pal EXIT")
    #        exit()
    #    else:
    #        return self.colorpal

    def __iter__(self):
        try:
            for v in self.colorpal:
                yield v
        except AttributeError:
            print("Specify length of color pal EXIT")
            exit()

    def random_rgb(self):
        return tuple(map(str, list(np.random.choice(range(255), size=3))))


# Utils func
def json_to_dict(jsonpath: str) -> dict:
    """
    Load json file and create a dict
    """
    with open(jsonpath) as json_file:
        return json.load(json_file)


def check_data_plot(dico, list_keys=None):
    try:
        var_numb = len(dico["chr_name"])
    except KeyError:
        print(dico.keys())
        exit()
    for d_fields in dico:
        assert len(dico[d_fields]) == var_numb, (
            "Missing values in "
            + d_fields
            + " fields\n\t...chrom: "
            + str(var_numb)
            + " "
            + d_fields
            + ": "
            + str(len(dico[d_fields]))
        )
    if "end" in dico:
        pass
    # else:
    #    if list_keys is None:
    #        list_keys = ["chr_name", "start", "val", "ref", "alt", "type", "color"]
    #    assert list_keys == list(dico.keys())[:7], (
    #        "ISSUES wrong list data order, leads to crash in mathematical operations \n\t..."
    #        + ", ".join(list(dico.keys())[:7])
    #        + " ,..."
    #    )


def delete_multiple_element(list_object, indices):
    """
    from https://thispointer.com/python-remove-elements-from-list-by-index/
    """
    indices = sorted(indices, reverse=True)
    for idx in indices:
        if idx < len(list_object):
            list_object.pop(idx)


def map_annotations(field_annot):
    if field_annot is not None:
        uniq = list(set(str(field_annot).split("|")))
        # uniq = list(set(str(field_annot)))
        # if len(uniq) == 1:
        #    return uniq[0]
        # else:
        #    return field_annot
        return uniq[0]
    else:
        return "."


def generate_hovertext_var(variants_list, full_annot=None, true_annot=None) -> Generator:
    # 30 longueur char
    # 15 hauteur annot
    def crop_annotations(iterable):
        # list(map(str, pairs[1]))[:-1]
        iterable_string = list(map(str, iterable))
        if len(iterable_string) == 1:
            return iterable_string
        elif len(iterable_string) == 2:
            return [iterable_string[1]]
        else:
            return iterable_string[:-1]

    for var in variants_list:
        tmp = []
        for i, pairs in enumerate(list(zip(var.keys(), var.values()))):
            # if pairs[0] == "OMIM_phenotype":
            # print(pairs[1])
            if true_annot:
                # If user want this annotations
                if pairs[0] not in true_annot:
                    # print(pairs[0] + " not in")
                    continue

            if full_annot is not None:
                if i == full_annot:
                    break
            else:
                if i == 15:
                    break
            if not isinstance(pairs[1], list):
                tmp.append(": ".join([pairs[0], str(pairs[1])]))
            else:
                tmp.append(": ".join([pairs[0], ",".join(crop_annotations(pairs[1]))]))
        to_add = []
        for items in tmp:
            if len(items) > 40:
                items = "".join(items[:40]) + "..."
                to_add.append(items)
            else:
                to_add.append(items)
        # exit()
        yield "<br>".join(to_add)
        # exit()
        # "SV_chrom",


# "SV_start",
# "SV_end",
# "FORMAT",
# "SpliceAI",
# "SPiP",
# "ACMG",
# "varankVarScore",
# "gene",
# "zygosity",
# "rsClinicalSignificance",
# "OMIM_ID",
# "OMIM_inheritance",


def cast_svtype(svtype):
    """
    Handle "<" in vcf
    """
    if isinstance(svtype, list):
        svtype = svtype[0]
    if "<" in svtype or ">" in svtype:
        svtype = svtype.replace("<", "")
        svtype = svtype.replace(">", "")
    if ":" in svtype:
        return svtype.split(":")[0]
    else:
        return svtype


def systemcall(command: str) -> list:
    """
    Call bash command
    https://github.com/JbaptisteLam/DPNI/blob/main/src/utils/utils.py

    In case of crash exit code 1 stop script
    """
    print("#[SYS] " + command)
    p = subprocess.Popen([command], stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    out, err = p.communicate()
    if not err:
        return out.decode("utf8").strip().split("\n")
    else:
        issues = err.decode("utf8").strip()
        try:
            re.search(r"(Warning|WARNING)", issues).group()
            print("--WARNING Systemcall--\n", err.decode("utf8").strip())
            return out.decode("utf8").strip().split("\n")
        except AttributeError:
            print("--ERROR Systemcall--\n", err.decode("utf8").strip())
            exit()


def formatted_refgene(refgene: str, assembly: str, ts=None) -> str:
    """
    In case of new version of refgene or new assembly\n
    transcripts: list in string format either 'NM_' or 'NR_ 'or both 'NM_,NR_'
    Took refgene raw file from ucsc curated and create proper exon refgene, WITHOUT UTR(default choice)
    https://genome.ucsc.edu/cgi-bin/hgTrackUi?g=refSeqComposite&db=hg38
    ex: python -c 'from foo import hello; print hello()'
     zcat transcripts.hg38.txt.gz.tmp | grep chr_name > transcripts.hg38.sorted.txt && zcat transcripts.hg38.txt.gz.tmp | grep -v "chr_name" | sort -k1,1V -k2,2n >> transcripts.hg38.sorted.txt
     same for both 3 files except exons sort after chr and pos by exons
     zcat exons.hg38.txt.gz.tmp | grep chr_name > exons.hg38.sorted.txt && zcat exons.hg38.txt.gz.tmp | grep -v "chr_name" | sort -k1,1V -k2,2n >> exons.hg38.sorted.txt
    """
    if ts is not None:
        transcripts = ts.split(",")
    else:
        transcripts = "NM_"
    cols = [
        "bin",
        "name",
        "chrom",
        "strand",
        "txStart",
        "txEnd",
        "cdsStart",
        "cdsEnd",
        "exonCount",
        "exonStarts",
        "exonEnds",
        "score",
        "name2",
        "cdsStartStat",
        "cdsEndStat",
        "exonFrames",
    ]
    df = pd.read_csv(refgene, sep="\t", header=None, compression="infer")
    assert (
        len(df.columns) == 16
    ), "Error in columns format more or less than 16 columns, expected\n\tcols: " + ",".join(cols)
    df.columns = cols
    output_genes = os.path.dirname(refgene) + "/genes." + assembly + ".txt"
    output_transcripts = os.path.dirname(refgene) + "/transcripts." + assembly + ".txt.gz"
    output_exons = os.path.dirname(refgene) + "/exons." + assembly + ".txt.gz"

    # og = os.path.dirname(refgene) + "/genes." + assembly + ".txt.gz"
    # ot = os.path.dirname(refgene) + "/transcripts." + assembly + ".txt.gz"
    # oe = os.path.dirname(refgene) + "/exons." + assembly + ".txt.gz"

    with gzip.open(output_exons, "wb+") as out_e:
        with gzip.open(output_transcripts, "wb+") as out_t:
            out_t.write(
                bytes(
                    "\t".join(
                        [
                            "chr_name",
                            "start",
                            "end",
                            "val",
                            "color",
                            "gene",
                            "transcript",
                        ]
                    )
                    + "\n",
                    "UTF-8",
                )
            )
            out_e.write(
                bytes(
                    "\t".join(
                        [
                            "chr_name",
                            "start",
                            "end",
                            "val",
                            "color",
                            "gene",
                            "exons",
                            "transcript",
                        ]
                    )
                    + "\n",
                    "UTF-8",
                )
            )
            for i, row in tqdm(
                df.iterrows(),
                total=len(df.index),
                desc="Formatting " + os.path.basename(refgene) + " file UCSC",
                leave=False,
            ):
                if [row["name"].startswith(ts) for ts in transcripts]:
                    out_t.write(
                        bytes(
                            "\t".join(
                                [
                                    row["chrom"],
                                    str(row["txStart"]),
                                    str(row["txEnd"]),
                                    "1",
                                    "lightgray",
                                    row["name2"],
                                    row["name"],
                                ]
                            )
                            + "\n",
                            "UTF-8",
                        )
                    )
                    for i in range(len(row["exonStarts"].split(",")[:-1])):
                        exons_start = row["exonStarts"].split(",")
                        exons_end = row["exonEnds"].split(",")
                        out_e.write(
                            bytes(
                                "\t".join(
                                    [
                                        row["chrom"],
                                        str(exons_start[i]),
                                        str(exons_end[i]),
                                        "1",
                                        "lightgray",
                                        row["name2"],
                                        "exon" + str(i + 1),
                                        row["name"],
                                    ]
                                )
                                + "\n",
                                "UTF-8",
                            )
                        )

    # genes, keep row with the most exon counts
    with gzip.open(output_genes + ".tmp", "wb+") as out_g:
        out_g.write(
            bytes(
                "\t".join(
                    [
                        "chr_name",
                        "start",
                        "end",
                        "val",
                        "color",
                        "gene",
                    ]
                )
                + "\n",
                "UTF-8",
            )
        )
        for gene, df_ in tqdm(
            df.groupby("name2"), total=len(df.index), desc="Formatting genes file UCSC", leave=False
        ):
            df_rows = df_.loc[df_["exonCount"] == df_["exonCount"].max()].loc[
                :, ["chrom", "cdsStart", "cdsEnd"]
            ]
            rows = df_rows.iloc[0, :]
            out_g.write(
                bytes(
                    "\t".join(
                        [
                            rows["chrom"],
                            str(rows["cdsStart"]),
                            str(rows["cdsEnd"]),
                            "1",
                            "lightgray",
                            gene,
                        ]
                    )
                    + "\n",
                    "UTF-8",
                )
            )
    systemcall(
        "zcat "
        + output_genes
        + ".tmp | grep 'chr_name' > "
        + output_genes
        + " && zcat "
        + output_genes
        + ".tmp | grep -v 'chr_name' | sort -k1,1V -k2,2n >> "
        + output_genes
        + " && /home1/TOOLS/tools/bgzip/current/bin/bgzip "
        + output_genes
    )
    if os.path.exists(output_genes + ".tmp"):
        os.remove(output_genes + ".tmp")
    print("#[INFO] Generated " + ",".join([output_genes, output_exons, output_transcripts]))
    return df


def timeit(func):
    """
    https://dev.to/kcdchennai/python-decorator-to-measure-execution-time-54hk
    """

    @wraps(func)
    def timeit_wrapper(*args, **kwargs):
        start_time = time.perf_counter()
        result = func(*args, **kwargs)
        end_time = time.perf_counter()
        total_time = end_time - start_time
        print(f"Function {func.__name__}{args} {kwargs} Took {total_time:.4f} seconds")
        return result

    return timeit_wrapper


def _gc(filedata, span, length, out):
    """
    need check in good chroms lists
    """
    dico = {}
    tmp = []
    val = []

    # 200 000
    lim = span
    threshold = span
    with open(out, "w+") as o:
        o.write("\t".join(["chr_name", "start", "end", "val", "color"]) + "\n")
        for i, lines in tqdm(
            enumerate(filedata),
            total=int(length[0]),
            desc="INFO GC percent, bins length: " + str(span) + "",
        ):
            if not isinstance(lines, str):
                lines = lines.decode("utf-8").strip()
            else:
                lines = lines.strip()
            # for each chrom
            if lines.startswith("variableStep"):
                chr = lines.split()[1].split("=")[1]
                lim = span
                tmp.clear()
                val.clear()
                continue
            else:
                # Process data
                start = int(lines.split()[0])  # 10001
                stop = int(lines.split()[0]) + 5  # 10006
                # POS
                tmp.append(start)
                # Cast to real value, 40 = 2  60 = 3 etc
                val.append(float(lines.split()[1]) / 20)
            if chr not in chr_valid():
                continue

            # if lim < start: #it Start chr 1 10001
            #    # first time only
            #    lim = start + lim
            if start >= lim or (lines.startswith("variableStep") and chr != "chr1"):
                gc_val = round(
                    (sum(val) / (stop - tmp[0])), 2
                )  # nombre de gc tout les 5 diviser par le nombres de bases tot

                # dico[chr + ":" + str(tmp[0]) + "-" + str(tmp[-1])] = gc_val
                if gc_val < 0:
                    color = "lightblue"
                else:
                    color = "blue"
                # If miss gc for example due to repeat regions
                if len(tmp) > 1:
                    o.write(
                        "\t".join(
                            [
                                str(chr),
                                str(tmp[0]),
                                str(tmp[-1]),
                                str(gc_val),
                                color,
                            ]
                        )
                        + "\n"
                    )
                tmp.clear()
                val.clear()
                lim += threshold
                # if chr == "chr3":
                #    print(dico)
                #    print(len(tmp))
                #    print(len(val))
                #    print(lim)
                #    print(lim)
                #    print(dict)
                #    exit()


def _gc_hg38(filedata, span, length, out):
    dico = {}
    tmp = []
    # 200 000
    val = []
    lim = span
    threshold = span
    chr_current = []
    df_size = pd.read_csv(
        "/home1/BAS/lamouchj/vcf2circos/Static/Assembly/hg38/GRCh38_chromFa.chromSizes",
        sep="\t",
        header=None,
    )
    df_size.columns = ["chrom", "size"]
    df_size["chrom"] = "chr" + df_size["chrom"]
    # chr_length = df_size.loc[df_size["chrom"] == chrom].loc[0].at["size"]
    chr_size = {}
    for j, k in df_size.iterrows():
        chr_size[k["chrom"]] = k["size"]
    with open(out, "w+") as o:
        o.write("\t".join(["chr_name", "start", "end", "val", "color"]) + "\n")
        for i, lines in tqdm(
            enumerate(filedata),
            total=length,
            desc="INFO GC percent, bins length: " + str(span) + "",
        ):
            if not isinstance(lines, str):
                lines = lines.decode("utf-8").strip()
            else:
                lines = lines.strip()

            chrom = lines.split()[0]

            if chrom not in dico:
                dico[chrom] = {"start": [], "stop": [], "val": []}
                lim = threshold
            if chrom not in chr_valid():
                continue

            # Process data
            start = int(lines.split()[1])
            stop = int(lines.split()[2])
            dico[chrom]["start"].append(start)
            dico[chrom]["stop"].append(stop)
            # POS
            # tmp.append(start)
            # Cast to real value, 40 = 2  60 = 3 etc
            # val.append(float(lines.split()[3]) / 20)
            dico[chrom]["val"].append(float(lines.split()[3]) / 20)
            if start >= lim:
                # print(dico[chrom]["stop"][-1])
                # print(dico[chrom]["start"][0])
                # print(sum(dico[chrom]["val"]))
                gc_val = round(
                    (sum(dico[chrom]["val"]) / (dico[chrom]["stop"][-1] - dico[chrom]["start"][0])),
                    2,
                )  # nombre de gc tout les stop-start diviser par le nombres de bases tot

                # dico[chr + ":" + str(tmp[0]) + "-" + str(tmp[-1])] = gc_val
                # if gc_val < 0:
                #    color = "lightblue"
                # else:
                #    color = "blue"
                # If miss gc for example due to repeat regions
                o.write(
                    "\t".join(
                        [
                            str(chrom),
                            str(dico[chrom]["start"][0]),
                            str(dico[chrom]["stop"][-1]),
                            str(gc_val),
                            "blue",
                        ]
                    )
                    + "\n"
                )
                dico[chrom]["start"].clear()
                dico[chrom]["stop"].clear()
                dico[chrom]["val"].clear()
                if lim >= chr_size[chrom]:
                    lim = chr_size[chrom]
                    print(lim)
                    print("End chr")
                else:
                    lim += threshold
            # chr_current.append(chrom)


def isdigit_float(string):
    if "." in string:
        numeric = string.split(".")[0]
        decimal = string.split(".")[1]
        if numeric.isdigit() and len(decimal) <= 2 and decimal.isdigit():
            return True
        else:
            return False
    else:
        return False


def full_split_annotSV():
    if len(val) > 1:
        length = len(val)
        tmp = "|".join(val)
        data[key] = [
            ",".join(tmp.split("|")[length:]),
            "|".join(tmp.split("|")[length:]),
        ]


def process_info_dict(data: dict) -> dict:
    """
    From a vcf.Record.INFO which is a dict if list contains only one string containing pipe split all values and cast string to int if values contains only digit
    \nreturn dict formatted
    """
    new_data = {}

    for key, val in data.items():
        try:
            if isinstance(val, list):
                if len(val) == 1 and isinstance(val[0], str):
                    if "|" in val[0]:
                        values = val[0].split("|")
                        new_values = []
                        for items in values:
                            if items.isdigit():
                                new_values.append(int(items))
                            elif isdigit_float(items):
                                new_values.append(float(items))
                            else:
                                new_values.append(items)
                        new_data[key] = new_values
                    elif val[0].isdigit():
                        new_data[key] = int(val[0])
                    else:
                        new_data[key] = val
                elif len(val) > 1 and "|" in val[-1]:
                    tmp = ",".join(val)
                    keep = tmp.split("|")[0].split(",")
                    keep.append("|".join(tmp.split("|")[1:]))
                    new_data[key] = keep
                else:
                    new_data[key] = val
            elif isinstance(val, str) and "|" in val:
                new_values = []
                values = val.split("|")
                for items in values:
                    if items.isdigit():
                        new_values.append(int(items))
                    elif isdigit_float(items):
                        new_values.append(float(items))
                    else:
                        new_values.append(items)
                new_data[key] = new_values
            else:
                if val.isdigit():
                    new_data[key] = int(val)
                elif isdigit_float(val):
                    new_data[key] = float(val)
                else:
                    new_data[key] = val
        except (TypeError, AttributeError):
            new_data[key] = val
    # pprint(new_data, sort_dicts=False)
    # exit()
    return new_data


def process_gc_percent(filename: str, out: str) -> str:
    """ """
    print(filename)
    print(out)
    if filename.endswith(".gz"):
        with gzip.open(filename, "rb") as bf:
            # length = systemcall("zcat " + filename + " | wc -l")
            _gc_hg38(bf, 5000000, 439385652, out)
    else:
        f = open(filename, "r")
        length = systemcall("wc -l " + filename)
        length = length[0].split()
        _gc_hg38(f, 5000000, length, out)


def pick(file):
    df = pd.read_csv(
        file,
        sep="\t",
        header=None,
        low_memory=False,
        compression="infer",
        names=["pos", "values"],
        nrows=60000000,
    )
    # filter = df.iloc[:3000000, :]
    # print(filter.iloc[:10, :])
    df.to_csv("pick_gc", header=False, index=False, sep="\t")


def repeatmasker(file, out):
    """
    Formatting repeatmasker file chr_name, start, stop, type
    """
    df = pd.read_csv(file, sep="\t", header=None, compression="infer")
    df.columns = ["chr_name", "start", "stop", "type", "length", "strand"]
    df["val"] = 0.5
    df = df.loc[
        :,
        [
            "chr_name",
            "start",
            "stop",
            "val",
        ],
    ]
    df.to_csv(out, sep="\t", header=True, index=False)
