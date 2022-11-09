from pprint import pprint
from typing import Generator
from vcf2circos.plotcategories.plotconfig import Plotconfig
from vcf2circos.utils import timeit

from os.path import join as osj
import pandas as pd
from itertools import chain, repeat
from collections import OrderedDict, Counter


class Histogram_(Plotconfig):
    """
    It need to create one histogram for each SV event FOR EACH SV height (from 0 copy number to 5), which will create the grey band between color dor
    """

    def __init__(self, plotconfig):
        self.plotconfig = plotconfig
        self.variants_position = self.config_ring["position"]
        self.variants_ring_space = self.config_ring["space"]
        self.variants_ring_height = self.config_ring["height"]
        # corresponding to SNV InDel height 7th ring (after 0 to 5 copy number height)
        self.radius = {
            "R0": self.variants_position
            + (max(self.rangescale) * self.variants_ring_space)
            + ((max(self.rangescale) + 1) * self.variants_ring_height),
            "R1": self.variants_position
            + (max(self.rangescale) * self.variants_ring_space)
            + ((max(self.rangescale) + 2) * self.variants_ring_height),
        }
        print("#Range", self.rangescale)
        self.hovertextformat = ' "<b>{}:{}-{}</b><br>{}<br><br>{}".format(a[i,0], a[i,1], a[i,2], a[i,6], a[i,8])'
        self.trace = {
            "hoverinfo": "text",
            "mode": "markers",
            "marker": {"size": 5, "symbol": 0, "color": "gray", "opacity": 0.1},
        }
        self.layout = {
            "type": "path",
            "opacity": 1,
            "fillcolor": "gray",
            "line": {"color": "gray", "width": 5},
        }
        # TODO tile same as cytobandinfo in vfreader
        self.cytoband_conf = pd.read_csv(
            osj(
                self.options["Static"],
                "Assembly",
                self.options["Assembly"],
                "cytoband_" + self.options["Assembly"] + "_chr_infos.txt.gz",
            ),
            sep="\t",
            header=0,
            compression="infer",
        )
        self.cytoband_data = {
            "show": "True",
            "file": {
                "path": "",
                "header": "infer",
                "sep": "\t",
                "dataframe": {"orient": "columns", "data": None},
            },
            "colorcolumn": 4,
            "radius": {"R0": 1, "R1": 1.1},
            "hovertextformat": " \"<b>{}:{}-{}<br>{}{}</b>\".format(a[i,0], a[i,1], a[i,2], a[i,0].replace('chr', ''), ''.join(a[i,5:]))",
            # "hovertextformat": " \"<b>{}</b>\".format(a[i,0])",
            "trace": {
                "uid": "cytoband_tile",
                "hoverinfo": "text",
                "mode": "markers",
                "marker": {
                    "size": 0,
                    "symbol": 0,
                    "color": None,
                    "opacity": 0,
                },  # 8
            },
            "layout": {
                "type": "path",
                "layer": "above",
                "opacity": 0,
                "line": {"color": None, "width": 0},
            },
        }
        self.genes = pd.read_csv(
            osj(
                self.options["Static"],
                "Assembly",
                self.options["Assembly"],
                "genes." + self.options["Assembly"] + ".sorted.txt",
            ),
            header=0,
            sep="\t",
        ).drop_duplicates(subset="gene", keep="first")
        self.df_data = pd.DataFrame.from_dict(self.data).astype(
            {
                "Chromosomes": str,
                "Genes": str,
                "Exons": str,
                "Variants": object,
                "Variants_type": str,
                "CopyNumber": int,
                "Color": str,
            }
        )
        self.df_morbid = pd.read_csv(
            osj(self.options["Static"], "morbid.txt"),
            header=None,
            sep="\t",
            names=["genes"],
        )

    def __getattr__(self, item):
        if hasattr(self.plotconfig, item):
            return getattr(self.plotconfig, item)

    def dict_to_str(self, info_field: list) -> Generator:
        for info_dict in info_field:
            yield ";".join(
                [str(key) + "=" + str(value) for key, value in info_dict.items()]
            )

    def adapt_data(self, cn: int) -> dict:
        d_file = {
            "path": "",
            "header": "infer",
            "sep": "\t",
            "dataframe": {"orient": "columns", "data": None},
        }
        data = {
            "chr_name": [],
            "start": [],
            "end": [],
            "val": [],
            "ref": [],
            "alt": [],
            "type": [],
            "color": [],
            "hovertext": [],
            "symbol": [],
            "genes": [],
            "exons": [],
        }

        df_data = self.df_data.loc[self.df_data["CopyNumber"] == cn]
        start = []
        stop = []
        ref = []
        alt = []
        for items in list(
            self.extract_start_stop_ref_alt(
                df_data["Record"].to_list(),
                df_data["Variants"].to_list(),
                df_data["Variants_type"].to_list(),
            )
        ):
            # DEBUGG
            # print(*items)
            # for val in items:
            #    if isinstance(val, list):
            #        print(type(val[0]))
            #    else:
            #        print(type(val))
            # exit()
            start.append(items[0])
            stop.append(items[1])
            ref.append(items[2])
            alt.append(str(items[3][0]))
        data["chr_name"].extend(df_data["Chromosomes"].to_list())
        data["start"].extend(start)
        data["end"].extend(stop)
        data["val"].extend(list(repeat(2, len(df_data.index))))
        data["ref"].extend(ref)
        data["alt"].extend(alt)
        data["type"].extend(df_data["Variants_type"].to_list())
        data["color"].extend(list(repeat("grey", len(df_data.index))))
        # data["hovertext"].extend(list(itertools.repeat("", len(df_data.index))))
        data["hovertext"].extend(list(self.generate_hovertext_var(df_data["Variants"])))
        # data["hovertext"].extend(
        #    [
        #        "Genes ("
        #        + str(len(record.split(",")))
        #        + "): "
        #        + ",".join(record.split(",")[:5])
        #        for record in df_data["Genes"].to_list()
        #    ]
        # )
        data["symbol"].extend(list(repeat(0, len(df_data.index))))
        data["genes"].extend(df_data["Genes"].to_list())
        data["exons"].extend(list(repeat("", len(df_data.index))))
        # data["info"].extend(list(self.dict_to_str(df_data["Variants"].to_list())))
        d_file["dataframe"]["data"] = data
        return d_file

    def histo_cnv_level(self, cn: int) -> dict:
        d = {}
        d_file = self.adapt_data(cn)
        d["show"] = "True"
        d["customfillcolor"] = "False"
        d["file"] = d_file
        d["sortbycolor"] = "False"
        d["colorcolumn"] = 7
        radius = (
            self.rangescale[cn]
            + self.rangescale[cn]
            + self.options["Variants"]["rings"]["height"]
        ) / 2
        d["radius"] = {
            "R0": radius,
            "R1": radius,
        }
        d["hovertextformat"] = self.hovertextformat
        d["trace"] = {
            "hoverinfo": "text",
            "mode": "markers",
            "marker": {
                "size": 5,
                "symbol": d_file["dataframe"]["data"]["symbol"],
                "color": "gray",
                "opacity": 0.1,
            },
            "uid": "cnv_scatter_level_" + str(cn),
        }
        d["layout"] = {
            "type": "path",
            "layer": "above",
            "opacity": 0.1,
            "fillcolor": "red",
            "line": {"color": "lightgray", "width": 5},
        }
        return d

    def generate_hovertext_var(self, variants_list) -> Generator:
        # print(self.data["Variants"])
        # print(len(self.data["Variants"]))
        # print(len(self.data["Chromosomes"]))
        # exit()
        # dict containing INFO field for each var
        for var in variants_list:
            yield "<br>".join(
                [
                    ": ".join(
                        [
                            str(value) if not isinstance(value, list) else str(value[0])
                            for value in pairs
                        ]
                    )
                    for pairs in list(zip(var.keys(), var.values()))
                ]
            )

    def merge_options(self, cytoband_data: dict) -> list:
        """
        func handle math and geometry need to take data in specific order
        chr start stop val OTHERWIS TROUBLE
        """
        cyto = {}
        cyto["chr_name"] = cytoband_data["chr_name"]
        cyto["start"] = cytoband_data["start"]
        cyto["end"] = cytoband_data["end"]
        # Remember to have val column in data otherwise it leads to crash]
        cyto["val"] = list(repeat(1, len(cytoband_data["chr_name"])))
        cyto["band_color"] = list(repeat("lightgray", len(cytoband_data["chr_name"])))
        cyto["band"] = cytoband_data["band"]
        # Cytoband tiles 3  need fill data
        self.cytoband_data["file"]["dataframe"]["data"] = cyto

        self.cytoband_data["layout"]["line"]["color"] = cyto["band_color"]
        self.cytoband_data["trace"]["marker"]["color"] = cyto["band_color"]
        whole_cn = []
        # Histo_cnv_level
        for cn in list(set(self.data["CopyNumber"])):
            res = self.histo_cnv_level(cn)
            whole_cn.append(res)
        # Genes plots
        whole_cn.append(self.histo_genes())

        # cytoband tiles
        whole_cn.append(self.cytoband_data)
        return whole_cn

        # def __call__(self):
        #    return pd.DataFrame.from_dict(self.data)

    def process_gene_list(self, genes_list: list) -> Generator:
        for record in genes_list:
            if record:
                yield record.split(",")

    def morbid_genes(self, genes: list) -> Generator:
        for g in genes:
            if g in self.df_morbid["genes"].to_list():
                yield "red"
            else:
                yield "lightgray"

    @timeit
    def histo_genes(self) -> dict:
        data = {}
        dico = {}
        # remove empty gene, df_data attribute of class basic data from plot config Parents class
        # gene_list = list(filter(lambda x: x != "", self.df_data["Genes"]))
        gene_list = list(
            set(
                list(
                    map(
                        str,
                        chain.from_iterable(
                            list(self.process_gene_list(self.df_data["Genes"]))
                        ),
                    )
                )
            )
        )
        print(*self.genes.columns)
        print(self.genes.head())
        # select genes in or batch of variations (from refeseq assembly)
        df_filter = self.genes.loc[self.genes["gene"].isin(gene_list)]
        print(df_filter.head())
        print(*df_filter.columns)
        print(*gene_list)
        for fields in df_filter.columns:
            if fields != "transcript":
                if fields == "color":
                    data[fields] = list(self.morbid_genes(df_filter["gene"]))
                else:
                    data[fields] = df_filter[fields].to_list()
        # pprint(data, sort_dicts=False)
        dico["file"] = {
            "path": "",
            "header": "infer",
            "sep": "\t",
            "dataframe": {"orient": "columns", "data": data},
        }
        dico["show"] = self.show
        dico["colorcolumn"] = 4
        dico["radius"] = {"R0": 0.98, "R1": 0.98}
        dico[
            "hovertextformat"
        ] = ' "<b>{}:{}-{}<br>Gene: {}</b><br>".format(a[i,0], a[i,1], a[i,2], a[i,5])'
        dico["trace"] = {
            "uid": "genes",
            "hoverinfo": "text",
            "mode": "markers",
            "marker": {
                "size": 3,
                "symbol": 0,
                "color": data["color"],
                "opacity": 1,
            },
        }
        dico["layout"] = {
            "type": "path",
            "layer": "above",
            "opacity": 0.2,
            "line": {"color": data["color"], "width": 3},
        }
        return dico

    def genes_omim_morbid(self):
        """ "
        If it's a morbid gene it will be colored in red in circos gene level
        done in static file in genes.<assembly>
        """
        pass

    def extract_start_stop_ref_alt(
        self, record: list, info_field: list, variant_type: list
    ) -> Generator:
        # infer type of var could be done before
        for i, info_dict in enumerate(info_field):
            if variant_type[i] not in ["OTHER", "SNV", "INDEL"]:
                if "SV_start" in record[i].INFO and "SV_end" in record[i].INFO:
                    yield (int(info_dict.get("SV_start")), int(info_dict.get("SV_end")))
                elif "END" in record[i].INFO:
                    yield (
                        int(record[i].POS),
                        int(record[i].INFO["END"]),
                        record[i].REF,
                        record[i].ALT,
                    )
                elif "SVLEN" in record[i].INFO:
                    yield (
                        int(record[i].POS),
                        int(abs(record[i].INFO["SVLEN"][0])) + int(record[i].POS),
                        record[i].REF,
                        record[i].ALT,
                    )
                else:
                    print("Can't establish SV length, annotations missing EXIT")
                    exit()
            # SNVINDEL
            else:
                alternate = int(str(max([len(alt) for alt in record[i].ALT])))
                yield (
                    int(str(record[i].POS)),
                    int(str(record[i].POS)) + alternate,
                    record[i].REF,
                    record[i].ALT,
                )
