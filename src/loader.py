import gzip
import pandas as pd
import os
from collections import defaultdict
from combat.pycombat import pycombat


def load_ann(path: str, keep_symbols: set, filtering: bool = True) -> tuple:
    id_name_map = {}
    # id_len_map = {}
    name_len_map = defaultdict(list)
    name_id_map = defaultdict(list)
    with gzip.open(path, 'r') as inf:
        for line in inf:
            split_line = line.strip().split(b"\t")
            if len(split_line) < 9:
                continue
            if split_line[2] == b"gene":
                gene_info = split_line[8].strip(b";").split(b"; ")
                # assert gene_info[3].startswith(b"gene_name \"")
                gene_name = str(gene_info[3][11:len(gene_info[3]) - 1], "ascii")
                if not filtering or gene_name in keep_symbols:
                    gene_id = str(gene_info[0][9:len(gene_info[0]) - 1], "ascii")
                    gene_len = int(split_line[4]) - int(split_line[3])
                    name_len_map[gene_name].append(gene_len)
                    name_id_map[gene_name].append(gene_id)
                    id_name_map[gene_id] = gene_name
                    # id_len_map[gene_id] = gene_len

    # Ne zelimo
    # dropping_names = filter(lambda x: any(it != name_len_map[x][0] for it in name_len_map[x]), name_len_map)
    dropping_names = filter(lambda x: len(name_len_map[x]) > 1, name_len_map)

    for it in list(dropping_names):
        for gene_id in name_id_map[it]:
            # del id_len_map[gene_id]
            del id_name_map[gene_id]
        del name_id_map[it]
        del name_len_map[it]

    name_len_df = pd.DataFrame.from_dict(
        {it: jt[0] for it, jt in name_len_map.items()},
        orient="index", columns=["length"]
    )

    return id_name_map, name_len_df


def load_tcga(path: str, id_name_map: dict, length_df: pd.DataFrame) -> pd.DataFrame:
    patients = pd.DataFrame(index=length_df.index.copy())
    # print(sorted(patients.index))
    # print("# of patients", len(patients.index))
    for patient in os.listdir(path):
        new_path = os.path.join(path, patient)
        if not os.path.isdir(new_path):
            continue
        patient_fn = next(filter(
            lambda x: x.endswith(".gz"),
            os.listdir(new_path)
        ))
        patient_path = os.path.join(new_path, patient_fn)
        df = pd.read_csv(
            patient_path, sep='\t', names=("gene_id", "raw_read_count")
        )
        remove_nameless_ids(df, id_name_map)
        df.dropna(axis=1, inplace=True)
        add_name_index(df, id_name_map)
        df.drop(columns="gene_id", inplace=True)
        tpm_normalize(df, length_df)
        patients[patient] = df["tpm"]
    return patients


def load_icgc(path: str, id_name_map: dict, length_df: pd.DataFrame, keep_symbols: set) -> pd.DataFrame:
    patients = pd.DataFrame(index=length_df.index.copy())
    keep_columns = [
        "icgc_sample_id",
        "gene_id",
        "raw_read_count"
    ]
    df = pd.read_csv(
        path,
        sep="\t",
        #  index_col=["icgc_sample_id", "gene_id"],
        usecols=keep_columns
    )
    df.rename(columns={"gene_id": "gene_name"}, inplace=True)
    # df.drop_duplicates(inplace=True, ignore_index=True)
    keep_symbols = set(length_df.index)
    filter_by_gene_name(df, keep_symbols)
    df.dropna(axis=1, inplace=True)
    dfs = [it for it in df.groupby(["icgc_sample_id"])]
    # pd.set_option('display.max_columns', None)
    print("First:", df.isnull().values.any())
    for i, (patient, df) in enumerate(dfs):
        # for gene_name, df2 in df.groupby(["gene_name"]):
        #     print(df2.head())
        # exit()
        df.set_index("gene_name", inplace=True)
        tpm_normalize(df, length_df)
        if df.isnull().values.any():
            print("patient")
        patients[patient] = df["tpm"]
    return patients


def remove_nameless_ids(df: pd.DataFrame, id_name_map: dict) -> None:
    df.drop(
        df[[it not in id_name_map for it in df["gene_id"]]].index, inplace=True
    )


def add_name_index(df: pd.DataFrame, id_name_map: dict) -> None:
    df["gene_name"] = df.apply(lambda x: id_name_map[x["gene_id"]], axis=1)
    df.set_index("gene_name", inplace=True)


def remove_id_column(df: pd.DataFrame) -> None:
    df.drop(columns="gene_id", inplace=True)


def tpm_normalize(df: pd.DataFrame, lengths_df: pd.DataFrame):
    rpk_df = df["raw_read_count"] / lengths_df["length"]
    scaling_factor = rpk_df.sum() / 1000000
    df["tpm"] = rpk_df * scaling_factor
    df.drop(columns="raw_read_count", inplace=True)
# TODO je li filtering usklađen sa anotacijama remove nameless ids


def filter_by_gene_name(df: pd.DataFrame, keep_symbols: set):
    print("Filtering by gene name")
    df.drop(
        df[[it not in keep_symbols for it in df["gene_name"]]].index, inplace=True
    )
    print("Done")


def load_metabolyc_genes(path: str) -> set:
    return set(
        pd.read_excel(path, sheet_name="All metabolic genes")["Gene Symbol"]
    )
    # zasto je SLC25A25 dupliciran?


if __name__ == "__main__":
    ANNOTATIONS_FILE = "../res/gencode.v22.annotation.gtf.gz"
    METABOLYC_FILE = "../res/41586_2011_BFnature10350_MOESM321_ESM.xls"
    TCGA_FOLDER = "../res/gdc_download_20211208_083143.349897"
    ICGC_FILE = "../res/icgc-dataset-1638970984337/exp_seq.tsv.gz"
    TCGA_CLEANED = "../res/tcga.csv"
    ICGC_CLEANED = "../res/icgc.csv"

    if os.path.isfile(TCGA_CLEANED) and os.path.isfile(ICGC_CLEANED):
        tcga = pd.read_csv(TCGA_CLEANED, index_col=0)
        icgc = pd.read_csv(ICGC_CLEANED, index_col=0)
    else:
        mg = load_metabolyc_genes(METABOLYC_FILE)
        inm, ld = load_ann(ANNOTATIONS_FILE, mg)

        tcga = load_tcga(TCGA_FOLDER, inm, ld)
        icgc = load_icgc(ICGC_FILE, inm, ld, mg)
        print(set(icgc.index) - set(tcga.index))
        tcga.to_csv(TCGA_CLEANED)
        icgc.to_csv(ICGC_CLEANED)

    tcga.dropna(axis=1, inplace=True)
    icgc.dropna(axis=1, inplace=True)
    # total = pd.concat([tcga, icgc], join="inner", axis=1)
    # batch = [0] * len(tcga.columns) + [1] * len(icgc.columns)
    # data_corrected = pycombat(total, batch)
