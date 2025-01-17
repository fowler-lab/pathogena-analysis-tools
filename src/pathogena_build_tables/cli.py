import pandas, pathlib, json

from tqdm import tqdm

from collections import defaultdict


def split_species(row):
    cols = row["name"].split("(")
    species = cols[0]
    lineage = cols[1][:-1]
    return pandas.Series([species, lineage])


def build(
    master_table: str = None,
    output_path: str = "data/",
    tables_path: str = None,
):
    master_file = pathlib.Path(master_table)
    master_table = pandas.read_csv(master_file)
    master_table["has_effects"] = False
    master_table["has_variants"] = False
    master_table["has_mutations"] = False
    master_table["has_predictions"] = False
    master_table["has_main_report"] = False
    master_table["in_mapping_file"] = False
    master_table.set_index("UNIQUEID", inplace=True)

    path = pathlib.Path(output_path)

    tables = defaultdict(list)

    for folder in ["effects", "mutations", "predictions", "variants"]:

        for i in tqdm((path / folder).glob("*.csv")):

            df = pandas.read_csv(i)
            uid = i.stem.split("." + folder)[0]
            df["uniqueid"] = uid
            master_table.at[uid, "has_" + folder] = True
            tables[folder].append(df)

    effects = pandas.concat(tables["effects"])
    mutations = pandas.concat(tables["mutations"])
    predictions = pandas.concat(tables["predictions"])
    variants = pandas.concat(tables["variants"])

    for folder in ["main_report"]:

        for i in tqdm((path / folder).glob("*.json")):

            uid = i.stem.split("." + folder)[0]
            master_table.at[uid, "has_" + folder] = True

            row = [uid]
            f = open(i)
            data = json.load(f)

            if (
                data["Pipeline Outcome"]
                != "Sufficient reads mapped to M. tuberculosis (H37Rv v3) for genome assembly, resistance prediction and relatedness assessment."
            ):
                continue

            row.append(data["Organism Identification"]["Mycobacterium Reads"])
            row.append(data["Mycobacterium Results"]["Summary"][0]["Name"])
            row.append(data["Mycobacterium Results"]["Summary"][0]["Num Reads"])
            row.append(data["Mycobacterium Results"]["Summary"][0]["Coverage"])
            row.append(data["Mycobacterium Results"]["Summary"][0]["Depth"])

            antibiogram = ""
            amr_results = data["Genomes"][0]["Resistance Prediction"][
                "Resistance Prediction Summary"
            ]
            for j in amr_results:
                for k in amr_results[j]:
                    antibiogram += amr_results[j][k]
            row.append(antibiogram)

            row.append(data["Metadata"]["Pipeline build"])

            tables["main_report"].append(row)

    genomes = pandas.DataFrame(
        tables["main_report"],
        columns=[
            "uniqueid",
            "mycobacterial_reads",
            "name",
            "tb_reads",
            "tb_coverage",
            "tb_depth",
            "antibiogram",
            "pipeline_build",
        ],
    )
    genomes[["species", "lineage"]] = genomes.apply(split_species, axis=1)
    genomes.drop(columns=["name"], inplace=True)
    genomes = genomes[
        [
            "uniqueid",
            "species",
            "lineage",
            "mycobacterial_reads",
            "tb_reads",
            "tb_coverage",
            "tb_depth",
            "antibiogram",
            "pipeline_build",
        ]
    ]

    tables_path = pathlib.Path(tables_path)

    genomes.to_csv(tables_path / "GENOMES.csv")
    variants.to_csv(tables_path / "VARIANTS.csv")
    mutations.to_csv(tables_path / "MUTATIONS.csv")
    effects.to_csv(tables_path / "EFFECTS.csv")
    predictions.to_csv(tables_path / "PREDICTIONS.csv")

    genomes.to_parquet(tables_path / "GENOMES.parquet")
    variants.to_parquet(tables_path / "VARIANTS.parquet")
    mutations.to_parquet(tables_path / "MUTATIONS.parquet")
    effects.to_parquet(tables_path / "EFFECTS.parquet")
    predictions.to_parquet(tables_path / "PREDICTIONS.parquet")


def main():
    import defopt

    defopt.run(
        [build],
        no_negated_flags=True,
        strict_kwonly=False,
        short={},
    )


if __name__ == "__main__":
    main()
