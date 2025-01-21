import pandas, pathlib, json

from tqdm import tqdm

from collections import defaultdict


def split_species(row):
    cols = row["name"].split(" (")
    species = cols[0]
    # sublineage = cols[1][:-1]
    # lineage = sublineage[:9]
    return species  # pandas.Series([species, lineage, sublineage])


def build_tables(
    lookup_table: str = None,
    source_files: str = "data/",
    output: str = None,
    filename: str = "all",
):
    master_file = pathlib.Path(lookup_table)
    master_table = pandas.read_csv(master_file)
    master_table["has_effects"] = False
    master_table["has_variants"] = False
    master_table["has_mutations"] = False
    master_table["has_predictions"] = False
    master_table["has_main_report"] = False
    master_table["in_mapping_file"] = False
    master_table.set_index("UNIQUEID", inplace=True)

    path = pathlib.Path(source_files)
    tables_path = pathlib.Path(output)

    assert filename in [
        "all",
        "effects",
        "variants",
        "mutations",
        "predictions",
        "main_report",
    ], "can only specify one from this list"

    if filename in ["effects", "mutations", "predictions", "variants"]:

        tables = []

        for i in tqdm((path).rglob("*" + filename + ".csv")):

            df = pandas.read_csv(i)
            uid = i.stem.split("." + filename)[0]
            df["uniqueid"] = uid
            master_table.at[uid, "has_" + filename] = True
            tables.append(df)

        df = pandas.concat(tables)
        df.to_csv(str(tables_path / filename) + ".csv", index=False)
        # if filename != "variants":
        df.to_parquet(str(tables_path / filename) + ".parquet")

    successful_genome = 0
    too_few_reads_genome = 0
    too_few_reads_id = 0
    if filename == "main_report":

        tables = []

        for i in tqdm((path).rglob("*" + filename + ".json")):

            uid = i.stem.split("." + filename)[0]
            master_table.at[uid, "has_" + filename] = True

            row = [uid]
            f = open(i)
            data = json.load(f)

            if (
                data["Pipeline Outcome"]
                == "Mycobacterial species identified. Reads mapped to M. tuberculosis (H37Rv v3) too low to proceed to M. tuberculosis complex genome assembly."
            ):
                too_few_reads_genome += 1
                continue
            elif (
                data["Pipeline Outcome"]
                == "Number of Mycobacterial reads is too low to proceed to Mycobacterial species identification."
            ):
                too_few_reads_id += 1
                continue
            elif (
                data["Pipeline Outcome"]
                == "Sufficient reads mapped to M. tuberculosis (H37Rv v3) for genome assembly, resistance prediction and relatedness assessment."
            ):
                successful_genome += 1
            else:
                print(uid, data["Pipeline Outcome"])
                continue

            row.append(data["Organism Identification"]["Mycobacterium Reads"])
            row.append(data["Mycobacterium Results"]["Summary"][0]["Name"])
            row.append(data["Mycobacterium Results"]["Summary"][0]["Num Reads"])
            row.append(data["Mycobacterium Results"]["Summary"][0]["Coverage"])
            row.append(data["Mycobacterium Results"]["Summary"][0]["Depth"])
            lineage_results = data["Mycobacterium Results"]["Lineage"]
            n_lineages = len(lineage_results)
            if n_lineages == 1:
                if lineage_results[0]["Name"][:7] == "lineage":
                    lineage = lineage_results[0]["Name"][:8]
                    sublineage = lineage_results[0]["Name"]
                    lineage_cov = lineage_results[0]["Coverage"]
                    lineage_depth = lineage_results[0]["Median Depth"]
                else:
                    lineage = lineage_results[0]["Name"]
                    sublineage = ""
                    lineage_cov = lineage_results[0]["Coverage"]
                    lineage_depth = lineage_results[0]["Median Depth"]

            else:
                lineage = "mixed"
                sublineage = ""
                lineage_cov = None
                lineage_depth = None
                for i in lineage_results:
                    sublineage += i["Name"] + "/"
                sublineage = sublineage[:-1]
            row.append(n_lineages)
            row.append(lineage)
            row.append(sublineage)
            row.append(lineage_cov)
            row.append(lineage_depth)

            antibiogram = ""
            amr_results = data["Genomes"][0]["Resistance Prediction"][
                "Resistance Prediction Summary"
            ]
            for j in amr_results:
                for k in amr_results[j]:
                    antibiogram += amr_results[j][k]
            row.append(antibiogram)

            pipeline_build = data["Metadata"]["Pipeline build"].replace("\n", "")
            row.append(pipeline_build)

            tables.append(row)

        genomes = pandas.DataFrame(
            tables,
            columns=[
                "uniqueid",
                "mycobacterial_reads",
                "name",
                "tb_reads",
                "tb_coverage",
                "tb_depth",
                "n_lineages",
                "lineage",
                "sublineage",
                "mykrobe_lineage_coverage",
                "mykrobe_lineage_depth",
                "antibiogram",
                "pipeline_build",
            ],
        )
        genomes["species"] = genomes.apply(split_species, axis=1)
        genomes.drop(columns=["name"], inplace=True)
        genomes = genomes[
            [
                "uniqueid",
                "species",
                "n_lineages",
                "lineage",
                "sublineage",
                "mycobacterial_reads",
                "tb_reads",
                "tb_coverage",
                "tb_depth",
                "mykrobe_lineage_coverage",
                "mykrobe_lineage_depth",
                "antibiogram",
                "pipeline_build",
            ]
        ]

        genomes.to_csv(tables_path / "genomes.csv")
        genomes.to_parquet(tables_path / "genomes.parquet")
        total_genomes = successful_genome + too_few_reads_genome + too_few_reads_id
        print(f"{total_genomes} samples were processed.")
        print(
            f"{successful_genome} successfully reached a genome, {too_few_reads_id} samples had too few reads for identification and {too_few_reads_genome} samples had too few reads for genome assembly"
        )

    master_table.to_csv(tables_path / "ENA_LOOKUP.csv", index=False)
    master_table.to_parquet(tables_path / "ENA_LOOKUP.parquet")


def main():
    import defopt

    defopt.run(
        [build_tables],
        no_negated_flags=True,
        strict_kwonly=False,
        short={},
    )


if __name__ == "__main__":
    main()
