import pandas, pathlib, json

from tqdm import tqdm

from collections import defaultdict


def split_species(row):
    cols = row["name"].split(" (")
    species = cols[0]
    return species


def parse_variants(row):

    variant = row.variant
    is_null = False
    is_minor = False
    minor_variant = None
    minor_reads = 0
    coverage = 0
    # frs = 0

    a = json.loads(row["vcf_evidence"])
    if "COV_TOTAL" in a.keys():
        coverage = a["COV_TOTAL"]
    # if "FRS" in a.keys():
    #     frs = a["FRS"]
    #     if frs == ".":
    #         frs = None

    if variant[-1] == "x":
        is_null = True
    elif ":" in variant:
        is_minor = True
        cols = variant.split(":")
        minor_variant = cols[0]
        minor_reads = int(cols[1])
        if "ins" in variant or "del" in variant:
            variant = minor_variant.split("_")[0] + "_minorindel"
        else:
            variant = minor_variant[:-1] + "z"

    return pandas.Series(
        [variant, is_null, is_minor, minor_variant, minor_reads, coverage]
    )


def parse_mutations(row):

    mutation = row.mutation
    is_null = False
    is_minor = False
    minor_mutation = None
    minor_reads = None

    if ":" in mutation:
        is_minor = True
        cols = mutation.split(":")
        minor_mutation = cols[0]
        minor_reads = int(cols[1])
        if "ins" in mutation or "del" in mutation:
            mutation = mutation.split("_")[0] + "_minorindel"
        else:
            if minor_mutation[0] in ["a", "t", "c", "g"]:
                mutation = minor_mutation[:-1] + "z"
            else:
                mutation = minor_mutation[:-1] + "Z"
    elif mutation[-1] in ["X", "x"]:
        is_null = True

    return pandas.Series([mutation, is_null, is_minor, minor_mutation, minor_reads])


def build_tables(
    lookup_table: str = None,
    source_files: str = "data/",
    max_samples: int = None,
    output: str = None,
    filename: str = None,
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
        "effects",
        "variants",
        "mutations",
        "predictions",
        "main_report",
    ], "can only specify one from this list"

    if filename in ["effects", "mutations", "predictions", "variants"]:

        tables = []
        n_samples = 0

        for i in tqdm((path).rglob("*" + filename + ".csv")):

            df = pandas.read_csv(i)
            n_samples += 1
            if max_samples is not None and n_samples > max_samples:
                break
            uid = i.stem.split("." + filename)[0]
            df["uniqueid"] = uid
            master_table.at[uid, "has_" + filename] = True
            tables.append(df)

        df = pandas.concat(tables)

        if filename == "effects":
            for col in [
                "gene",
                "drug",
                "prediction",
                "catalogue_name",
                "prediction_values",
            ]:
                df[col] = df[col].astype("category")
            df.set_index(
                [
                    "uniqueid",
                    "catalogue_name",
                    "catalogue_version",
                    "prediction_values",
                    "drug",
                    "gene",
                    "mutation",
                ],
                inplace=True,
            )

        elif filename == "predictions":
            for col in [
                "drug",
                "prediction",
                "catalogue_name",
                "catalogue_values",
            ]:
                df[col] = df[col].astype("category")
            df.set_index(
                [
                    "uniqueid",
                    "catalogue_name",
                    "catalogue_version",
                    "catalogue_values",
                    "drug",
                ],
                inplace=True,
            )
        elif filename == "variants":
            df[
                [
                    "var",
                    "is_null",
                    "is_minor",
                    "minor_variant",
                    "minor_reads",
                    "coverage",
                ]
            ] = df.apply(parse_variants, axis=1)
            df.drop(columns=["variant"], inplace=True)
            df.rename(columns={"var": "variant"}, inplace=True)
            for col in [
                "gene",
            ]:
                df[col] = df[col].astype("category")
            df.set_index(["uniqueid", "gene", "variant"], inplace=True)

        elif filename == "mutations":
            df[["mut", "is_null", "is_minor", "minor_mutation", "minor_reads"]] = (
                df.apply(parse_mutations, axis=1)
            )
            df.drop(columns=["mutation"], inplace=True)
            df.rename(columns={"mut": "mutation"}, inplace=True)
            for col in ["gene", "ref", "alt", "amino_acid_sequence", "minor_mutation"]:
                df[col] = df[col].astype("category")
            df.set_index(["uniqueid", "gene", "mutation"], inplace=True)

        df.to_csv(str(tables_path / filename) + ".csv", index=False)
        # if filename != "variants":
        df.to_parquet(str(tables_path / filename) + ".parquet")

    successful_genome = 0
    too_few_reads_genome = 0
    too_few_reads_id = 0
    if filename == "main_report":

        tables = []
        n_samples = 0

        for i in tqdm((path).rglob("*" + filename + ".json")):

            uid = i.stem.split("." + filename)[0]
            master_table.at[uid, "has_" + filename] = True

            row = [uid]
            f = open(i)
            data = json.load(f)
            n_samples += 1
            if max_samples is not None and n_samples > max_samples:
                break

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
        genomes["lineage"] = genomes["lineage"].astype("category")
        genomes["sublineage"] = genomes["sublineage"].astype("category")
        genomes["antibiogram"] = genomes["antibiogram"].astype("category")
        genomes["pipeline_build"] = genomes["pipeline_build"].astype("category")
        genomes.set_index("uniqueid", inplace=True)

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
