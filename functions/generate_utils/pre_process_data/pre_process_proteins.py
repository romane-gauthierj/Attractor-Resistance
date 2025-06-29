# Pre-process proteins data (only id in the patients id list  and proteins in montagud list)


def process_proteins(proteins_data, montagud_nodes, patients_id):
    proteins_data["sample"] = proteins_data["sample"].str.rsplit("-", n=2).str[0]
    proteins_data["sample"] = proteins_data["sample"].str.upper()
    proteins_data["sample"] = proteins_data["sample"].str.replace("_", "", regex=False)
    proteins_data["sample"] = proteins_data["sample"].str.replace("-", "", regex=False)

    pattern = "|".join(montagud_nodes)

    # Filter rows where 'sample' contains any name from the list
    proteins_data = proteins_data[
        proteins_data["sample"].str.contains(pattern, case=False, na=False)
    ]
    proteins_data = proteins_data.dropna(how="all", subset=proteins_data.columns[1:])

    mods = ["_PS", "_PT", "_PY"]
    proteins_data = proteins_data[
        ~proteins_data["sample"].apply(lambda p: any(mod in p for mod in mods))
    ]

    # Keep only the proteins present in the montagud list

    proteins_data_col = list(proteins_data.columns)
    common_col = list(set(proteins_data_col) & set(patients_id))
    col_keep = ["sample"] + common_col
    proteins_data_filtered = proteins_data[col_keep]

    df_melted_protein = proteins_data_filtered.melt(
        id_vars=["sample"],  # columns to keep fixed
        var_name="samples_id",  # name for the variable column (sample IDs)
        value_name="expression_value",  # name for the values
    )

    df_melted_protein["sample"] = df_melted_protein["sample"].str.split("|").str[0]

    df_melted_protein = df_melted_protein.rename(
        columns={
            "samples_id": "model_id",
            "sample": "protein_symbol",
            "expression_value": "rsem_tpm",
        }
    )
    df_melted_protein["protein_symbol"] = df_melted_protein[
        "protein_symbol"
    ].str.upper()

    def replace_with_base_name(protein_name):
        for base in montagud_nodes:
            if protein_name.startswith(base):
                return base
        return protein_name  # if no match found, keep original

    # Assuming your dataframe is df and column to replace is 'protein_symbol'
    df_melted_protein["protein_symbol"] = df_melted_protein["protein_symbol"].apply(
        replace_with_base_name
    )

    df_melted_protein = df_melted_protein[
        df_melted_protein["protein_symbol"].isin(montagud_nodes)
    ]
    df_melted_protein["protein_symbol"] = df_melted_protein[
        "protein_symbol"
    ].str.replace("_", "", regex=False)
    df_melted_protein = df_melted_protein[df_melted_protein["rsem_tpm"].notna()]
    df_melted_protein.to_csv(
        "data/TCGA_data/prostate/filtered_data/proteins_samples_table.csv"
    )
    return df_melted_protein


def classify_expression(z):
    if z > 1:
        return "high"
    elif z < -1:
        return "low"
    else:
        return "normal"


def create_table_proteins_patients(proteins_data):
    # filter only data with patient id and the nodes of the montagud_data

    proteins_data["z_score"] = proteins_data.groupby("protein_symbol")[
        "rsem_tpm"
    ].transform(lambda x: (x - x.mean()) / x.std())

    proteins_data["protein_expression_level"] = proteins_data["z_score"].apply(
        classify_expression
    )

    proteins_data = proteins_data[
        ["model_id", "protein_symbol", "protein_expression_level"]
    ]
    proteins_data.rename(columns={"protein_symbol": "protein_name"}, inplace=True)
    proteins_data = proteins_data[
        proteins_data["protein_expression_level"].isin(["low", "high"])
    ]

    table_proteins_patients = proteins_data.pivot_table(
        index="model_id",
        columns="protein_expression_level",
        values="protein_name",
        aggfunc=lambda x: ", ".join(sorted(set(x))),
    ).fillna("-")

    table_proteins_patients = table_proteins_patients.rename(
        columns={"low": "Low Protein Abundance", "high": "High Protein Abundance"}
    )

    return table_proteins_patients
