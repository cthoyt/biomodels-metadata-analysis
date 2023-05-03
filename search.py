import pandas as pd
from collections import defaultdict

if __name__ == '__main__':
    triples_df = pd.read_csv("triples.tsv", sep='\t')
    triples_df = triples_df[triples_df["object_curie"].str.startswith("mamo:")]
    sdf = triples_df.groupby(["object_curie", "object_label"]).count()["predicate_uri"].reset_index()
    sdf.sort_values("predicate_uri", inplace=True, ascending=False)
    sdf.rename(columns={"predicate_uri": "count"}, inplace=True)
    sdf.to_csv("mamo_summary.tsv", sep='\t', index=False)
    # print(sdf)

    metadata_df = pd.read_csv("metadata.tsv", sep='\t')
    model_names = dict(metadata_df[["biomodels_luid", "name"]].values)
    model_years = dict(metadata_df[["biomodels_luid", "publication_year"]].values)

    # filter to population models
    target_columns = ["biomodels_luid", "name", "publication_year"]
    population_model_df: pd.DataFrame = triples_df[triples_df["object_curie"] == "mamo:0000028"].copy()
    population_model_df["name"] = population_model_df["biomodels_luid"].map(model_names)
    population_model_df["publication_year"] = population_model_df["biomodels_luid"].map(model_years)
    population_model_df = population_model_df[target_columns]
    print(population_model_df)
    population_model_df.to_csv("population_models.tsv", sep='\t', index=False)

    names = {
        "sir", "sird", "seir", "sei", "epide", "sihrd", "seiahrd",}
    name_idx = metadata_df["name"].map(lambda s: any(x in s.lower() for x in names))
    sdf = metadata_df.loc[name_idx][target_columns]
    sdf = sdf[~sdf.index.isin(population_model_df.index)]
    sdf.to_csv("name_search.tsv", sep='\t', index=False)

    hits = [
        "BIOMD0000000715",  # SEIS epidemic model with the impact of media
        "BIOMD0000001045",  # hong kong flu
        "MODEL1805220001",  # Human/Mosquito SEIR/SEI Mode
        "MODEL1805230001",  # Model for HIV-Malaria co-infection
        "MODEL1808280006",  # SIRWS model with immune boosting and cross-immunity between two pathogens
        "MODEL1008060002",  # zombie infection toy model (lol)
    ]
