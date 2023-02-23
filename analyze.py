from collections import Counter
from functools import lru_cache
from pathlib import Path

import bioregistry
import pandas as pd
import pyobo
import pystow
import requests
from libsbml import SBMLReader
from lxml import etree
from tabulate import tabulate
from tqdm import tqdm
from tqdm.contrib.concurrent import process_map, thread_map

HERE = Path(__file__).parent.resolve()
BIOMODELS = pystow.module("mira", "biomodels")

DOWNLOAD_URL = "https://www.ebi.ac.uk/biomodels/search/download"

#: Skip looking up names from these resources
SKIP_NAMES = {
    # Looking up KEGG names is a problem
    "kegg",
    # Looking up Biomodels names is redundant
    "biomodels.db",
    # The following are publications - we don't necessarily need their titles
    "pubmed",
    "doi",
    "arxiv",
}

#: Use this prefix map when parsing XML
PREFIX_MAP = {
    "rdf": "http://www.w3.org/1999/02/22-rdf-syntax-ns#",
    "dcterms": "http://purl.org/dc/terms/",
    "vCard": "http://www.w3.org/2001/vcard-rdf/3.0#",
    "vCard4": "http://www.w3.org/2006/vcard/ns#",
    "bqbiol": "http://biomodels.net/biology-qualifiers/",
    "bqmodel": "http://biomodels.net/model-qualifiers/",
    "CopasiMT": "http://www.copasi.org/RDF/MiriamTerms#",
    "copasi": "http://www.copasi.org/static/sbml",
    "jd": "http://www.sys-bio.org/sbml",
}

#: This key is used to extract the RDF content from a given annotation
RESOURCE_KEY = "{http://www.w3.org/1999/02/22-rdf-syntax-ns#}resource"


@lru_cache(maxsize=1)
def get_converter():
    return bioregistry.get_converter(include_prefixes=True)


def parse_uri(uri):
    return get_converter().parse_uri(uri)


def _get_mamo_luid_to_name() -> dict[str, str]:
    mamo_properties_df = pyobo.get_properties_df("mamo")
    label_idx = mamo_properties_df["property"] == "prefLabel"
    mamo_properties_df = mamo_properties_df[label_idx]
    return dict(mamo_properties_df[["mamo_id", "value"]].values)


MAMO_NAMES = _get_mamo_luid_to_name()
del _get_mamo_luid_to_name


def get_name_by_curie(curie: str) -> str:
    if any(curie.startswith(x) for x in SKIP_NAMES):
        return ""
    if curie.startswith("mamo:"):
        return MAMO_NAMES.get(curie.removeprefix("mamo:")) or ""
    try:
        rv = pyobo.get_name_by_curie(curie) or ""
    except Exception:
        return ""
    else:
        return rv


def get_sbml_from_file(file):
    """Return an SBML model object from an SBML file."""
    sbml_string = file.read().decode("utf-8")
    sbml_document = SBMLReader().readSBMLFromString(sbml_string)
    return sbml_document.getModel()


def _download_model_metadata(model_id: str):
    model_module = BIOMODELS.module("models", model_id)
    url = f"https://www.ebi.ac.uk/biomodels/{model_id}?format=json"
    data = model_module.ensure_json(
        url=url,
        name=f"{model_id}_metadata.json",
        download_kwargs={"progress_bar": False},
    )
    return model_id, data


def _download_model(model_id: str):
    url = f"{DOWNLOAD_URL}?models={model_id}"
    return BIOMODELS.ensure(
        "models",
        model_id,
        url=url,
        name=f"{model_id}.zip",
        download_kwargs={"progress_bar": False},
    )


def main(log_failure: bool = False):
    # Get list of all models
    identifiers_res = requests.get(
        "https://www.ebi.ac.uk/biomodels/model/identifiers",
        params={"format": "json"},
    )
    model_ids = sorted(identifiers_res.json()["models"])

    # Download all model metadata
    model_info = thread_map(
        _download_model_metadata,
        model_ids,
        desc="Downloading model metadata",
        unit_scale=True,
        unit="model",
    )

    skip = {
        "MODEL1311110001",  # annoyingly big
        "MODEL1703310000",  # annoyingly big
    }

    # Only keep models in SBML format from now on
    model_ids = [
        model_id
        for model_id, data in model_info
        if model_id in skip or data.get("format", {}).get("name") == "SBML"
    ]
    thread_map(
        _download_model,
        model_ids,
        desc="Downloading model zip",
        unit_scale=True,
        unit="model",
    )
    parsed_models = process_map(
        _cache_annotations_str,
        model_ids,
        desc="Parsing models",
        unit_scale=True,
        unit="model",
        chunksize=300,
    )

    triples = []
    example_tag_to_model = {}
    example_tag_to_curie = {}
    example_tag_prefix_to_identifier = {}
    counter = Counter()
    tag_prefix_counter = Counter()

    for model_id, annotations_string in tqdm(
        parsed_models, desc="Analyzing models", unit_scale=True, unit="model"
    ):
        if not annotations_string:
            continue

        try:
            et = etree.fromstring(annotations_string)
        except etree.XMLSyntaxError as e:
            if log_failure:
                tqdm.write(f"[{model_id}] failed to parse: {e}")
            continue

        description = et.find("rdf:RDF/rdf:Description", namespaces=PREFIX_MAP)
        if description is None:
            continue

        for element in description:
            counter[element.tag] += 1
            if element.tag not in example_tag_to_model:
                example_tag_to_model[element.tag] = model_id
                # tqdm.write(f"{model_id}: {element.tag}")
            for subelement in element.findall("rdf:Bag/rdf:li", namespaces=PREFIX_MAP):
                uri = subelement.attrib.get(RESOURCE_KEY)
                if uri:
                    prefix, identifier = parse_uri(uri)
                    if prefix == "kegg":
                        prefix = "kegg.pathway"
                    if not prefix:
                        tqdm.write(f"{model_id} unparsable URI: {uri}")
                    elif prefix == "idot":
                        tqdm.write(f"parsed {uri} wrong: {prefix}:{identifier}")
                    else:
                        tag_prefix_counter[element.tag, prefix] += 1
                        curie = bioregistry.curie_to_str(prefix, identifier)
                        name = get_name_by_curie(curie)
                        triples.append(
                            (
                                model_id,
                                str(element.tag).replace("{", "").replace("}", ""),
                                curie,
                                name,
                            )
                        )
                        if (
                            element.tag,
                            prefix,
                        ) not in example_tag_prefix_to_identifier:
                            example_tag_prefix_to_identifier[element.tag, prefix] = (
                                model_id,
                                identifier,
                            )
                            tqdm.write(
                                f"{model_id}: {element.tag} example: {curie} ({name})"
                            )
                        if element.tag not in example_tag_to_curie:
                            example_tag_to_curie[element.tag] = curie

    rows1 = [
        (tag, count, example_tag_to_model[tag], example_tag_to_curie.get(tag, ""))
        for tag, count in counter.most_common()
    ]
    header1 = ["tag", "count", "example_model", "example_value"]
    df1 = pd.DataFrame(rows1, columns=header1)
    df1.to_csv(HERE.joinpath("tag_summary.tsv"), index=False, sep="\t")
    print(
        tabulate(
            rows1,
            headers=header1,
            tablefmt="github",
        )
    )

    rows2 = [
        (
            tag,
            prefix,
            count,
            *example_tag_prefix_to_identifier.get((tag, prefix), ("", "")),
        )
        for (tag, prefix), count in tag_prefix_counter.most_common()
    ]
    headers2 = ["tag", "prefix", "count", "example_model", "example_luid"]
    df2 = pd.DataFrame(rows2, columns=headers2)
    df2.to_csv(HERE.joinpath("tag_prefix_summary.tsv"), index=False, sep="\t")
    print(
        tabulate(
            rows2,
            headers=headers2,
            tablefmt="github",
        )
    )

    triples_df = pd.DataFrame(
        triples, columns=["biomodels_luid", "predicate_uri", "object_curie", "object_label"]
    )
    triples_df.to_csv(HERE.joinpath("triples.tsv"), index=False, sep="\t")


def _cache_annotations_str(model_id: str, log_failure: bool = False) -> tuple[str, str]:
    path: Path = BIOMODELS.join("models", model_id, name=f"{model_id}_annotations.xml")
    if path.is_file():
        return model_id, path.read_text().strip()

    url = f"{DOWNLOAD_URL}?models={model_id}"
    try:
        with BIOMODELS.ensure_open_zip(
            "models",
            model_id,
            url=url,
            name=f"{model_id}.zip",
            inner_path=f"{model_id}.xml",
        ) as file:
            try:
                sbml_model = get_sbml_from_file(file)
            except Exception as e:
                tqdm.write(f"[{model_id}] failed to parse: {e}")
                sbml_model = None
    except KeyError:
        sbml_model = None

    if sbml_model is None:
        if log_failure:
            tqdm.write(f"[{model_id}] did not parse")
        path.write_text("")
        return model_id, ""

    annotations_string = sbml_model.getAnnotationString() or ""
    path.write_text(annotations_string)
    return model_id, annotations_string


if __name__ == "__main__":
    main()
