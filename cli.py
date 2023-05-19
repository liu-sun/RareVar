import argparse
import pprint
from collections import defaultdict
from functools import singledispatch
from urllib.parse import urljoin

import requests
from requests.adapters import HTTPAdapter, Retry

assembly = "GRCh38"
scheme = "https"
match assembly, scheme:
    case "GRCh38", "http":
        server = "http://rest.ensembl.org"
    case "GRCh37", "http":
        server = "http://grch37.rest.ensembl.org"
    case "GRCh38", "https":
        server = "https://rest.ensembl.org"
    case "GRCh37", "https":
        server = "https://grch37.rest.ensembl.org"


media_type = dict(
    json="application/json",
    xml="text/xml",
    nh="text/x-nh",
    phyloxml="text/x-phyloxml+xml",
    orthoxml="text/x-orthoxml+xml",
    gff3="text/x-gff3",
    fasta="text/x-fasta",
    bed="text/x-bed",
    seqxml="text/x-seqxml+xml",
    text="text/plain",
    yaml="text/x-yaml",
    jsonp="text/javascript")


session = requests.Session()
adapter = HTTPAdapter(max_retries=Retry(backoff_factor=3600/55000,
                                        respect_retry_after_header=True, status_forcelist=[429], allowed_methods=["GET", "POST"]))
session.mount(server, adapter)


def get(endpoint, params, format):
    headers = defaultdict(str)
    headers["Content-Type"] = media_type[format]
    response = session.get(urljoin(server, endpoint),
                           headers=headers, params=params)
    if response.ok:
        if headers["Content-Type"] == "application/json":
            return response.json()
        else:
            return response.text
    else:
        response.raise_for_status()


def post(endpoint, params, json, format):
    headers = defaultdict(str)
    headers["Content-Type"] = media_type[format]
    headers['Accept'] = media_type[format]
    response = session.post(urljoin(server, endpoint),
                            headers=headers, params=params, json=json)
    if response.ok:
        if headers["Accept"] == "application/json":
            return response.json()
        else:
            return response.text
    else:
        response.raise_for_status()


@singledispatch
def variant_recoder(id: str, species='human', format="json", **kwargs):
    """Translate a variant identifier, HGVS notation or genomic SPDI notation to all possible variant IDs, HGVS and genomic SPDI"""
    return get(endpoint=f"variant_recoder/{species}/{id}", params=kwargs, format=format)


@variant_recoder.register
def _(id: list, species='human', format="json", **kwargs):
    """Translate a list of variant identifiers, HGVS notations or genomic SPDI notations to all possible variant IDs, HGVS and genomic SPDI"""
    return post(endpoint=f"variant_recoder/{species}", params=kwargs, json={"ids": id}, format=format)


@singledispatch
def variation(id: str, species='human', format="json", **kwargs):
    """Uses a variant identifier (e.g. rsID) to return the variation features including optional genotype, phenotype and population data"""
    return get(endpoint=f"variation/{species}/{id}", params=kwargs, format=format)


@variation.register
def _(id: list, species='human', format="json", **kwargs):
    """Uses a list of variant identifiers (e.g. rsID) to return the variation features including optional genotype, phenotype and population data"""
    return post(endpoint=f"variation/{species}", params=kwargs, json={"ids": id}, format=format)


def variation_pmcid(pmcid, species='human', format="json", **kwargs):
    """Fetch variants by publication using PubMed Central reference number (PMCID)"""
    return get(endpoint=f"variation/{species}/pmcid/{pmcid}", params=kwargs, format=format)


def variation_pmid(pmid, species='human', format="json", **kwargs):
    """Fetch variants by publication using PubMed reference number (PMID)"""
    return get(endpoint=f"variation/{species}/pmid/{pmid}", params=kwargs, format=format)


@singledispatch
def vep_hgvs(hgvs_notation: str, species='human', format="json", **kwargs):
    """Fetch variant consequences based on a HGVS notation"""
    return get(endpoint=f"vep/{species}/hgvs/{hgvs_notation}", params=kwargs, format=format)


@vep_hgvs.register
def _(hgvs: list, species='human', format="json", **kwargs):
    """Fetch variant consequences for multiple HGVS notations"""
    return post(endpoint=f"vep/{species}/hgvs", params=kwargs, json={"hgvs_notations": hgvs}, format=format)


@singledispatch
def vep_id(id: str, species='human', format="json", **kwargs):
    """Fetch variant consequences based on a variant identifier"""
    return get(endpoint=f"vep/{species}/id/{id}", params=kwargs, format=format)


@vep_id.register
def _(id: list, species='human', format="json", **kwargs):
    """Fetch variant consequences for multiple ids"""
    return post(endpoint=f"vep/{species}/id", params=kwargs, json={"ids": id}, format=format)


parser = argparse.ArgumentParser(description='Ensembl REST API client', fromfile_prefix_chars='@',
                                 prog="ensembl_rest", epilog="For further help on a command, type ensembl_rest <command> -h")
parser.add_argument('--assembly', default="GRCh38", help='Genome assembly')
parser.add_argument('--scheme', default="https", help='Scheme')
subparsers = parser.add_subparsers(dest="command")
parser_variant_recoder = subparsers.add_parser(
    'variant_recoder', help='variant_recoder help')
parser_variant_recoder.add_argument(
    'id', help='Variant identifier, HGVS notation or genomic SPDI notation', action="extend", nargs="+")
parser_variant_recoder.add_argument(
    '--species', default="human", help='Species')
parser_variant_recoder.add_argument('--format', default="json", help='Format')
parser_variation = subparsers.add_parser('variation', help='variation help')
parser_variation.add_argument(
    'id', help='Variant identifier', action="extend", nargs="+")
parser_variation.add_argument('--species', default="human", help='Species')
parser_variation.add_argument('--format', default="json", help='Format')
parser_variation_pmcid = subparsers.add_parser("variation_pmcid")
parser_variation_pmcid.add_argument('pmcid', help='PMCID')
parser_variation_pmcid.add_argument(
    '--species', default="human", help='Species')
parser_variation_pmcid.add_argument('--format', default="json", help='Format')
parser_variation_pmid = subparsers.add_parser("variation_pmid")
parser_variation_pmid.add_argument('pmid', help='PMID')
parser_variation_pmid.add_argument(
    '--species', default="human", help='Species')
parser_variation_pmid.add_argument('--format', default="json", help='Format')
parser_vep_hgvs = subparsers.add_parser("vep_hgvs")
parser_vep_hgvs.add_argument(
    'hgvs_notation', help='HGVS notation', action="extend", nargs="+")
parser_vep_hgvs.add_argument('--species', default="human", help='Species')
parser_vep_hgvs.add_argument('--format', default="json", help='Format')
parser_vep_id = subparsers.add_parser("vep_id")
parser_vep_id.add_argument(
    'id', help='Variant identifier', action="extend", nargs="+")
parser_vep_id.add_argument('--species', default="human", help='Species')
parser_vep_id.add_argument('--format', default="json", help='Format')
parser_variation.set_defaults(func=variation)
parser_variation_pmcid.set_defaults(func=variation_pmcid)
parser_variation_pmid.set_defaults(func=variation_pmid)
parser_vep_hgvs.set_defaults(func=vep_hgvs)
parser_vep_id.set_defaults(func=vep_id)
parser_variant_recoder.set_defaults(func=variant_recoder)
args = parser.parse_args()
assembly = args.assembly
scheme = args.scheme
match assembly, scheme:
    case "GRCh38", "http":
        server = "http://rest.ensembl.org"
    case "GRCh37", "http":
        server = "http://grch37.rest.ensembl.org"
    case "GRCh38", "https":
        server = "https://rest.ensembl.org"
    case "GRCh37", "https":
        server = "https://grch37.rest.ensembl.org"
match args.command:
    case "variant_recoder":
        pprint.pprint(
            args.func(args.id, species=args.species, format=args.format))
    case "variation":
        pprint.pprint(
            args.func(args.id, species=args.species, format=args.format))
    case "variation_pmcid":
        pprint.pprint(
            args.func(args.pmid, species=args.species, format=args.format))
    case "variation_pmid":
        pprint.pprint(
            args.func(args.pmcid, species=args.species, format=args.format))
    case "vep_hgvs":
        pprint.pprint(args.func(args.hgvs_notation,
                      species=args.species, format=args.format))
    case "vep_id":
        pprint.pprint(
            args.func(args.id, species=args.species, format=args.format))
