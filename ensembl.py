from collections import defaultdict
from functools import singledispatch, singledispatchmethod
from urllib.parse import urljoin

import requests
from requests.adapters import HTTPAdapter, Retry

assembly = "GRCh38"
scheme = "https"


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


match assembly, scheme:
    case "GRCh38", "http":
        server = "http://rest.ensembl.org"
    case "GRCh37", "http":
        server = "http://grch37.rest.ensembl.org"
    case "GRCh38", "https":
        server = "https://rest.ensembl.org"
    case "GRCh37", "https":
        server = "https://grch37.rest.ensembl.org"

session = requests.Session()
adapter = HTTPAdapter(max_retries=Retry(backoff_factor=3600/55000,
                                        respect_retry_after_header=True, status_forcelist=[429], allowed_methods=["GET", "POST"]))
session.mount(server, adapter)


def get(endpoint, params, format):
    headers = defaultdict(str)
    headers["Content-Type"] = media_type[format]
    response = session.get(urljoin(server, endpoint), headers=headers, params=params)
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
    response = session.post(urljoin(server, endpoint), headers=headers, params=params, json=json)
    if response.ok:
        if headers["Accept"] == "application/json":
            return response.json()
        else:
            return response.text
    else:
        response.raise_for_status()


class Ensembl:
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

    def __init__(self, assembly="GRCh38", scheme="https"):
        self.session = requests.Session()
        adapter = HTTPAdapter(max_retries=Retry(backoff_factor=3600/55000,
                                                respect_retry_after_header=True, status_forcelist=[429], allowed_methods=["GET", "POST"]))
        match assembly, scheme:
            case "GRCh38", "http":
                self.server = "http://rest.ensembl.org"
            case "GRCh37", "http":
                self.server = "http://grch37.rest.ensembl.org"
            case "GRCh38", "https":
                self.server = "https://rest.ensembl.org"
            case "GRCh37", "https":
                self.server = "https://grch37.rest.ensembl.org"
        self.session.mount(self.server, adapter)

    def get(self, endpoint, params, format):
        headers = defaultdict(str)
        headers["Content-Type"] = self.media_type[format]
        response = self.session.get(urljoin(self.server, endpoint), headers=headers, params=params)
        if response.ok:
            if headers["Content-Type"] == "application/json":
                return response.json()
            else:
                return response.text
        else:
            response.raise_for_status()

    def post(self, endpoint, params, json, format):
        headers = defaultdict(str)
        headers["Content-Type"] = self.media_type[format]
        headers['Accept'] = self.media_type[format]
        response = self.session.post(urljoin(self.server, endpoint), headers=headers, params=params, json=json)
        if response.ok:
            if headers["Accept"] == "application/json":
                return response.json()
            else:
                return response.text
        else:
            response.raise_for_status()

    @singledispatchmethod
    def variant_recoder(self, id: str, species="human", format='json', **kwargs):
        return self.get(endpoint=f"variant_recoder/{species}/{id}", format=format, params=kwargs)

    @variant_recoder.register
    def _(self, id: list, species="human", format='json', **kwargs):
        return self.post(endpoint=f"variant_recoder/{species}", format=format, params=kwargs, json={"ids": id})

    @singledispatchmethod
    def variation(self, id: str, species="human", format='json', **kwargs):
        return self.get(endpoint=f"variation/{species}/{id}", format=format, params=kwargs)

    @variation.register
    def _(self, id: list, species="human", format='json', **kwargs):
        return self.post(endpoint=f"variation/{species}", format=format, params=kwargs, json={"ids": id})

    def variation_pmcid(self, pmcid, species="human", format='json'):
        return self.get(endpoint=f"variation/{species}/pmcid/{pmcid}", format=format)

    def variation_pmid(self, pmid, species="human", format='json'):
        return self.get(endpoint=f"variation/{species}/pmid/{pmid}", format=format)

    @singledispatchmethod
    def vep_hgvs(self, hgvs: str, species="human", format='json', **kwargs):
        return self.get(endpoint=f"vep/{species}/hgvs/{hgvs}", params=kwargs, format=format)

    @vep_hgvs.register
    def _(self, hgvs: list, species="human", format='json', **kwargs):
        return self.post(endpoint=f"vep/{species}/hgvs", params=kwargs, format=format, json={"hgvs_notations": hgvs})

    @singledispatchmethod
    def vep_id(self, id: str, species="human", format='json', **kwargs):
        return self.get(endpoint=f"vep/{species}/id/{id}", params=kwargs, format=format)

    @vep_id.register
    def _(self, id: list, species="human", format='json', **kwargs):
        return self.post(endpoint=f"vep/{species}/id", params=kwargs, json={"ids": id}, format=format)


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
def vep_hgvs(hgvs: str, species='human', format="json", **kwargs):
    """Fetch variant consequences based on a HGVS notation"""
    return get(endpoint=f"vep/{species}/hgvs/{hgvs}", params=kwargs, format=format)


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


@singledispatch
def vep_region(region: str, species='human', format="json", **kwargs):
    """Fetch variant consequences based on a genomic region"""
    return get(endpoint=f"vep/{species}/region/{region}", params=kwargs, format=format)


@vep_region.register
def _(region: list, species='human', format="json", **kwargs):
    """Fetch variant consequences for multiple genomic regions"""
    return post(endpoint=f"vep/{species}/region", params=kwargs, json={"variants": region}, format=format)


if __name__ == "__main__":
    import pprint

    ensembl = Ensembl()
    pprint.pprint(ensembl.variant_recoder(["rs137853119", "rs137853120"], fields="id", vcf_string=True))
    pprint.pprint(vep_hgvs(["NP_001361433.1:p.Asp512Asn", "NP_001361433.1:p.Gly433Arg"]))
    pprint.pprint(vep_region("22:37075180-37075180:1/T"))
    pprint.pprint(vep_region("22:37073553-37073553:1/T"))
    pprint.pprint(vep_region(["22 37075180 rs137853119 C T . . .", "22 37073553 rs137853120 C T . . ."]))
