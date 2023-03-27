from collections import defaultdict
from functools import singledispatchmethod
from urllib.parse import urljoin

import requests


class EnsemblRestClient:
    headers = defaultdict(str)

    def __init__(self, assembly="GRCh38"):
        self.session = requests.Session()
        if assembly == "GRCh38":
            self.server = "https://rest.ensembl.org"
        elif assembly == "GRCh37":
            self.server = "https://grch37.rest.ensembl.org"

    def get(self, endpoint, params, format):
        if format == 'json':
            self.headers['content-type'] = 'application/json'
            response = self.session.get(urljoin(self.server, endpoint), headers=self.headers, params=params)
            return response.json()
        elif format == 'xml':
            self.headers['content-type'] = 'text/xml'
            response = self.session.get(urljoin(self.server, endpoint), headers=self.headers, params=params)
            return response.text

    def post(self, endpoint, params, json, format):
        if format == 'json':
            self.headers['content-type'] = 'application/json'
            response = self.session.post(urljoin(self.server, endpoint), params=params, json=json, headers=self.headers)
            return response.json()
        elif format == "xml":
            self.headers['content-type'] = 'text/xml'
            response = self.session.post(urljoin(self.server, endpoint), params=params, json=json, headers=self.headers)
            return response.text

    @ singledispatchmethod
    def variant_recoder(self, id: str, species="human", format='json', **kwargs):
        return self.get(endpoint=f"variant_recoder/{species}/{id}", format=format, params=kwargs)

    @ variant_recoder.register
    def _(self, id: list, species="human", format='json', **kwargs):
        return self.post(endpoint=f"variant_recoder/{species}", format=format, params=kwargs, json={"ids": id})

    @ singledispatchmethod
    def variation(self, id: str, species="human", format='json', **kwargs):
        return self.get(endpoint=f"variation/{species}/{id}", format=format, params=kwargs)

    @ variation.register
    def _(self, id: list, species="human", format='json', **kwargs):
        return self.post(endpoint=f"variation/{species}", format=format, params=kwargs, json={"ids": id})

    def variation_pmcid(self, pmcid, species="human", format='json'):
        return self.get(endpoint=f"variation/{species}/pmcid/{pmcid}", format=format)

    def variation_pmid(self, pmid, species="human", format='json'):
        return self.get(endpoint=f"variation/{species}/pmid/{pmid}", format=format)

    @ singledispatchmethod
    def vep_hgvs(self, hgvs: str, species="human", format='json', **kwargs):
        return self.get(endpoint=f"vep/{species}/hgvs/{hgvs}", params=kwargs, format=format)

    @ vep_hgvs.register
    def _(self, hgvs: list, species="human", format='json', **kwargs):
        return self.post(endpoint=f"vep/{species}/hgvs", params=kwargs, format=format, json={"hgvs_notations": hgvs})

    @ singledispatchmethod
    def vep_id(self, id: str, species="human", format='json', **kwargs):
        return self.get(endpoint=f"vep/{species}/id/{id}", params=kwargs, format=format)

    @ vep_id.register
    def _(self, id: list, species="human", format='json', **kwargs):
        return self.post(endpoint=f"vep/{species}/id", params=kwargs, json={"ids": id}, format=format)


if __name__ == "__main__":
    import pprint

    ensembl = EnsemblRestClient()
    pprint.pprint(ensembl.vep_id(["rs137853119", "rs137853120"]))
