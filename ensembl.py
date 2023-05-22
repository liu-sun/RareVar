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
    headers = {}
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
    headers = {}
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
        headers = {}
        headers["Content-Type"] = self.media_type[format]
        response = self.session.get(
            urljoin(self.server, endpoint), headers=headers, params=params)
        if response.ok:
            if headers["Content-Type"] == "application/json":
                return response.json()
            else:
                return response.text
        else:
            response.raise_for_status()

    def post(self, endpoint, params, json, format):
        headers = {}
        headers["Content-Type"] = self.media_type[format]
        headers['Accept'] = self.media_type[format]
        response = self.session.post(
            urljoin(self.server, endpoint), headers=headers, params=params, json=json)
        if response.ok:
            if headers["Accept"] == "application/json":
                return response.json()
            else:
                return response.text
        else:
            response.raise_for_status()

    @singledispatchmethod
    def variant_recoder(self, id: str, species='human', format="json", fields=None, var_synonyms=None, vcf_string=None):
        """Translate a variant identifier, HGVS notation or genomic SPDI notation to all possible variant IDs, HGVS and genomic SPDI"""
        return get(endpoint=f"variant_recoder/{species}/{id}", params=dict(fields=fields, var_synonyms=var_synonyms, vcf_string=vcf_string), format=format)

    @variant_recoder.register
    def _(self, id: list, species='human', format="json", fields=None, var_synonyms=None, vcf_string=None):
        """Translate a list of variant identifiers, HGVS notations or genomic SPDI notations to all possible variant IDs, HGVS and genomic SPDI"""
        return post(endpoint=f"variant_recoder/{species}", params=dict(fields=fields, var_synonyms=var_synonyms, vcf_string=vcf_string), json={"ids": id}, format=format)

    @singledispatchmethod
    def variation(self, id: str, species='human', format="json", pops=None, genotypes=None, genotyping_chips=None, phenotypes=None, population_genotypes=None):
        """Uses a variant identifier (e.g. rsID) to return the variation features including optional genotype, phenotype and population data"""
        return get(endpoint=f"variation/{species}/{id}", params=dict(pops=pops, genotypes=genotypes, genotyping_chips=genotyping_chips, phenotypes=phenotypes, population_genotypes=population_genotypes), format=format)

    @variation.register
    def _(self, id: list, species='human', format="json", pops=None, genotypes=None, phenotypes=None, population_genotypes=None):
        """Uses a list of variant identifiers (e.g. rsID) to return the variation features including optional genotype, phenotype and population data"""
        return post(endpoint=f"variation/{species}", params=dict(pops=pops, genotypes=genotypes, phenotypes=phenotypes, population_genotypes=population_genotypes), json={"ids": id}, format=format)

    def variation_pmcid(self, pmcid, species='human', format="json"):
        """Fetch variants by publication using PubMed Central reference number (PMCID)"""
        return get(endpoint=f"variation/{species}/pmcid/{pmcid}", params=None, format=format)

    def variation_pmid(self, pmid, species='human', format="json"):
        """Fetch variants by publication using PubMed reference number (PMID)"""
        return get(endpoint=f"variation/{species}/pmid/{pmid}", params=None, format=format)

    @singledispatch
    def vep_hgvs(self, hgvs_notation: str, species='human', AncestralAllele=None, Blosum62=None, CADD=None, Conservation=None, DisGeNET=None, EVE=None, GO=None, GeneSplicer=None, IntAct=None, LoF=None, Mastermind=None, MaxEntScan=None,
                 NMD=None, Phenotypes=None, SpliceAI=None, UTRAnnotator=None, ambiguous_hgvs=None, appris=None, callback=None, canonical=None, ccds=None, dbNSFP=None, dbscSNV=None, distance=None, domains=None, failed=None, hgvs=None, mane=None, merged=None, minimal=None,
                 mirna=None, mutfunc=None, numbers=None, protein=None, refseq=None, shift_3prime=None, shift_genomic=None, transcript_id=None, transcript_version=None, tsl=None, uniprot=None, variant_class=None, vcf_string=None, xref_refseq=None, format='json'):
        """Fetch variant consequences based on a HGVS notation"""
        return get(endpoint=f"vep/{species}/hgvs/{hgvs_notation}", params=dict(AncestralAllele=AncestralAllele, Blosum62=Blosum62, CADD=CADD, Conservation=Conservation, DisGeNET=DisGeNET, EVE=EVE, GO=GO, GeneSplicer=GeneSplicer, IntAct=IntAct, LoF=LoF, Mastermind=Mastermind, MaxEntScan=MaxEntScan,
                                                                               NMD=NMD, Phenotypes=Phenotypes, SpliceAI=SpliceAI, UTRAnnotator=UTRAnnotator, ambiguous_hgvs=ambiguous_hgvs, appris=appris, callback=callback, canonical=canonical, ccds=ccds, dbNSFP=dbNSFP, dbscSNV=dbscSNV, distance=distance, domains=domains, failed=failed, hgvs=hgvs, mane=mane, merged=merged, minimal=minimal,
                                                                               mirna=mirna, mutfunc=mutfunc, numbers=numbers, protein=protein, refseq=refseq, shift_3prime=shift_3prime, shift_genomic=shift_genomic, transcript_id=transcript_id, transcript_version=transcript_version, tsl=tsl, uniprot=uniprot, variant_class=variant_class, vcf_string=vcf_string, xref_refseq=xref_refseq), format=format)

    @vep_hgvs.register
    def _(self, hgvs_notation: list, species='human', AncestralAllele=None, Blosum62=None, CADD=None, DisGeNET=None, EVE=None, GO=None, GeneSplicer=None, IntAct=None, LoF=None, Mastermind=None, MaxEntScan=None,
          NMD=None, Phenotypes=None, SpliceAI=None, UTRAnnotator=None, ambiguous_hgvs=None, appris=None, callback=None, canonical=None, ccds=None, dbNSFP=None, dbscSNV=None, distance=None, domains=None, failed=None, hgvs=None, mane=None, merged=None, minimal=None,
          mirna=None, mutfunc=None, numbers=None, protein=None, refseq=None, shift_3prime=None, shift_genomic=None, transcript_id=None, transcript_version=None, tsl=None, uniprot=None, variant_class=None, vcf_string=None, xref_refseq=None, format='json'):
        """Fetch variant consequences for multiple HGVS notations"""
        return post(endpoint=f"vep/{species}/hgvs", params=dict(AncestralAllele=AncestralAllele, Blosum62=Blosum62, CADD=CADD, DisGeNET=DisGeNET, EVE=EVE, GO=GO, GeneSplicer=GeneSplicer, IntAct=IntAct, LoF=LoF, Mastermind=Mastermind, MaxEntScan=MaxEntScan,
                                                                NMD=NMD, Phenotypes=Phenotypes, SpliceAI=SpliceAI, UTRAnnotator=UTRAnnotator, ambiguous_hgvs=ambiguous_hgvs, appris=appris, callback=callback, canonical=canonical, ccds=ccds, dbNSFP=dbNSFP, dbscSNV=dbscSNV, distance=distance, domains=domains, failed=failed, hgvs=hgvs, mane=mane, merged=merged, minimal=minimal,
                                                                mirna=mirna, mutfunc=mutfunc, numbers=numbers, protein=protein, refseq=refseq, shift_3prime=shift_3prime, shift_genomic=shift_genomic, transcript_id=transcript_id, transcript_version=transcript_version, tsl=tsl, uniprot=uniprot, variant_class=variant_class, vcf_string=vcf_string, xref_refseq=xref_refseq), json={"hgvs_notations": hgvs_notation}, format=format)

    @singledispatchmethod
    def vep_id(self, id: str, species='human', AncestralAllele=None, Blosum62=None, CADD=None, Conservation=None, DisGeNET=None, EVE=None, GO=None, GeneSplicer=None, IntAct=None, LoF=None, Mastermind=None, MaxEntScan=None,
               NMD=None, Phenotypes=None, SpliceAI=None, UTRAnnotator=None, appris=None, callback=None, canonical=None, ccds=None, dbNSFP=None, dbscSNV=None, distance=None, domains=None, failed=None, hgvs=None, mane=None, merged=None, minimal=None,
               mirna=None, mutfunc=None, numbers=None, protein=None, refseq=None, shift_3prime=None, shift_genomic=None, transcript_id=None, transcript_version=None, tsl=None, uniprot=None, variant_class=None, vcf_string=None, xref_refseq=None, format='json'):
        """Fetch variant consequences based on a variant identifier"""
        return get(endpoint=f"vep/{species}/id/{id}", params=dict(AncestralAllele=AncestralAllele, Blosum62=Blosum62, CADD=CADD, Conservation=Conservation, DisGeNET=DisGeNET, EVE=EVE, GO=GO, GeneSplicer=GeneSplicer, IntAct=IntAct, LoF=LoF, Mastermind=Mastermind, MaxEntScan=MaxEntScan,
                                                                  NMD=NMD, Phenotypes=Phenotypes, SpliceAI=SpliceAI, UTRAnnotator=UTRAnnotator, appris=appris, callback=callback, canonical=canonical, ccds=ccds, dbNSFP=dbNSFP, dbscSNV=dbscSNV, distance=distance, domains=domains, failed=failed, hgvs=hgvs, mane=mane, merged=merged, minimal=minimal,
                                                                  mirna=mirna, mutfunc=mutfunc, numbers=numbers, protein=protein, refseq=refseq, shift_3prime=shift_3prime, shift_genomic=shift_genomic, transcript_id=transcript_id, transcript_version=transcript_version, tsl=tsl, uniprot=uniprot, variant_class=variant_class, vcf_string=vcf_string, xref_refseq=xref_refseq), format=format)

    @vep_id.register
    def _(self, id: list, species='human', AncestralAllele=None, Blosum62=None, CADD=None, DisGeNET=None, EVE=None, GO=None, GeneSplicer=None, IntAct=None, LoF=None, Mastermind=None, MaxEntScan=None,
          NMD=None, Phenotypes=None, SpliceAI=None, UTRAnnotator=None, appris=None, callback=None, canonical=None, ccds=None, dbNSFP=None, dbscSNV=None, distance=None, domains=None, failed=None, hgvs=None, mane=None, merged=None, minimal=None,
          mirna=None, mutfunc=None, numbers=None, protein=None, refseq=None, shift_3prime=None, shift_genomic=None, transcript_id=None, transcript_version=None, tsl=None, uniprot=None, variant_class=None, vcf_string=None, xref_refseq=None):
        """Fetch variant consequences for multiple ids"""
        return post(endpoint=f"vep/{species}/id", params=dict(AncestralAllele=AncestralAllele, Blosum62=Blosum62, CADD=CADD, DisGeNET=DisGeNET, EVE=EVE, GO=GO, GeneSplicer=GeneSplicer, IntAct=IntAct, LoF=LoF, Mastermind=Mastermind, MaxEntScan=MaxEntScan,
                                                              NMD=NMD, Phenotypes=Phenotypes, SpliceAI=SpliceAI, UTRAnnotator=UTRAnnotator, appris=appris, callback=callback, canonical=canonical, ccds=ccds, dbNSFP=dbNSFP, dbscSNV=dbscSNV, distance=distance, domains=domains, failed=failed, hgvs=hgvs, mane=mane, merged=merged, minimal=minimal,
                                                              mirna=mirna, mutfunc=mutfunc, numbers=numbers, protein=protein, refseq=refseq, shift_3prime=shift_3prime, shift_genomic=shift_genomic, transcript_id=transcript_id, transcript_version=transcript_version, tsl=tsl, uniprot=uniprot, variant_class=variant_class, vcf_string=vcf_string, xref_refseq=xref_refseq), json={"ids": id}, format=format)

    @singledispatchmethod
    def vep_region(region: str, allele: str, species='human', AncestralAllele=None, Blosum62=None, CADD=None, Conservation=None, DisGeNET=None, EVE=None, GO=None, GeneSplicer=None, IntAct=None, LoF=None, Mastermind=None, MaxEntScan=None,
                   NMD=None, Phenotypes=None, SpliceAI=None, UTRAnnotator=None, appris=None, callback=None, canonical=None, ccds=None, dbNSFP=None, dbscSNV=None, distance=None, domains=None, failed=None, hgvs=None, mane=None, merged=None, minimal=None,
                   mirna=None, mutfunc=None, numbers=None, protein=None, refseq=None, shift_3prime=None, shift_genomic=None, transcript_id=None, transcript_version=None, tsl=None, uniprot=None, variant_class=None, vcf_string=None, xref_refseq=None):
        """Fetch variant consequences based on a region"""
        return get(endpoint=f"vep/{species}/region/{region}/{allele}", params=dict(AncestralAllele=AncestralAllele, Blosum62=Blosum62, CADD=CADD, Conservation=Conservation, DisGeNET=DisGeNET, EVE=EVE, GO=GO, GeneSplicer=GeneSplicer, IntAct=IntAct, LoF=LoF, Mastermind=Mastermind, MaxEntScan=MaxEntScan,
                                                                                   NMD=NMD, Phenotypes=Phenotypes, SpliceAI=SpliceAI, UTRAnnotator=UTRAnnotator, appris=appris, callback=callback, canonical=canonical, ccds=ccds, dbNSFP=dbNSFP, dbscSNV=dbscSNV, distance=distance, domains=domains, failed=failed, hgvs=hgvs, mane=mane, merged=merged, minimal=minimal,
                                                                                   mirna=mirna, mutfunc=mutfunc, numbers=numbers, protein=protein, refseq=refseq, shift_3prime=shift_3prime, shift_genomic=shift_genomic, transcript_id=transcript_id, transcript_version=transcript_version, tsl=tsl, uniprot=uniprot, variant_class=variant_class, vcf_string=vcf_string, xref_refseq=xref_refseq), format=format)

    @vep_region.register
    def _(region: list, species='human', AncestralAllele=None, Blosum62=None, CADD=None, DisGeNET=None, EVE=None, GO=None, GeneSplicer=None, IntAct=None, LoF=None, Mastermind=None, MaxEntScan=None,
          NMD=None, Phenotypes=None, SpliceAI=None, UTRAnnotator=None, appris=None, callback=None, canonical=None, ccds=None, dbNSFP=None, dbscSNV=None, distance=None, domains=None, failed=None, hgvs=None, mane=None, merged=None, minimal=None,
          mirna=None, mutfunc=None, numbers=None, protein=None, refseq=None, shift_3prime=None, shift_genomic=None, transcript_id=None, transcript_version=None, tsl=None, uniprot=None, variant_class=None, vcf_string=None, xref_refseq=None):
        """Fetch variant consequences based on a region"""
        return post(endpoint=f"vep/{species}/region", params=dict(AncestralAllele=AncestralAllele, Blosum62=Blosum62, CADD=CADD, DisGeNET=DisGeNET, EVE=EVE, GO=GO, GeneSplicer=GeneSplicer, IntAct=IntAct, LoF=LoF, Mastermind=Mastermind, MaxEntScan=MaxEntScan,
                                                                  NMD=NMD, Phenotypes=Phenotypes, SpliceAI=SpliceAI, UTRAnnotator=UTRAnnotator, appris=appris, callback=callback, canonical=canonical, ccds=ccds, dbNSFP=dbNSFP, dbscSNV=dbscSNV, distance=distance, domains=domains, failed=failed, hgvs=hgvs, mane=mane, merged=merged, minimal=minimal,
                                                                  mirna=mirna, mutfunc=mutfunc, numbers=numbers, protein=protein, refseq=refseq, shift_3prime=shift_3prime, shift_genomic=shift_genomic, transcript_id=transcript_id, transcript_version=transcript_version, tsl=tsl, uniprot=uniprot, variant_class=variant_class, vcf_string=vcf_string, xref_refseq=xref_refseq), json={"variants": region}, format=format)


@singledispatch
def variant_recoder(id: str, species='human', format="json", fields=None, var_synonyms=None, vcf_string=None):
    """Translate a variant identifier, HGVS notation or genomic SPDI notation to all possible variant IDs, HGVS and genomic SPDI"""
    return get(endpoint=f"variant_recoder/{species}/{id}", params=dict(fields=fields, var_synonyms=var_synonyms, vcf_string=vcf_string), format=format)


@variant_recoder.register
def _(id: list, species='human', format="json", fields=None, var_synonyms=None, vcf_string=None):
    """Translate a list of variant identifiers, HGVS notations or genomic SPDI notations to all possible variant IDs, HGVS and genomic SPDI"""
    return post(endpoint=f"variant_recoder/{species}", params=dict(fields=fields, var_synonyms=var_synonyms, vcf_string=vcf_string), json={"ids": id}, format=format)


@singledispatch
def variation(id: str, species='human', format="json", pops=None, genotypes=None, genotyping_chips=None, phenotypes=None, population_genotypes=None):
    """Uses a variant identifier (e.g. rsID) to return the variation features including optional genotype, phenotype and population data"""
    return get(endpoint=f"variation/{species}/{id}", params=dict(pops=pops, genotypes=genotypes, genotyping_chips=genotyping_chips, phenotypes=phenotypes, population_genotypes=population_genotypes), format=format)


@variation.register
def _(id: list, species='human', format="json", pops=None, genotypes=None, phenotypes=None, population_genotypes=None):
    """Uses a list of variant identifiers (e.g. rsID) to return the variation features including optional genotype, phenotype and population data"""
    return post(endpoint=f"variation/{species}", params=dict(pops=pops, genotypes=genotypes, phenotypes=phenotypes, population_genotypes=population_genotypes), json={"ids": id}, format=format)


def variation_pmcid(pmcid, species='human', format="json"):
    """Fetch variants by publication using PubMed Central reference number (PMCID)"""
    return get(endpoint=f"variation/{species}/pmcid/{pmcid}", params=None, format=format)


def variation_pmid(pmid, species='human', format="json"):
    """Fetch variants by publication using PubMed reference number (PMID)"""
    return get(endpoint=f"variation/{species}/pmid/{pmid}", params=None, format=format)


@singledispatch
def vep_hgvs(hgvs_notation: str, species='human', AncestralAllele=None, Blosum62=None, CADD=None, Conservation=None, DisGeNET=None, EVE=None, GO=None, GeneSplicer=None, IntAct=None, LoF=None, Mastermind=None, MaxEntScan=None,
             NMD=None, Phenotypes=None, SpliceAI=None, UTRAnnotator=None, ambiguous_hgvs=None, appris=None, callback=None, canonical=None, ccds=None, dbNSFP=None, dbscSNV=None, distance=None, domains=None, failed=None, hgvs=None, mane=None, merged=None, minimal=None,
             mirna=None, mutfunc=None, numbers=None, protein=None, refseq=None, shift_3prime=None, shift_genomic=None, transcript_id=None, transcript_version=None, tsl=None, uniprot=None, variant_class=None, vcf_string=None, xref_refseq=None, format='json'):
    """Fetch variant consequences based on a HGVS notation"""
    return get(endpoint=f"vep/{species}/hgvs/{hgvs_notation}", params=dict(AncestralAllele=AncestralAllele, Blosum62=Blosum62, CADD=CADD, Conservation=Conservation, DisGeNET=DisGeNET, EVE=EVE, GO=GO, GeneSplicer=GeneSplicer, IntAct=IntAct, LoF=LoF, Mastermind=Mastermind, MaxEntScan=MaxEntScan,
                                                                           NMD=NMD, Phenotypes=Phenotypes, SpliceAI=SpliceAI, UTRAnnotator=UTRAnnotator, ambiguous_hgvs=ambiguous_hgvs, appris=appris, callback=callback, canonical=canonical, ccds=ccds, dbNSFP=dbNSFP, dbscSNV=dbscSNV, distance=distance, domains=domains, failed=failed, hgvs=hgvs, mane=mane, merged=merged, minimal=minimal,
                                                                           mirna=mirna, mutfunc=mutfunc, numbers=numbers, protein=protein, refseq=refseq, shift_3prime=shift_3prime, shift_genomic=shift_genomic, transcript_id=transcript_id, transcript_version=transcript_version, tsl=tsl, uniprot=uniprot, variant_class=variant_class, vcf_string=vcf_string, xref_refseq=xref_refseq), format=format)


@vep_hgvs.register
def _(hgvs_notation: list, species='human', AncestralAllele=None, Blosum62=None, CADD=None, DisGeNET=None, EVE=None, GO=None, GeneSplicer=None, IntAct=None, LoF=None, Mastermind=None, MaxEntScan=None,
      NMD=None, Phenotypes=None, SpliceAI=None, UTRAnnotator=None, ambiguous_hgvs=None, appris=None, callback=None, canonical=None, ccds=None, dbNSFP=None, dbscSNV=None, distance=None, domains=None, failed=None, hgvs=None, mane=None, merged=None, minimal=None,
      mirna=None, mutfunc=None, numbers=None, protein=None, refseq=None, shift_3prime=None, shift_genomic=None, transcript_id=None, transcript_version=None, tsl=None, uniprot=None, variant_class=None, vcf_string=None, xref_refseq=None, format='json'):
    """Fetch variant consequences for multiple HGVS notations"""
    return post(endpoint=f"vep/{species}/hgvs", params=dict(AncestralAllele=AncestralAllele, Blosum62=Blosum62, CADD=CADD, DisGeNET=DisGeNET, EVE=EVE, GO=GO, GeneSplicer=GeneSplicer, IntAct=IntAct, LoF=LoF, Mastermind=Mastermind, MaxEntScan=MaxEntScan,
                                                            NMD=NMD, Phenotypes=Phenotypes, SpliceAI=SpliceAI, UTRAnnotator=UTRAnnotator, ambiguous_hgvs=ambiguous_hgvs, appris=appris, callback=callback, canonical=canonical, ccds=ccds, dbNSFP=dbNSFP, dbscSNV=dbscSNV, distance=distance, domains=domains, failed=failed, hgvs=hgvs, mane=mane, merged=merged, minimal=minimal,
                                                            mirna=mirna, mutfunc=mutfunc, numbers=numbers, protein=protein, refseq=refseq, shift_3prime=shift_3prime, shift_genomic=shift_genomic, transcript_id=transcript_id, transcript_version=transcript_version, tsl=tsl, uniprot=uniprot, variant_class=variant_class, vcf_string=vcf_string, xref_refseq=xref_refseq), json={"hgvs_notations": hgvs_notation}, format=format)


@singledispatch
def vep_id(id: str, species='human', AncestralAllele=None, Blosum62=None, CADD=None, Conservation=None, DisGeNET=None, EVE=None, GO=None, GeneSplicer=None, IntAct=None, LoF=None, Mastermind=None, MaxEntScan=None,
           NMD=None, Phenotypes=None, SpliceAI=None, UTRAnnotator=None, appris=None, callback=None, canonical=None, ccds=None, dbNSFP=None, dbscSNV=None, distance=None, domains=None, failed=None, hgvs=None, mane=None, merged=None, minimal=None,
           mirna=None, mutfunc=None, numbers=None, protein=None, refseq=None, shift_3prime=None, shift_genomic=None, transcript_id=None, transcript_version=None, tsl=None, uniprot=None, variant_class=None, vcf_string=None, xref_refseq=None, format='json'):
    """Fetch variant consequences based on a variant identifier"""
    return get(endpoint=f"vep/{species}/id/{id}", params=dict(AncestralAllele=AncestralAllele, Blosum62=Blosum62, CADD=CADD, Conservation=Conservation, DisGeNET=DisGeNET, EVE=EVE, GO=GO, GeneSplicer=GeneSplicer, IntAct=IntAct, LoF=LoF, Mastermind=Mastermind, MaxEntScan=MaxEntScan,
                                                              NMD=NMD, Phenotypes=Phenotypes, SpliceAI=SpliceAI, UTRAnnotator=UTRAnnotator, appris=appris, callback=callback, canonical=canonical, ccds=ccds, dbNSFP=dbNSFP, dbscSNV=dbscSNV, distance=distance, domains=domains, failed=failed, hgvs=hgvs, mane=mane, merged=merged, minimal=minimal,
                                                              mirna=mirna, mutfunc=mutfunc, numbers=numbers, protein=protein, refseq=refseq, shift_3prime=shift_3prime, shift_genomic=shift_genomic, transcript_id=transcript_id, transcript_version=transcript_version, tsl=tsl, uniprot=uniprot, variant_class=variant_class, vcf_string=vcf_string, xref_refseq=xref_refseq), format=format)


@vep_id.register
def _(id: list, species='human', AncestralAllele=None, Blosum62=None, CADD=None, DisGeNET=None, EVE=None, GO=None, GeneSplicer=None, IntAct=None, LoF=None, Mastermind=None, MaxEntScan=None,
      NMD=None, Phenotypes=None, SpliceAI=None, UTRAnnotator=None, appris=None, callback=None, canonical=None, ccds=None, dbNSFP=None, dbscSNV=None, distance=None, domains=None, failed=None, hgvs=None, mane=None, merged=None, minimal=None,
      mirna=None, mutfunc=None, numbers=None, protein=None, refseq=None, shift_3prime=None, shift_genomic=None, transcript_id=None, transcript_version=None, tsl=None, uniprot=None, variant_class=None, vcf_string=None, xref_refseq=None):
    """Fetch variant consequences for multiple ids"""
    return post(endpoint=f"vep/{species}/id", params=dict(AncestralAllele=AncestralAllele, Blosum62=Blosum62, CADD=CADD, DisGeNET=DisGeNET, EVE=EVE, GO=GO, GeneSplicer=GeneSplicer, IntAct=IntAct, LoF=LoF, Mastermind=Mastermind, MaxEntScan=MaxEntScan,
                                                          NMD=NMD, Phenotypes=Phenotypes, SpliceAI=SpliceAI, UTRAnnotator=UTRAnnotator, appris=appris, callback=callback, canonical=canonical, ccds=ccds, dbNSFP=dbNSFP, dbscSNV=dbscSNV, distance=distance, domains=domains, failed=failed, hgvs=hgvs, mane=mane, merged=merged, minimal=minimal,
                                                          mirna=mirna, mutfunc=mutfunc, numbers=numbers, protein=protein, refseq=refseq, shift_3prime=shift_3prime, shift_genomic=shift_genomic, transcript_id=transcript_id, transcript_version=transcript_version, tsl=tsl, uniprot=uniprot, variant_class=variant_class, vcf_string=vcf_string, xref_refseq=xref_refseq), json={"ids": id}, format=format)


@singledispatch
def vep_region(region: str, allele: str, species='human', AncestralAllele=None, Blosum62=None, CADD=None, Conservation=None, DisGeNET=None, EVE=None, GO=None, GeneSplicer=None, IntAct=None, LoF=None, Mastermind=None, MaxEntScan=None,
               NMD=None, Phenotypes=None, SpliceAI=None, UTRAnnotator=None, appris=None, callback=None, canonical=None, ccds=None, dbNSFP=None, dbscSNV=None, distance=None, domains=None, failed=None, hgvs=None, mane=None, merged=None, minimal=None,
               mirna=None, mutfunc=None, numbers=None, protein=None, refseq=None, shift_3prime=None, shift_genomic=None, transcript_id=None, transcript_version=None, tsl=None, uniprot=None, variant_class=None, vcf_string=None, xref_refseq=None):
    """Fetch variant consequences based on a region"""
    return get(endpoint=f"vep/{species}/region/{region}/{allele}", params=dict(AncestralAllele=AncestralAllele, Blosum62=Blosum62, CADD=CADD, Conservation=Conservation, DisGeNET=DisGeNET, EVE=EVE, GO=GO, GeneSplicer=GeneSplicer, IntAct=IntAct, LoF=LoF, Mastermind=Mastermind, MaxEntScan=MaxEntScan,
                                                                               NMD=NMD, Phenotypes=Phenotypes, SpliceAI=SpliceAI, UTRAnnotator=UTRAnnotator, appris=appris, callback=callback, canonical=canonical, ccds=ccds, dbNSFP=dbNSFP, dbscSNV=dbscSNV, distance=distance, domains=domains, failed=failed, hgvs=hgvs, mane=mane, merged=merged, minimal=minimal,
                                                                               mirna=mirna, mutfunc=mutfunc, numbers=numbers, protein=protein, refseq=refseq, shift_3prime=shift_3prime, shift_genomic=shift_genomic, transcript_id=transcript_id, transcript_version=transcript_version, tsl=tsl, uniprot=uniprot, variant_class=variant_class, vcf_string=vcf_string, xref_refseq=xref_refseq), format=format)


@vep_region.register
def _(region: list, species='human', AncestralAllele=None, Blosum62=None, CADD=None, DisGeNET=None, EVE=None, GO=None, GeneSplicer=None, IntAct=None, LoF=None, Mastermind=None, MaxEntScan=None,
      NMD=None, Phenotypes=None, SpliceAI=None, UTRAnnotator=None, appris=None, callback=None, canonical=None, ccds=None, dbNSFP=None, dbscSNV=None, distance=None, domains=None, failed=None, hgvs=None, mane=None, merged=None, minimal=None,
      mirna=None, mutfunc=None, numbers=None, protein=None, refseq=None, shift_3prime=None, shift_genomic=None, transcript_id=None, transcript_version=None, tsl=None, uniprot=None, variant_class=None, vcf_string=None, xref_refseq=None):
    """Fetch variant consequences based on a region"""
    return post(endpoint=f"vep/{species}/region", params=dict(AncestralAllele=AncestralAllele, Blosum62=Blosum62, CADD=CADD, DisGeNET=DisGeNET, EVE=EVE, GO=GO, GeneSplicer=GeneSplicer, IntAct=IntAct, LoF=LoF, Mastermind=Mastermind, MaxEntScan=MaxEntScan,
                                                              NMD=NMD, Phenotypes=Phenotypes, SpliceAI=SpliceAI, UTRAnnotator=UTRAnnotator, appris=appris, callback=callback, canonical=canonical, ccds=ccds, dbNSFP=dbNSFP, dbscSNV=dbscSNV, distance=distance, domains=domains, failed=failed, hgvs=hgvs, mane=mane, merged=merged, minimal=minimal,
                                                              mirna=mirna, mutfunc=mutfunc, numbers=numbers, protein=protein, refseq=refseq, shift_3prime=shift_3prime, shift_genomic=shift_genomic, transcript_id=transcript_id, transcript_version=transcript_version, tsl=tsl, uniprot=uniprot, variant_class=variant_class, vcf_string=vcf_string, xref_refseq=xref_refseq), json={"variants": region}, format=format)
