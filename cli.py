import argparse
import pprint
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


@singledispatch
def variant_recoder(id: str, species='human', format="json", fields=None, var_synonyms=None, vcf_string=None):
    """Translate a variant identifier, HGVS notation or genomic SPDI notation to all possible variant IDs, HGVS and genomic SPDI"""
    return get(endpoint=f"variant_recoder/{species}/{id}", params=dict(fields=fields, var_synonyms=var_synonyms, vcf_string=vcf_string), format=format)


@variant_recoder.register
def _(id: list, species='human', format="json", fields=None, var_synonyms=None, vcf_string=None):
    """Translate a list of variant identifiers, HGVS notations or genomic SPDI notations to all possible variant IDs, HGVS and genomic SPDI"""
    return post(endpoint=f"variant_recoder/{species}", params=dict(fields=fields, var_synonyms=var_synonyms, vcf_string=vcf_string), json={"ids": id}, format=format)


@singledispatch
def variation(id: str, species='human', format="json", pops=False, genotypes=False, genotyping_chips=False, phenotypes=False, population_genotypes=False):
    """Uses a variant identifier (e.g. rsID) to return the variation features including optional genotype, phenotype and population data"""
    return get(endpoint=f"variation/{species}/{id}", params=dict(pops=pops, genotypes=genotypes, genotyping_chips=genotyping_chips, phenotypes=phenotypes, population_genotypes=population_genotypes), format=format)


@variation.register
def _(id: list, species='human', format="json", pops=False, genotypes=False, phenotypes=False, population_genotypes=False):
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


parser = argparse.ArgumentParser(description='Ensembl REST API client', fromfile_prefix_chars='@',
                                 prog="ensembl_rest", epilog="For further help on a command, type ensembl_rest <command> -h")
parser.add_argument('--assembly', default="GRCh38",
                    help='Genome assembly', choices=["GRCh38", "GRCh37"])
parser.add_argument('--scheme', default="https",
                    help='Scheme', choices=["http", "https"])
subparsers = parser.add_subparsers(dest="command")
# ------------------------------------------------------------------------------------------------------------------------
# variant_recoder_get
# ------------------------------------------------------------------------------------------------------------------------
parser_variant_recoder_get = subparsers.add_parser(
    'variant_recoder_get', help='variant_recoder help')
parser_variant_recoder_get.add_argument(
    'id', help='Variant identifier, HGVS notation or genomic SPDI notation')
parser_variant_recoder_get.add_argument(
    '--species', default="human", help='Species')
parser_variant_recoder_get.add_argument(
    "--vcf_string", help="VCF string", action="store_const", const=1)
parser_variant_recoder_get.add_argument('--fields', help='Fields')
parser_variant_recoder_get.add_argument(
    '--var_synonyms', help='Var synonyms', action="store_const", const=1)
parser_variant_recoder_get.add_argument(
    '--format', default="json", help='Format', choices=["json", "xml", "jsonp"])
# ------------------------------------------------------------------------------------------------------------------------
# variant_recoder_post
# -------------------------------------------------------------------------------------------------
parser_variant_recoder_post = subparsers.add_parser(
    'variant_recoder_post', help='variant_recoder help')
parser_variant_recoder_post.add_argument(
    'id', help='Variant identifier, HGVS notation or genomic SPDI notation', action="extend", nargs="+")
parser_variant_recoder_post.add_argument(
    '--species', default="human", help='Species')
parser_variant_recoder_post.add_argument(
    "--vcf_string", help="VCF string", action="store_const", const=1)
parser_variant_recoder_post.add_argument('--fields', help='Fields')
parser_variant_recoder_post.add_argument(
    '--var_synonyms', help='Var synonyms', action="store_const", const=1)
parser_variant_recoder_post.add_argument(
    '--format', default="json", help='Format', choices=["json", "xml", "jsonp"])
# -------------------------------------------------------------------------------------------------
# variation_get
# -------------------------------------------------------------------
parser_variation_get = subparsers.add_parser(
    'variation_get', help='variation help')
parser_variation_get.add_argument(
    'id', help='Variant identifier')
parser_variation_get.add_argument(
    '--species', default="human", help='Species')
parser_variation_get.add_argument('--format', default="json", help='Format')
parser_variation_get.add_argument(
    '--phenotypes', action="store_const", const=1, help='Phenotypes')
parser_variation_get.add_argument(
    '--pops', help='Pops', action="store_const", const=1)
parser_variation_get.add_argument(
    '--genotypes', help='Genotypes', action="store_const", const=1)
parser_variation_get.add_argument(
    '--genotyping_chips', help='Genotyping chips', action="store_const", const=1)
parser_variation_get.add_argument(
    '--population_genotypes', help='Population genotypes', action="store_const", const=1)
# -------------------------------------------------------------------
# variation_post
# -------------------------------------------------------------------
parser_variation_post = subparsers.add_parser(
    'variation_post', help='variation help')
parser_variation_post.add_argument(
    'id', help='Variant identifier', action="extend", nargs="+")
parser_variation_post.add_argument(
    '--species', default="human", help='Species')
parser_variation_post.add_argument('--format', default="json", help='Format')
parser_variation_post.add_argument(
    '--phenotypes', action="store_const", const=1, help='Phenotypes')
parser_variation_post.add_argument(
    '--pops', help='Pops', action="store_const", const=1)
parser_variation_post.add_argument(
    '--genotypes', help='Genotypes', action="store_const", const=1)
parser_variation_post.add_argument(
    '--population_genotypes', help='Population genotypes', action="store_const", const=1)
# -------------------------------------------------------------------
# variation_pmcid
# -------------------------------------------------------------------
parser_variation_pmcid = subparsers.add_parser("variation_pmcid")
parser_variation_pmcid.add_argument('pmcid', help='PMCID')
parser_variation_pmcid.add_argument(
    '--species', default="human", help='Species')
parser_variation_pmcid.add_argument('--format', default="json", help='Format')
# -------------------------------------------------------------------
# variation_pmid
# -------------------------------------------------------------------
parser_variation_pmid = subparsers.add_parser("variation_pmid")
parser_variation_pmid.add_argument('pmid', help='PMID')
parser_variation_pmid.add_argument(
    '--species', default="human", help='Species')
parser_variation_pmid.add_argument('--format', default="json", help='Format')
# -------------------------------------------------------------------
# vep_hgvs_get
# -------------------------------------------------------------------
parser_vep_hgvs_get = subparsers.add_parser("vep_hgvs_get")
parser_vep_hgvs_get.add_argument('hgvs_notation', help='HGVS notation')
parser_vep_hgvs_get.add_argument('--species', default="human", help='Species')
parser_vep_hgvs_get.add_argument('--format', default="json", help='Format')
parser_vep_hgvs_get.add_argument(
    '--AncestralAllele', action="store_const", const=1)
parser_vep_hgvs_get.add_argument('--Blosum62', action="store_const", const=1)
parser_vep_hgvs_get.add_argument('--CADD', action="store_const", const=1)
parser_vep_hgvs_get.add_argument(
    '--Conservation', action="store_const", const=1)
parser_vep_hgvs_get.add_argument('--DisGeNET', action="store_const", const=1)
parser_vep_hgvs_get.add_argument('--EVE', action="store_const", const=1)
parser_vep_hgvs_get.add_argument('--GO', action="store_const", const=1)
parser_vep_hgvs_get.add_argument(
    '--GeneSplicer', action="store_const", const=1)
parser_vep_hgvs_get.add_argument('--IntAct', action="store_const", const=1)
parser_vep_hgvs_get.add_argument('--LoF', action="store_const", const=1)
parser_vep_hgvs_get.add_argument('--Mastermind', action="store_const", const=1)
parser_vep_hgvs_get.add_argument('--MaxEntScan', action="store_const", const=1)
parser_vep_hgvs_get.add_argument('--NMD', action="store_const", const=1)
parser_vep_hgvs_get.add_argument('--Phenotypes', action="store_const", const=1)
parser_vep_hgvs_get.add_argument('--SpliceAI', default=0, type=int)
parser_vep_hgvs_get.add_argument(
    '--UTRAnnotator', action="store_const", const=1)
parser_vep_hgvs_get.add_argument(
    '--ambiguous_hgvs', action="store_const", const=1)
parser_vep_hgvs_get.add_argument('--appris', action="store_const", const=1)
parser_vep_hgvs_get.add_argument('--callback')
parser_vep_hgvs_get.add_argument('--canonical', action="store_const", const=1)
parser_vep_hgvs_get.add_argument('--ccds', action="store_const", const=1)
parser_vep_hgvs_get.add_argument('--dbNSFP')
parser_vep_hgvs_get.add_argument('--dbscSNV', action="store_const", const=1)
parser_vep_hgvs_get.add_argument('--distance', type=int)
parser_vep_hgvs_get.add_argument('--domains', action="store_const", const=1)
parser_vep_hgvs_get.add_argument('--failed', action="store_const", const=1)
parser_vep_hgvs_get.add_argument('--hgvs', action="store_const", const=1)
parser_vep_hgvs_get.add_argument('--mane', action="store_const", const=1)
parser_vep_hgvs_get.add_argument('--merged', action="store_const", const=1)
parser_vep_hgvs_get.add_argument('--minimal', action="store_const", const=1)
parser_vep_hgvs_get.add_argument('--mirna', action="store_const", const=1)
parser_vep_hgvs_get.add_argument('--mutfunc', action="store_const", const=1)
parser_vep_hgvs_get.add_argument('--numbers', action="store_const", const=1)
parser_vep_hgvs_get.add_argument('--protein', action="store_const", const=1)
parser_vep_hgvs_get.add_argument('--refseq', action="store_const", const=1)
parser_vep_hgvs_get.add_argument(
    '--shift_3prime', action="store_const", const=1)
parser_vep_hgvs_get.add_argument(
    '--shift_genomic', action="store_const", const=1)
parser_vep_hgvs_get.add_argument('--transcript_id')
parser_vep_hgvs_get.add_argument(
    '--transcript_version', action="store_const", const=1)
parser_vep_hgvs_get.add_argument('--tsl', action="store_const", const=1)
parser_vep_hgvs_get.add_argument('--uniprot', action="store_const", const=1)
parser_vep_hgvs_get.add_argument(
    '--variant_class', action="store_const", const=1)
parser_vep_hgvs_get.add_argument('--vcf_string', action="store_const", const=1)
parser_vep_hgvs_get.add_argument(
    '--xref_refseq', action="store_const", const=1)
# -----------------------------------------------------------------------------------------------------------------------
# vep_hgvs_post
# -----------------------------------------------------------------------------------------------------------------------
parser_vep_hgvs_post = subparsers.add_parser("vep_hgvs_post")
parser_vep_hgvs_post.add_argument(
    'hgvs_notation', help='HGVS notation', action="extend", nargs="+")
parser_vep_hgvs_post.add_argument('--species', default="human", help='Species')
parser_vep_hgvs_post.add_argument('--format', default="json", help='Format')
parser_vep_hgvs_post.add_argument(
    '--AncestralAllele', action="store_const", const=1)
parser_vep_hgvs_post.add_argument('--Blosum62', action="store_const", const=1)
parser_vep_hgvs_post.add_argument('--CADD', action="store_const", const=1)
parser_vep_hgvs_post.add_argument('--DisGeNET', action="store_const", const=1)
parser_vep_hgvs_post.add_argument('--EVE', action="store_const", const=1)
parser_vep_hgvs_post.add_argument('--GO', action="store_const", const=1)
parser_vep_hgvs_post.add_argument(
    '--GeneSplicer', action="store_const", const=1)
parser_vep_hgvs_post.add_argument('--IntAct', action="store_const", const=1)
parser_vep_hgvs_post.add_argument('--LoF', action="store_const", const=1)
parser_vep_hgvs_post.add_argument(
    '--Mastermind', action="store_const", const=1)
parser_vep_hgvs_post.add_argument(
    '--MaxEntScan', action="store_const", const=1)
parser_vep_hgvs_post.add_argument('--NMD', action="store_const", const=1)
parser_vep_hgvs_post.add_argument(
    '--Phenotypes', action="store_const", const=1)
parser_vep_hgvs_post.add_argument('--SpliceAI', default=0, type=int)
parser_vep_hgvs_post.add_argument(
    '--UTRAnnotator', action="store_const", const=1)
parser_vep_hgvs_post.add_argument(
    '--ambiguous_hgvs', action="store_const", const=1)
parser_vep_hgvs_post.add_argument('--appris', action="store_const", const=1)
parser_vep_hgvs_post.add_argument('--callback')
parser_vep_hgvs_post.add_argument('--canonical', action="store_const", const=1)
parser_vep_hgvs_post.add_argument('--ccds', action="store_const", const=1)
parser_vep_hgvs_post.add_argument('--dbNSFP')
parser_vep_hgvs_post.add_argument('--dbscSNV', action="store_const", const=1)
parser_vep_hgvs_post.add_argument('--distance', type=int)
parser_vep_hgvs_post.add_argument('--domains', action="store_const", const=1)
parser_vep_hgvs_post.add_argument('--failed', action="store_const", const=1)
parser_vep_hgvs_post.add_argument('--hgvs', action="store_const", const=1)
parser_vep_hgvs_post.add_argument('--mane', action="store_const", const=1)
parser_vep_hgvs_post.add_argument('--merged', action="store_const", const=1)
parser_vep_hgvs_post.add_argument('--minimal', action="store_const", const=1)
parser_vep_hgvs_post.add_argument('--mirna', action="store_const", const=1)
parser_vep_hgvs_post.add_argument('--mutfunc', action="store_const", const=1)
parser_vep_hgvs_post.add_argument('--numbers', action="store_const", const=1)
parser_vep_hgvs_post.add_argument('--protein', action="store_const", const=1)
parser_vep_hgvs_post.add_argument('--refseq', action="store_const", const=1)
parser_vep_hgvs_post.add_argument(
    '--shift_3prime', action="store_const", const=1)
parser_vep_hgvs_post.add_argument(
    '--shift_genomic', action="store_const", const=1)
parser_vep_hgvs_post.add_argument('--transcript_id')
parser_vep_hgvs_post.add_argument(
    '--transcript_version', action="store_const", const=1)
parser_vep_hgvs_post.add_argument('--tsl', action="store_const", const=1)
parser_vep_hgvs_post.add_argument('--uniprot', action="store_const", const=1)
parser_vep_hgvs_post.add_argument(
    '--variant_class', action="store_const", const=1)
parser_vep_hgvs_post.add_argument(
    '--vcf_string', action="store_const", const=1)
parser_vep_hgvs_post.add_argument(
    '--xref_refseq', action="store_const", const=1)
# ---------------------------------------------------------------------------------------------
# vep_id_get
# ---------------------------------------------------------------------------------------------
parser_vep_id_get = subparsers.add_parser("vep_id_get")
parser_vep_id_get.add_argument(
    'id', help='Variant identifier')
parser_vep_id_get.add_argument('--species', default="human", help='Species')
parser_vep_id_get.add_argument('--format', default="json", help='Format')
parser_vep_id_get.add_argument(
    '--AncestralAllele', action="store_const", const=1)
parser_vep_id_get.add_argument('--Blosum62', action="store_const", const=1)
parser_vep_id_get.add_argument('--CADD', action="store_const", const=1)
parser_vep_id_get.add_argument('--Conservation', action="store_const", const=1)
parser_vep_id_get.add_argument('--DisGeNET', action="store_const", const=1)
parser_vep_id_get.add_argument('--EVE', action="store_const", const=1)
parser_vep_id_get.add_argument('--GO', action="store_const", const=1)
parser_vep_id_get.add_argument('--GeneSplicer', action="store_const", const=1)
parser_vep_id_get.add_argument('--IntAct', action="store_const", const=1)
parser_vep_id_get.add_argument('--LoF', action="store_const", const=1)
parser_vep_id_get.add_argument('--Mastermind', action="store_const", const=1)
parser_vep_id_get.add_argument('--MaxEntScan', action="store_const", const=1)
parser_vep_id_get.add_argument('--NMD', action="store_const", const=1)
parser_vep_id_get.add_argument('--Phenotypes', action="store_const", const=1)
parser_vep_id_get.add_argument('--SpliceAI', type=int)
parser_vep_id_get.add_argument('--UTRAnnotator', action="store_const", const=1)
parser_vep_id_get.add_argument('--appris', action="store_const", const=1)
parser_vep_id_get.add_argument('--callback')
parser_vep_id_get.add_argument('--canonical', action="store_const", const=1)
parser_vep_id_get.add_argument('--ccds', action="store_const", const=1)
parser_vep_id_get.add_argument('--dbNSFP')
parser_vep_id_get.add_argument('--dbscSNV', action="store_const", const=1)
parser_vep_id_get.add_argument('--distance', type=int)
parser_vep_id_get.add_argument('--domains', action="store_const", const=1)
parser_vep_id_get.add_argument('--failed', action="store_const", const=1)
parser_vep_id_get.add_argument('--hgvs', action="store_const", const=1)
parser_vep_id_get.add_argument('--mane', action="store_const", const=1)
parser_vep_id_get.add_argument('--merged', action="store_const", const=1)
parser_vep_id_get.add_argument('--minimal', action="store_const", const=1)
parser_vep_id_get.add_argument('--mirna', action="store_const", const=1)
parser_vep_id_get.add_argument('--mutfunc', action="store_const", const=1)
parser_vep_id_get.add_argument('--numbers', action="store_const", const=1)
parser_vep_id_get.add_argument('--protein', action="store_const", const=1)
parser_vep_id_get.add_argument('--refseq', action="store_const", const=1)
parser_vep_id_get.add_argument('--shift_3prime', action="store_const", const=1)
parser_vep_id_get.add_argument(
    '--shift_genomic', action="store_const", const=1)
parser_vep_id_get.add_argument('--transcript_id')
parser_vep_id_get.add_argument('--transcript_version',
                               action="store_const", const=1)
parser_vep_id_get.add_argument('--tsl', action="store_const", const=1)
parser_vep_id_get.add_argument('--uniprot', action="store_const", const=1)
parser_vep_id_get.add_argument(
    '--variant_class', action="store_const", const=1)
parser_vep_id_get.add_argument('--vcf_string', action="store_const", const=1)
parser_vep_id_get.add_argument('--xref_refseq', action="store_const", const=1)
# ----------------------------------------------------------------------------------------------
# vep_id_post
# ---------------------------------------------------------------------------------------------
parser_vep_id_post = subparsers.add_parser("vep_id_post")
parser_vep_id_post.add_argument(
    'id', help='Variant identifier', action="extend", nargs="+")
parser_vep_id_post.add_argument('--species', default="human", help='Species')
parser_vep_id_post.add_argument('--format', default="json", help='Format')
parser_vep_id_post.add_argument(
    '--AncestralAllele', action="store_const", const=1)
parser_vep_id_post.add_argument('--Blosum62', action="store_const", const=1)
parser_vep_id_post.add_argument('--CADD', action="store_const", const=1)
parser_vep_id_post.add_argument('--DisGeNET', action="store_const", const=1)
parser_vep_id_post.add_argument('--EVE', action="store_const", const=1)
parser_vep_id_post.add_argument('--GO', action="store_const", const=1)
parser_vep_id_post.add_argument('--GeneSplicer', action="store_const", const=1)
parser_vep_id_post.add_argument('--IntAct', action="store_const", const=1)
parser_vep_id_post.add_argument('--LoF', action="store_const", const=1)
parser_vep_id_post.add_argument('--Mastermind', action="store_const", const=1)
parser_vep_id_post.add_argument('--MaxEntScan', action="store_const", const=1)
parser_vep_id_post.add_argument('--NMD', action="store_const", const=1)
parser_vep_id_post.add_argument('--Phenotypes', action="store_const", const=1)
parser_vep_id_post.add_argument('--SpliceAI', type=int)
parser_vep_id_post.add_argument(
    '--UTRAnnotator', action="store_const", const=1)
parser_vep_id_post.add_argument('--appris', action="store_const", const=1)
parser_vep_id_post.add_argument('--callback')
parser_vep_id_post.add_argument('--canonical', action="store_const", const=1)
parser_vep_id_post.add_argument('--ccds', action="store_const", const=1)
parser_vep_id_post.add_argument('--dbNSFP')
parser_vep_id_post.add_argument('--dbscSNV', action="store_const", const=1)
parser_vep_id_post.add_argument('--distance', type=int)
parser_vep_id_post.add_argument('--domains', action="store_const", const=1)
parser_vep_id_post.add_argument('--failed', action="store_const", const=1)
parser_vep_id_post.add_argument('--hgvs', action="store_const", const=1)
parser_vep_id_post.add_argument('--mane', action="store_const", const=1)
parser_vep_id_post.add_argument('--merged', action="store_const", const=1)
parser_vep_id_post.add_argument('--minimal', action="store_const", const=1)
parser_vep_id_post.add_argument('--mirna', action="store_const", const=1)
parser_vep_id_post.add_argument('--mutfunc', action="store_const", const=1)
parser_vep_id_post.add_argument('--numbers', action="store_const", const=1)
parser_vep_id_post.add_argument('--protein', action="store_const", const=1)
parser_vep_id_post.add_argument('--refseq', action="store_const", const=1)
parser_vep_id_post.add_argument(
    '--shift_3prime', action="store_const", const=1)
parser_vep_id_post.add_argument(
    '--shift_genomic', action="store_const", const=1)
parser_vep_id_post.add_argument('--transcript_id')
parser_vep_id_post.add_argument('--transcript_version',
                                action="store_const", const=1)
parser_vep_id_post.add_argument('--tsl', action="store_const", const=1)
parser_vep_id_post.add_argument('--uniprot', action="store_const", const=1)
parser_vep_id_post.add_argument(
    '--variant_class', action="store_const", const=1)
parser_vep_id_post.add_argument('--vcf_string', action="store_const", const=1)
parser_vep_id_post.add_argument('--xref_refseq', action="store_const", const=1)
# ---------------------------------------------------------------------------------------------
# vep_region_get
# ---------------------------------------------------------------------------------------------
parser_vep_region_get = subparsers.add_parser("vep_region_get")
parser_vep_region_get.add_argument('region', help='Ensembl VEP format')
parser_vep_region_get.add_argument('allele', help='Allele')
parser_vep_region_get.add_argument(
    '--species', default="human", help='Species')
parser_vep_region_get.add_argument('--format', default="json", help='Format')
parser_vep_region_get.add_argument(
    '--AncestralAllele', action="store_const", const=1)
parser_vep_region_get.add_argument('--Blosum62', action="store_const", const=1)
parser_vep_region_get.add_argument('--CADD', action="store_const", const=1)
parser_vep_region_get.add_argument(
    '--Conservation', action="store_const", const=1)
parser_vep_region_get.add_argument('--DisGeNET', action="store_const", const=1)
parser_vep_region_get.add_argument('--EVE', action="store_const", const=1)
parser_vep_region_get.add_argument('--GO', action="store_const", const=1)
parser_vep_region_get.add_argument(
    '--GeneSplicer', action="store_const", const=1)
parser_vep_region_get.add_argument('--IntAct', action="store_const", const=1)
parser_vep_region_get.add_argument('--LoF', action="store_const", const=1)
parser_vep_region_get.add_argument(
    '--Mastermind', action="store_const", const=1)
parser_vep_region_get.add_argument(
    '--MaxEntScan', action="store_const", const=1)
parser_vep_region_get.add_argument('--NMD', action="store_const", const=1)
parser_vep_region_get.add_argument(
    '--Phenotypes', action="store_const", const=1)
parser_vep_region_get.add_argument('--SpliceAI', type=int)
parser_vep_region_get.add_argument(
    '--UTRAnnotator', action="store_const", const=1)
parser_vep_region_get.add_argument('--appris', action="store_const", const=1)
parser_vep_region_get.add_argument('--callback')
parser_vep_region_get.add_argument(
    '--canonical', action="store_const", const=1)
parser_vep_region_get.add_argument('--ccds', action="store_const", const=1)
parser_vep_region_get.add_argument('--dbNSFP')
parser_vep_region_get.add_argument('--dbscSNV', action="store_const", const=1)
parser_vep_region_get.add_argument('--distance', type=int)
parser_vep_region_get.add_argument('--domains', action="store_const", const=1)
parser_vep_region_get.add_argument('--failed', action="store_const", const=1)
parser_vep_region_get.add_argument('--hgvs', action="store_const", const=1)
parser_vep_region_get.add_argument('--mane', action="store_const", const=1)
parser_vep_region_get.add_argument('--merged', action="store_const", const=1)
parser_vep_region_get.add_argument('--minimal', action="store_const", const=1)
parser_vep_region_get.add_argument('--mirna', action="store_const", const=1)
parser_vep_region_get.add_argument('--mutfunc', action="store_const", const=1)
parser_vep_region_get.add_argument('--numbers', action="store_const", const=1)
parser_vep_region_get.add_argument('--protein', action="store_const", const=1)
parser_vep_region_get.add_argument('--refseq', action="store_const", const=1)
parser_vep_region_get.add_argument(
    '--shift_3prime', action="store_const", const=1)
parser_vep_region_get.add_argument(
    '--shift_genomic', action="store_const", const=1)
parser_vep_region_get.add_argument('--transcript_id')
parser_vep_region_get.add_argument('--transcript_version',
                                   action="store_const", const=1)
parser_vep_region_get.add_argument('--tsl', action="store_const", const=1)
parser_vep_region_get.add_argument('--uniprot', action="store_const", const=1)
parser_vep_region_get.add_argument(
    '--variant_class', action="store_const", const=1)
parser_vep_region_get.add_argument(
    '--vcf_string', action="store_const", const=1)
parser_vep_region_get.add_argument(
    '--xref_refseq', action="store_const", const=1)
# ---------------------------------------------------------------------------------------------
# vep_region_post
# ---------------------------------------------------------------------------------------------
parser_vep_region_post = subparsers.add_parser("vep_region_post")
parser_vep_region_post.add_argument(
    'region', help='VCF format', nargs='+', action='extend')
parser_vep_region_post.add_argument(
    '--species', default="human", help='Species')
parser_vep_region_post.add_argument('--format', default="json", help='Format')
parser_vep_region_post.add_argument(
    '--AncestralAllele', action="store_const", const=1)
parser_vep_region_post.add_argument(
    '--Blosum62', action="store_const", const=1)
parser_vep_region_post.add_argument('--CADD', action="store_const", const=1)
parser_vep_region_post.add_argument(
    '--DisGeNET', action="store_const", const=1)
parser_vep_region_post.add_argument('--EVE', action="store_const", const=1)
parser_vep_region_post.add_argument('--GO', action="store_const", const=1)
parser_vep_region_post.add_argument(
    '--GeneSplicer', action="store_const", const=1)
parser_vep_region_post.add_argument('--IntAct', action="store_const", const=1)
parser_vep_region_post.add_argument('--LoF', action="store_const", const=1)
parser_vep_region_post.add_argument(
    '--Mastermind', action="store_const", const=1)
parser_vep_region_post.add_argument(
    '--MaxEntScan', action="store_const", const=1)
parser_vep_region_post.add_argument('--NMD', action="store_const", const=1)
parser_vep_region_post.add_argument(
    '--Phenotypes', action="store_const", const=1)
parser_vep_region_post.add_argument('--SpliceAI', type=int)
parser_vep_region_post.add_argument(
    '--UTRAnnotator', action="store_const", const=1)
parser_vep_region_post.add_argument('--appris', action="store_const", const=1)
parser_vep_region_post.add_argument('--callback')
parser_vep_region_post.add_argument(
    '--canonical', action="store_const", const=1)
parser_vep_region_post.add_argument('--ccds', action="store_const", const=1)
parser_vep_region_post.add_argument('--dbNSFP')
parser_vep_region_post.add_argument('--dbscSNV', action="store_const", const=1)
parser_vep_region_post.add_argument('--distance', type=int)
parser_vep_region_post.add_argument('--domains', action="store_const", const=1)
parser_vep_region_post.add_argument('--failed', action="store_const", const=1)
parser_vep_region_post.add_argument('--hgvs', action="store_const", const=1)
parser_vep_region_post.add_argument('--mane', action="store_const", const=1)
parser_vep_region_post.add_argument('--merged', action="store_const", const=1)
parser_vep_region_post.add_argument('--minimal', action="store_const", const=1)
parser_vep_region_post.add_argument('--mirna', action="store_const", const=1)
parser_vep_region_post.add_argument('--mutfunc', action="store_const", const=1)
parser_vep_region_post.add_argument('--numbers', action="store_const", const=1)
parser_vep_region_post.add_argument('--protein', action="store_const", const=1)
parser_vep_region_post.add_argument('--refseq', action="store_const", const=1)
parser_vep_region_post.add_argument(
    '--shift_3prime', action="store_const", const=1)
parser_vep_region_post.add_argument(
    '--shift_genomic', action="store_const", const=1)
parser_vep_region_post.add_argument('--transcript_id')
parser_vep_region_post.add_argument('--transcript_version',
                                    action="store_const", const=1)
parser_vep_region_post.add_argument('--tsl', action="store_const", const=1)
parser_vep_region_post.add_argument('--uniprot', action="store_const", const=1)
parser_vep_region_post.add_argument(
    '--variant_class', action="store_const", const=1)
parser_vep_region_post.add_argument(
    '--vcf_string', action="store_const", const=1)
parser_vep_region_post.add_argument(
    '--xref_refseq', action="store_const", const=1)
# ---------------------------------------------------------------------------------------------
parser_variation_get.set_defaults(func=variation)
parser_variation_post.set_defaults(func=variation)
parser_variation_pmcid.set_defaults(func=variation_pmcid)
parser_variation_pmid.set_defaults(func=variation_pmid)
parser_vep_hgvs_post.set_defaults(func=vep_hgvs)
parser_vep_hgvs_get.set_defaults(func=vep_hgvs)
parser_vep_id_get.set_defaults(func=vep_id)
parser_vep_id_post.set_defaults(func=vep_id)
parser_vep_region_get.set_defaults(func=vep_region)
parser_vep_region_post.set_defaults(func=vep_region)
parser_variant_recoder_get.set_defaults(func=variant_recoder)
parser_variant_recoder_post.set_defaults(func=variant_recoder)
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
    case "variant_recoder_get":
        pprint.pprint(args.func(args.id, species=args.species, fields=args.fields,
                      var_synonyms=args.var_synonyms, vcf_string=args.vcf_string, format=args.format))
    case "variant_recoder_post":
        pprint.pprint(args.func(args.id, species=args.species, fields=args.fields,
                      var_synonyms=args.var_synonyms, vcf_string=args.vcf_string, format=args.format))
    case "variation_get":
        pprint.pprint(args.func(args.id, species=args.species, pops=args.pops, genotypes=args.genotypes,
                      genotyping_chips=args.genotyping_chips, phenotypes=args.phenotypes, population_genotypes=args.population_genotypes))
    case "variation_post":
        pprint.pprint(args.func(args.id, species=args.species, pops=args.pops, genotypes=args.genotypes,
                      phenotypes=args.phenotypes, population_genotypes=args.population_genotypes))
    case "variation_pmcid":
        pprint.pprint(
            args.func(args.pmcid, species=args.species, format=args.format))
    case "variation_pmid":
        pprint.pprint(
            args.func(args.pmid, species=args.species, format=args.format))
    case "vep_hgvs_post":
        pprint.pprint(args.func(args.hgvs_notation, species=args.species, format=args.format, AncestralAllele=args.AncestralAllele, Blosum62=args.Blosum62, CADD=args.CADD, DisGeNET=args.DisGeNET, EVE=args.EVE, GO=args.GO, GeneSplicer=args.GeneSplicer, IntAct=args.IntAct, LoF=args.LoF, Mastermind=args.Mastermind, MaxEntScan=args.MaxEntScan,
                                NMD=args.NMD, Phenotypes=args.Phenotypes, SpliceAI=args.SpliceAI, UTRAnnotator=args.UTRAnnotator, ambiguous_hgvs=args.ambiguous_hgvs, appris=args.appris, callback=args.callback, canonical=args.canonical, ccds=args.ccds, dbNSFP=args.dbNSFP, dbscSNV=args.dbscSNV, distance=args.distance, domains=args.domains, failed=args.failed, hgvs=args.hgvs, mane=args.mane, merged=args.merged, minimal=args.minimal,
                                mirna=args.mirna, mutfunc=args.mutfunc, numbers=args.numbers, protein=args.protein, refseq=args.refseq, shift_3prime=args.shift_3prime, shift_genomic=args.shift_genomic, transcript_id=args.transcript_id, transcript_version=args.transcript_version, tsl=args.tsl, uniprot=args.uniprot, variant_class=args.variant_class, vcf_string=args.vcf_string, xref_refseq=args.xref_refseq))
    case "vep_id_post":
        pprint.pprint(args.func(args.id, species=args.species, format=args.format, AncestralAllele=args.AncestralAllele, Blosum62=args.Blosum62, CADD=args.CADD, DisGeNET=args.DisGeNET, EVE=args.EVE, GO=args.GO, GeneSplicer=args.GeneSplicer, IntAct=args.IntAct, LoF=args.LoF, Mastermind=args.Mastermind, MaxEntScan=args.MaxEntScan,
                      NMD=args.NMD, Phenotypes=args.Phenotypes, SpliceAI=args.SpliceAI, UTRAnnotator=args.UTRAnnotator, appris=args.appris, callback=args.callback, canonical=args.canonical, ccds=args.ccds, dbNSFP=args.dbNSFP, dbscSNV=args.dbscSNV, distance=args.distance, domains=args.domains, failed=args.failed, hgvs=args.hgvs, mane=args.mane, merged=args.merged, minimal=args.minimal,
                      mirna=args.mirna, mutfunc=args.mutfunc, numbers=args.numbers, protein=args.protein, refseq=args.refseq, shift_3prime=args.shift_3prime, shift_genomic=args.shift_genomic, transcript_id=args.transcript_id, transcript_version=args.transcript_version, tsl=args.tsl, uniprot=args.uniprot, variant_class=args.variant_class, vcf_string=args.vcf_string, xref_refseq=args.xref_refseq))
    case "vep_region_post":
        pprint.pprint(args.func(args.region, species=args.species, format=args.format, AncestralAllele=args.AncestralAllele, Blosum62=args.Blosum62, CADD=args.CADD, DisGeNET=args.DisGeNET, EVE=args.EVE, GO=args.GO, GeneSplicer=args.GeneSplicer, IntAct=args.IntAct, LoF=args.LoF, Mastermind=args.Mastermind, MaxEntScan=args.MaxEntScan,
                      NMD=args.NMD, Phenotypes=args.Phenotypes, SpliceAI=args.SpliceAI, UTRAnnotator=args.UTRAnnotator, appris=args.appris, callback=args.callback, canonical=args.canonical, ccds=args.ccds, dbNSFP=args.dbNSFP, dbscSNV=args.dbscSNV, distance=args.distance, domains=args.domains, failed=args.failed, hgvs=args.hgvs, mane=args.mane, merged=args.merged, minimal=args.minimal,
                      mirna=args.mirna, mutfunc=args.mutfunc, numbers=args.numbers, protein=args.protein, refseq=args.refseq, shift_3prime=args.shift_3prime, shift_genomic=args.shift_genomic, transcript_id=args.transcript_id, transcript_version=args.transcript_version, tsl=args.tsl, uniprot=args.uniprot, variant_class=args.variant_class, vcf_string=args.vcf_string, xref_refseq=args.xref_refseq))
    case "vep_hgvs_get":
        pprint.pprint(args.func(args.hgvs_notation, species=args.species, format=args.format, AncestralAllele=args.AncestralAllele, Blosum62=args.Blosum62, CADD=args.CADD, Conservation=args.Conservation, DisGeNET=args.DisGeNET, EVE=args.EVE, GO=args.GO, GeneSplicer=args.GeneSplicer, IntAct=args.IntAct, LoF=args.LoF, Mastermind=args.Mastermind, MaxEntScan=args.MaxEntScan,
                      NMD=args.NMD, Phenotypes=args.Phenotypes, SpliceAI=args.SpliceAI, UTRAnnotator=args.UTRAnnotator, ambiguous_hgvs=args.ambiguous_hgvs, appris=args.appris, callback=args.callback, canonical=args.canonical, ccds=args.ccds, dbNSFP=args.dbNSFP, dbscSNV=args.dbscSNV, distance=args.distance, domains=args.domains, failed=args.failed, hgvs=args.hgvs, mane=args.mane, merged=args.merged, minimal=args.minimal,
                      mirna=args.mirna, mutfunc=args.mutfunc, numbers=args.numbers, protein=args.protein, refseq=args.refseq, shift_3prime=args.shift_3prime, shift_genomic=args.shift_genomic, transcript_id=args.transcript_id, transcript_version=args.transcript_version, tsl=args.tsl, uniprot=args.uniprot, variant_class=args.variant_class, vcf_string=args.vcf_string, xref_refseq=args.xref_refseq))
    case "vep_id_get":
        pprint.pprint(args.func(args.id, species=args.species, format=args.format, AncestralAllele=args.AncestralAllele, Blosum62=args.Blosum62, CADD=args.CADD, Conservation=args.Conservation, DisGeNET=args.DisGeNET, EVE=args.EVE, GO=args.GO, GeneSplicer=args.GeneSplicer, IntAct=args.IntAct, LoF=args.LoF, Mastermind=args.Mastermind, MaxEntScan=args.MaxEntScan,
                      NMD=args.NMD, Phenotypes=args.Phenotypes, SpliceAI=args.SpliceAI, UTRAnnotator=args.UTRAnnotator, appris=args.appris, callback=args.callback, canonical=args.canonical, ccds=args.ccds, dbNSFP=args.dbNSFP, dbscSNV=args.dbscSNV, distance=args.distance, domains=args.domains, failed=args.failed, hgvs=args.hgvs, mane=args.mane, merged=args.merged, minimal=args.minimal,
                      mirna=args.mirna, mutfunc=args.mutfunc, numbers=args.numbers, protein=args.protein, refseq=args.refseq, shift_3prime=args.shift_3prime, shift_genomic=args.shift_genomic, transcript_id=args.transcript_id, transcript_version=args.transcript_version, tsl=args.tsl, uniprot=args.uniprot, variant_class=args.variant_class, vcf_string=args.vcf_string, xref_refseq=args.xref_refseq))
    case "vep_region_get":
        pprint.pprint(args.func(args.region, args.allele, species=args.species, format=args.format, AncestralAllele=args.AncestralAllele, Blosum62=args.Blosum62, CADD=args.CADD, Conservation=args.Conservation, DisGeNET=args.DisGeNET, EVE=args.EVE, GO=args.GO, GeneSplicer=args.GeneSplicer, IntAct=args.IntAct, LoF=args.LoF, Mastermind=args.Mastermind, MaxEntScan=args.MaxEntScan,
                      NMD=args.NMD, Phenotypes=args.Phenotypes, SpliceAI=args.SpliceAI, UTRAnnotator=args.UTRAnnotator, appris=args.appris, callback=args.callback, canonical=args.canonical, ccds=args.ccds, dbNSFP=args.dbNSFP, dbscSNV=args.dbscSNV, distance=args.distance, domains=args.domains, failed=args.failed, hgvs=args.hgvs, mane=args.mane, merged=args.merged, minimal=args.minimal,
                      mirna=args.mirna, mutfunc=args.mutfunc, numbers=args.numbers, protein=args.protein, refseq=args.refseq, shift_3prime=args.shift_3prime, shift_genomic=args.shift_genomic, transcript_id=args.transcript_id, transcript_version=args.transcript_version, tsl=args.tsl, uniprot=args.uniprot, variant_class=args.variant_class, vcf_string=args.vcf_string, xref_refseq=args.xref_refseq))
