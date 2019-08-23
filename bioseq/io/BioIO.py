import os
from tqdm import tqdm
from datetime import datetime
from typing import Iterable

from Bio.SeqRecord import SeqRecord

from bioseq.models.Taxon import Taxon, TaxonName, TaxIdx
from bioseq.models.Ontology import Ontology
from bioseq.models.Term import Term, TermRelationship, TermDbxref, TermIdx
from bioseq.models.Biodatabase import Biodatabase
from bioseq.models.Bioentry import Bioentry, BioentryQualifierValue, BioentryDbxref
from bioseq.models.Biosequence import Biosequence
from bioseq.models.Seqfeature import Seqfeature, SeqfeatureDbxref, SeqfeatureQualifierValue
from bioseq.models.Location import Location
from bioseq.models.Dbxref import Dbxref, DbxrefQualifierValue

from bioresources.models.Resource import Resource
from bioresources.models.Person import Person
from bioresources.models.Organization import Organization
from bioresources.models.Affiliation import Affiliation
from bioresources.models.Identity import Identity
from bioresources.models.ExternalId import ExternalId
from bioresources.models.Publication import Publication
from bioresources.models.Tool import Tool
from bioresources.models.Assembly import Assembly
from bioresources.models.ReadsArchive import ReadsArchive
from bioresources.models.BioProject import BioProject
from bioresources.models.Sample import Sample
from bioresources.models.ResourceRelation import ResourceRelation
from bioresources.models.ResourceProperty import ResourceProperty, ResourcePropertyValue


class BioIO:
    def __init__(self, biodb_name: str,ncbi_tax:int):
        self.biodb_name = biodb_name
        self.ncbi_tax = ncbi_tax

    def create_db(self):
        self.genomedb = Biodatabase(self.biodb_name)
        self.genomedb.save()
        self.proteindb = Biodatabase(self.biodb_name  + "_prots")
        self.genomedb.save()
        self.rnadb = Biodatabase(self.biodb_name + "_rnas")
        self.genomedb.save()

    def process_record_list(self, seq_record_iterator: Iterable[SeqRecord],contig_count:int):
        for seqrecord in tqdm(seq_record_iterator,total=contig_count):
            self.process_seqrecord(seqrecord)
