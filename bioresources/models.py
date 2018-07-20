# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.contrib.auth.models import User
from django.db import models
from biosql.models import Taxon


class Person(models.Model):
    surname = models.CharField(max_length=200, blank=False)
    name = models.CharField(max_length=200, default="")
    scopus_id = models.CharField(max_length=200, null=True)
    scopus_names = models.TextField(null=True)
    email = models.EmailField()
    arg_affiliation = models.BooleanField(default=False)

    def complete_name(self):
        return self.surname + ", " + self.name

    def __str__(self):
        return self.name + " " + self.surname


class Identity(models.Model):
    person = models.ForeignKey(Person, on_delete=models.CASCADE)
    identifier = models.CharField(max_length=200)
    email = models.EmailField(null=True)
    url = models.URLField(null=True)
    authority = models.CharField(max_length=200)
    created_at = models.DateTimeField(auto_now_add=True)
    updated_at = models.DateTimeField()
    ends = models.DateTimeField()


class Organization(models.Model):
    name = models.CharField(max_length=200)
    url = models.URLField(null=True)
    country = models.CharField(max_length=200, null=True)
    city = models.CharField(max_length=200, null=True)
    scopus_id = models.CharField(max_length=200, null=True)
    scopus_names = models.TextField(null=True)

    def __str__(self):
        return " ".join([x for x in [self.name, "|", self.country, self.city] if x])


class Resource(models.Model):
    PUBLICATION = 'pmc'
    BIOPROJECT = 'bioproject'
    SEQUENCE = 'nuccore'
    ASSEMBLY = 'assembly'
    GENOME = 'genome'
    READS = 'sra'
    STRUCTURE = 'structure'
    EXPRESSION = 'expression'

    RESOURCE_TYPES = (
        (PUBLICATION, PUBLICATION),
        (BIOPROJECT, BIOPROJECT),
        (SEQUENCE, SEQUENCE),
        (ASSEMBLY, ASSEMBLY),
        (GENOME, GENOME),
        (READS, READS),
        (STRUCTURE, STRUCTURE),
    )

    type = models.CharField(max_length=10, choices=RESOURCE_TYPES)
    name = models.CharField(max_length=350, blank=False)
    description = models.TextField(blank=True)

    ncbi_tax = models.ForeignKey(Taxon, db_column="ncbi_tax", to_field="ncbi_taxon_id",
                                 on_delete=models.SET_NULL, null=True, related_name="bioresources")

    deprecated = models.BooleanField(default=False)
    index_updated = models.BooleanField(default=False)

    created_at = models.DateTimeField(auto_now_add=True)
    updated_at = models.DateTimeField(auto_now=True)

    class Meta:
        unique_together = ('type', 'name',)

    def taxon_name(self):
        if self.ncbi_tax:
            return [x for x in self.ncbi_tax.names.all() if x.name_class == "scientific name"][0].name
        return None


    def __str__(self):
        return self.name


class ResourceRelation(models.Model):
    source = models.ForeignKey(Resource, on_delete=models.CASCADE, related_name="targets")
    target = models.ForeignKey(Resource, on_delete=models.CASCADE, related_name="sources")
    role = models.CharField(max_length=200, blank=False)
    deprecated = models.BooleanField(default=False)

    class Meta:
        unique_together = (('source', 'target', 'deprecated'),)

    def __str__(self):
        return ("(" + self.source.type + ":" + str(self.source.id) +
                ") -> " + "(" + self.target.type + ":" + str(self.target.id) + ")")


class BioProject(Resource):
    """
    https://www.ncbi.nlm.nih.gov/books/NBK169438/
    """
    ss_monoisolate = 'monoisolate'
    ss_multispecies = 'multi-species'
    ss_environment = 'environment'
    ss_synthetic = 'synthetic'
    ss_other = 'other'

    SAMPLE_SCOPE_TYPES = (
        (ss_monoisolate, ss_monoisolate),
        (ss_multispecies, ss_multispecies),
        (ss_environment, ss_environment),
        (ss_synthetic, ss_synthetic),
        (ss_other, ss_other),

    )

    m_genome = 'genome'
    m_metagenome = 'metagenome'
    m_chromosome = 'chromosome'
    m_transcriptome = 'transcriptome'
    m_reagent = 'reagent'
    m_proteome = 'proteome'

    MATERIAL_TYPES = (
        (m_genome, m_genome),
        (m_metagenome, m_metagenome),
        (m_chromosome, m_chromosome),
        (m_transcriptome, m_transcriptome),
        (m_reagent, m_reagent),
        (m_proteome, m_proteome),

    )

    c_whole = "whole"
    c_exome = "exome"
    c_barcode = "barcode"
    c_TargetedLocusLoci = "TargetedLocusLoci"

    CAPTURE_TYPES = (
        (c_whole, c_whole),
        (c_exome, c_exome),
        (c_barcode, c_barcode),
        (c_TargetedLocusLoci, c_TargetedLocusLoci),
    )

    sample_scope = models.CharField(max_length=20, choices=SAMPLE_SCOPE_TYPES, null=True)
    material = models.CharField(max_length=20, choices=MATERIAL_TYPES, null=True)
    capture = models.CharField(max_length=200, choices=CAPTURE_TYPES, null=True)

    target = models.CharField(max_length=200, null=True)
    submitters = models.TextField(null=True)

    method = models.CharField(max_length=200, null=True)

    # objetive = models.CharField(max_length=200, blank=False)

    def related_author_names(self):
        ran = []
        for affiliation in self.targets.all():
            ran += affiliation.source.author_names()
        return list(set(ran))

    def related_org_names(self):
        ran = []
        for affiliation in self.targets.all():
            ran += affiliation.source.affiliation_names()
        return list(set(ran))


class Publication(Resource):
    doi = models.CharField(max_length=100)
    date_of_publication = models.DateField(max_length=200)
    pubmed_id = models.CharField(max_length=50, null=True)
    electronic_id = models.CharField(max_length=50, null=True)
    scopus_id = models.CharField(max_length=50, null=True)
    issn = models.CharField(max_length=50, null=True)

    def affiliation_names(self):
        affs = []
        for affiliation in self.affiliations.all():
            for org in affiliation.organizations.all():
                if org.name not in affs:
                    affs.append(org.name)
        return affs

    def author_names(self):
        authors = []
        for affiliation in self.affiliations.all():
            if affiliation.author.complete_name() not in authors:
                authors.append(affiliation.author.complete_name())
        return authors

    def author_names_idx(self):
        return " ".join(self.author_names())

    def affiliation_names_idx(self):
        return " ".join(self.affiliation_names())


class Affiliation(models.Model):
    publication = models.ForeignKey(Publication, on_delete=models.CASCADE, related_name="affiliations")
    author = models.ForeignKey(Person, on_delete=models.CASCADE, related_name="authors")
    organizations = models.ManyToManyField(Organization)


class ExternalId(models.Model):
    organization = models.ForeignKey(Organization, on_delete=models.CASCADE)
    identifier = models.CharField(max_length=20)
    type = models.CharField(max_length=20)
    resource = models.ForeignKey(Resource, on_delete=models.CASCADE, related_name="external_ids")

    def __str__(self):
        return self.organization.name + ":" + self.identifier


class Structure(Resource):
    pdbClass = models.CharField(max_length=50, null=True)
    deposit_date = models.DateField(null=True)
    method = models.CharField(max_length=50, null=True)
    org_list = models.TextField(null=True)

    def __str__(self):
        return self.name

    def related_author_names(self):
        ran = []
        for affiliation in self.targets.all():
            ran += affiliation.source.author_names()
        return list(set(ran))

    def related_org_names(self):
        ran = []
        for affiliation in self.targets.all():
            ran += affiliation.source.affiliation_names()
        return list(set(ran))


class Assembly(Resource):
    intraspecific_name = models.CharField(max_length=250, null=True)
    species_name = models.CharField(max_length=200, null=True)
    level = models.CharField(max_length=50, null=True)
    ncbi_org = models.CharField(max_length=200, null=True)
    release_date = models.DateField(null=True)
    update_date = models.DateField(null=True)
    assembly_type = models.CharField(max_length=50, null=True, choices=(
        ("haploid", "haploid"), ("diploid", "diploid"), ("other", "other")))

    def related_author_names(self):
        ran = []
        for affiliation in self.targets.all():
            ran += affiliation.source.author_names()
        return list(set(ran))

    def related_org_names(self):
        ran = []
        for affiliation in self.targets.all():
            ran += affiliation.source.affiliation_names()
        return list(set(ran))


class Expression(Resource):
    pdat = models.DateField(null=True)
    gdstype = models.CharField(max_length=250, null=True)
    submitters = models.TextField(null=True)

    def related_author_names(self):
        ran = []
        for affiliation in self.targets.all():
            ran += affiliation.source.author_names()
        return list(set(ran))

    def related_org_names(self):
        ran = []
        for affiliation in self.targets.all():
            ran += affiliation.source.affiliation_names()
        return list(set(ran))
