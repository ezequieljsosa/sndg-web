# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.contrib.auth.models import User
from django.db import models

from bioresources.managers import BioResourceManager
from biosql.models import Taxon, Term
from . import compile_data
from django.urls import reverse

from django.utils.translation import gettext as _
from model_utils import Choices


def get_class(kls):
    parts = kls.split('.')
    module = ".".join(parts[:-1])
    m = __import__(module)
    for comp in parts[1:]:
        m = getattr(m, comp)
    return m


class ProcessStatus(models.Model):
    name = models.CharField(max_length=200)
    created_at = models.DateTimeField(auto_now_add=True)

    def __str__(self):
        return "ProcessStatus(%s)" % self.name

    def step(self, key):
        for step in self.steps.all():
            if step.name == key:
                return step
        raise IndexError("'%s' step not found" % key)

    def __repr__(self):
        return self.__str__()


class ProcessStatusStep(models.Model):
    name = models.CharField(max_length=200)
    class_name = models.CharField(max_length=200)
    class_identifier = models.CharField(max_length=200)
    cast_int_identifier = models.BooleanField(default=True)
    created_at = models.DateTimeField(auto_now_add=True)
    completed = models.BooleanField(default=False)
    process_status = models.ForeignKey(ProcessStatus, related_name="steps", on_delete=models.CASCADE)

    def __contains__(self, key):
        return len(self.units.filter(process_identifier=key))

    def append(self, db_identifier, process_identifier=None):
        if not process_identifier:
            process_identifier = db_identifier
        ProcessStatusStepProcessUnit(process_status_step=self,
                                     db_identifier=db_identifier, process_identifier=process_identifier).save()

    def results(self):
        clazz = get_class(self.class_name)
        for unit in self.units.all():
            if unit.db_identifier:
                identifier = int(unit.db_identifier) if self.cast_int_identifier else unit.db_identifier
                yield clazz.objects.get(**{self.class_identifier: identifier})

    def __str__(self):
        return "ProcessStatusStep(%s)" % self.name

    def __repr__(self):
        return self.__str__()


class ProcessStatusStepProcessUnit(models.Model):
    db_identifier = models.CharField(max_length=30, null=True)
    process_identifier = models.CharField(max_length=50)
    created_at = models.DateTimeField(auto_now_add=True)
    process_status_step = models.ForeignKey(ProcessStatusStep, related_name="units", on_delete=models.CASCADE)

    def __str__(self):
        return "PSStepPUnit('%s','%s')" % (self.db_identifier, self.process_identifier)

    def __repr__(self):
        return self.__str__()


class Person(models.Model):
    surname = models.CharField(max_length=200, blank=False)
    name = models.CharField(max_length=200, default="")
    scopus_id = models.CharField(max_length=200, null=True)
    scopus_names = models.TextField(null=True)
    email = models.EmailField()
    arg_affiliation = models.BooleanField(default=False)

    def rtype(self):
        return "person"

    def organizations(self):
        organizations = []
        for aff in self.affiliations.all():
            for org in aff.organizations.all():
                if org not in organizations:
                    organizations.append(org)
        return organizations

    def related_org_names(self):
        return [x.name for x in self.organizations()]

    def complete_name(self):
        return self.surname + ", " + self.name

    def type(self):
        return "person"

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

    def rtype(self):
        return "org"

    def __str__(self):
        return " ".join([x for x in [self.name, "|", self.country, self.city] if x])


class Resource(models.Model):
    RESOURCE_TYPES = Choices(
        [(i, x, _(x)) for i, x in enumerate([
            "PUBLICATION", "BIOPROJECT", "SEQUENCE", "ASSEMBLY", "GENOME", "READS",
            "STRUCTURE", "EXPRESSION", "BARCODE", "SAMPLE",
        ])]
    )

    facet_dict = {
        "assembly": ["species_name", "level", "assembly_type"],
        "gds": ["pdat", "gdstype"],
        "bioproject": ["sample_scope", "material"],  # , "capture_target", "method"
        "barcode": ["subdivision", "marker"],
    }

    id = models.AutoField(primary_key=True)

    type = models.CharField(max_length=10, choices=RESOURCE_TYPES)
    name = models.CharField(max_length=350, blank=False)
    description = models.TextField(blank=True)

    ncbi_tax = models.ForeignKey(Taxon, db_column="ncbi_tax", to_field="ncbi_taxon_id",
                                 on_delete=models.SET_NULL, null=True, related_name="bioresources")

    deprecated = models.BooleanField(default=False)
    index_updated = models.BooleanField(default=False)

    created_at = models.DateTimeField(auto_now_add=True)
    updated_at = models.DateTimeField(auto_now=True)

    objects = BioResourceManager()

    class Meta:
        unique_together = ('type', 'name',)

    def __str__(self):
        return self.name

    def get_absolute_url(self):
        return reverse('people.views.details', args=[str(self.id)])

    def compile(self):
        return compile_data

    def ncbi_tax_keywords(self):
        if self.ncbi_tax:
            return self.ncbi_tax.keywords.first().text
        return None

    def taxon_name(self):
        if self.ncbi_tax:
            return [x for x in self.ncbi_tax.names.all() if x.name_class == "scientific name"][0].name
        return None

    def oai_public(self):
        return True

    def oai_collections(self):
        return ["collection1"]

    def oai_communities(self):
        return ["community1"]

    def oai_submitter(self):
        return "sndg"

    def permalink(self):
        return "oai:" + str(self.id)

    def metadata_dc_language(self):
        return ["en"]

    def metadata_dc_rights(self):
        return ["info:eu-repo/semantics/openAccess"]

    def metadata_dc_format(self):
        return ["application/pdf"]

    def metadata_dc_publisher(self):
        return ["Faculdade de Medicina de São José do Rio Preto",
                "Programa de Pós-Graduação em Ciências da Saúde",
                "FAMERP",
                "BR",
                "Medicina Interna; Medicina e Ciências Correlatas"]


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

    SAMPLE_SCOPE_TYPES = Choices(
        [(i, x, _(x)) for i, x in enumerate([
            "monoisolate", 'multi-species', "environment", "synthetic", "other",
        ])]
    )

    MATERIAL_TYPES =    Choices([(i, x, _(x)) for i, x in enumerate([
            "genome", 'metagenome', "chromosome", "transcriptome", "reagent", "proteome",
        ])])

    CAPTURE_TYPES = Choices(
        [(i, x, _(x)) for i, x in enumerate([
            "whole", 'exome', "barcode", "TargetedLocusLoci",
        ])]
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
        for affiliation in self.targets.filter(source__type="pmc").all():
            ran += affiliation.source.author_names()
        return list(set(ran))

    def related_org_names(self):
        ran = []
        for affiliation in self.targets.filter(source__type="pmc").all():
            ran += affiliation.source.affiliation_names()
        return list(set(ran))


class Publication(Resource):
    doi = models.CharField(max_length=100)
    date_of_publication = models.DateField(max_length=200)
    pubmed_id = models.CharField(max_length=50, null=True)
    electronic_id = models.CharField(max_length=50, null=True)
    scopus_id = models.CharField(max_length=50, null=True)
    issn = models.CharField(max_length=50, null=True)

    def affiliation_names(self, all=False):
        affs = []
        for affiliation in self.affiliations.all():
            qs = affiliation.organizations
            if not all:
                qs = qs.filter(country="Argentina")
            for org in qs.all():
                if org.name not in affs:
                    affs.append(org.name)
        return affs

    def author_names(self):
        authors = []
        for affiliation in self.affiliations.filter(organizations__country="Argentina").all():
            if affiliation.author.complete_name() not in authors:
                authors.append(affiliation.author.complete_name())
        return authors

    def author_names_idx(self):
        return " ".join(self.author_names())

    def affiliation_names_idx(self):
        return " ".join(self.affiliation_names())


class Affiliation(models.Model):
    publication = models.ForeignKey(Publication, on_delete=models.CASCADE, related_name="affiliations")
    author = models.ForeignKey(Person, on_delete=models.CASCADE, related_name="affiliations")
    organizations = models.ManyToManyField(Organization)

    def __str__(self):
        return ("Affiliation: (%s) (%s) (%s) " % [str(x) for x in [self.author,self.publication] + self.organizations])


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

    ASSEMBLY_LEVEL = Choices(
        [(i, x, _(x)) for i, x in enumerate([
            "complete", 'chromosome', "scaffold", "contig",
        ])]
    )
    ASSEMBLY_TYPES = Choices(
        [(i, x, _(x)) for i, x in enumerate([
            "haploid", 'diploid', "other",
        ])]
    )

    intraspecific_name = models.CharField(max_length=250, null=True)
    species_name = models.CharField(max_length=200, null=True)
    level = models.CharField(max_length=50, null=True, choices=ASSEMBLY_LEVEL)
    ncbi_org = models.CharField(max_length=200, null=True)
    release_date = models.DateField(null=True)
    update_date = models.DateField(null=True)
    assembly_type = models.CharField(max_length=50, null=True, choices=ASSEMBLY_TYPES)

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


class Barcode(Resource):
    country = models.CharField(max_length=100)
    subdivision = models.CharField(max_length=150)
    marker = models.CharField(max_length=50, null=True)
    image_url = models.URLField(null=True)
    bold_org = models.CharField(max_length=255, null=True)
    collectors = models.CharField(max_length=255, null=True)


class ResourceProperty(models.Model):
    organization = models.ForeignKey(Organization, on_delete=models.DO_NOTHING)
    resource = models.ForeignKey(Resource, on_delete=models.CASCADE)
    term = models.ForeignKey(Term, on_delete=models.DO_NOTHING)


class ResourcePropertyValue(models.Model):
    property = models.ForeignKey(ResourceProperty, on_delete=models.DO_NOTHING)
    value = models.CharField(max_length=200)


class Sample(Resource):
    country = models.CharField(max_length=100)
    subdivision = models.CharField(max_length=150)
    collection_date = models.DateField(null=True)
    publication_date = models.DateField(null=True)
    update_date = models.DateField(null=True)


class ReadsArchive(Resource):
    release_date = models.DateField(null=True)
    update_date = models.DateField(null=True)
