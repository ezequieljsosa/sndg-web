# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.contrib.auth.models import User
from django.db import models

from bioresources.managers import BioResourceManager
from biosql.models import Taxon, Term, Ontology
from django.db.models import Q
from . import compile_data
from django.urls import reverse

from django.utils.translation import gettext as _
from django.utils.translation import gettext_lazy as __
from model_utils import Choices
from django.template.loader import get_template

from django.conf import settings


def get_class(kls):
    parts = kls.split('.')
    module = ".".join(parts[:-1])
    m = __import__(module)
    for comp in parts[1:]:
        m = getattr(m, comp)
    return m


from django.contrib.auth.models import AbstractUser


# class User(AbstractUser):
#     pass

class ProcessStatus(models.Model):
    name = models.CharField(max_length=200, help_text=__('Process name'))
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
    # TODO add manager to query affiliations+organizations
    surname = models.CharField(max_length=200, blank=False)
    name = models.CharField(max_length=200, default="")
    scopus_id = models.CharField(max_length=200, null=True)
    scopus_names = models.TextField(null=True)
    email = models.EmailField()

    def __str__(self):
        return self.name + " " + self.surname

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


class Identity(models.Model):
    person = models.ForeignKey(Person, on_delete=models.CASCADE)
    identifier = models.CharField(max_length=200)
    email = models.EmailField(null=True)
    url = models.URLField(null=True)
    authority = models.CharField(max_length=200)
    created_at = models.DateTimeField(auto_now_add=True)
    updated_at = models.DateTimeField()
    ends = models.DateTimeField()

    def __str__(self):
        return self.identifier


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


class RKeyword(models.Model):
    name = models.CharField(max_length=200, unique=True)

    def save(self, *args, **kwargs):
        self.name = self.name.lower()
        super(RKeyword, self).save(*args, **kwargs)

    def __str__(self):
        return self.name


class Resource(models.Model):
    RESOURCE_TYPES = Choices(
        *[(i, x, _(x)) for i, x in enumerate([
            "PUBLICATION", "BIOPROJECT", "SEQUENCE", "ASSEMBLY", "GENOME", "READS",
            "STRUCTURE", "EXPRESSION", "BARCODE", "SAMPLE", "TOOL", "PROTEIN",
        ])]
    )

    facet_dict = {
        "assembly": ["species_name", "level", "assembly_type"],
        "gds": ["pdat", "gdstype"],
        "bioproject": ["sample_scope", "material"],  # , "capture_target", "method"
        "barcode": ["subdivision", "marker"],
    }

    id = models.AutoField(primary_key=True)

    type = models.PositiveIntegerField(choices=RESOURCE_TYPES)
    name = models.CharField(max_length=350, blank=False)
    description = models.TextField(blank=True)

    creators = models.ManyToManyField(Organization, related_name="created_resources",blank=True)
    publishers = models.ManyToManyField(Organization, related_name="published_resources",blank=True)
    keywords = models.ManyToManyField(RKeyword, related_name="associated_resources")

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

    def save(self, *args, **kwargs):
        self.type = self.__class__.TYPE
        super(Resource, self).save(*args, **kwargs)

    def get_absolute_url(self):
        return reverse('bioresources:%s_view' % Resource.RESOURCE_TYPES[self.type].lower() , args=[str(self.id)])

    def compile(self):
        templ = get_template("resources/xoai_resource.xml")
        return templ.render({"r": self})

        #return compile_data

    def ncbi_tax_keywords(self):
        return self.ncbi_tax.keywords.text if self.ncbi_tax else None

    def taxon_name(self):
        return self.ncbi_tax.scientific_name() if self.ncbi_tax else None

    def oai_public(self):
        return True

    def oai_collections(self):
        return [ "sndg." +  Resource.RESOURCE_TYPES[self.type].lower() ]

    def oai_communities(self):
        return ["sndg"]

    def oai_submitter(self):
        return "sndg"

    def handle(self):
        return Resource.RESOURCE_TYPES[self.type].lower() + "/" + str(self.id)

    def permalink(self):
        return "oai:" + settings.OAIPMH_DOMAIN + ":" + self.handle()

    def metadata_dc_language(self):
        return ["en"]

    def metadata_dc_rights(self):
        return ["info:eu-repo/semantics/openAccess"]

    def metadata_dc_format(self):
        return []

    def metadata_dc_creator(self):
        return [x.name for x in self.creators.all()]

    def metadata_dc_publisher(self):
        return [x.name for x in self.publishers.all()]


class ResourceRelation(models.Model):
    source = models.ForeignKey(Resource, on_delete=models.CASCADE, related_name="targets")
    target = models.ForeignKey(Resource, on_delete=models.CASCADE, related_name="sources")
    role = models.CharField(max_length=200, blank=False)
    deprecated = models.BooleanField(default=False)

    class Meta:
        unique_together = (('source', 'target', 'deprecated'),)
        verbose_name_plural = __("Resource Relations")

    def __str__(self):
        return ("(" + self.source.type + ":" + str(self.source.id) +
                ") -> " + "(" + self.target.type + ":" + str(self.target.id) + ")")


class BioProject(Resource):
    """
    https://www.ncbi.nlm.nih.gov/books/NBK169438/
    """

    SAMPLE_SCOPE_TYPES = Choices(
        *[(i, x, _(x)) for i, x in enumerate([
            "monoisolate", 'multi-species', "environment", "synthetic", "other",
        ])]
    )

    MATERIAL_TYPES = Choices(
        *[(i, x, _(x)) for i, x in enumerate([
            "genome", 'metagenome', "chromosome", "transcriptome", "reagent", "proteome",
        ])])

    CAPTURE_TYPES = Choices(
        *[(i, x, _(x)) for i, x in enumerate([
            "whole", 'exome', "barcode", "TargetedLocusLoci",
        ])]
    )

    TYPE = Resource.RESOURCE_TYPES.BIOPROJECT

    sample_scope = models.PositiveIntegerField(choices=SAMPLE_SCOPE_TYPES, null=True)
    material = models.PositiveIntegerField(choices=MATERIAL_TYPES, null=True)
    capture = models.PositiveIntegerField(choices=CAPTURE_TYPES, null=True)

    target = models.CharField(max_length=200, null=True)
    submitters = models.TextField(null=True)

    method = models.CharField(max_length=200, null=True)

    # objetive = models.CharField(max_length=200, blank=False)

    class Meta:
        verbose_name_plural = __("BioProjects")

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
    TYPE = Resource.RESOURCE_TYPES.PUBLICATION

    doi = models.CharField(max_length=100)
    date_of_publication = models.DateField(max_length=200)
    pubmed_id = models.CharField(max_length=50, null=True)
    electronic_id = models.CharField(max_length=50, null=True)
    scopus_id = models.CharField(max_length=50, null=True)
    issn = models.CharField(max_length=50, null=True)

    class Meta:
        verbose_name_plural = __("Publications")

    def affiliation_names(self, country=False):
        affs = []
        for affiliation in self.affiliations.all():
            qs = affiliation.organizations
            if country:
                qs = qs.filter(country=country)
            for org in qs.all():
                if org.name not in affs:
                    affs.append(org.name)
        return affs

    def author_names(self, country=None):
        authors = []
        if country:
            sqs = self.affiliations.filter(organizations__country=country)
        else:
            sqs = self.affiliations

        for affiliation in sqs.all():
            if affiliation.author.complete_name() not in authors:
                authors.append(affiliation.author.complete_name())
        return authors

    def author_names_idx(self, country=None):
        return " ".join(self.author_names(country))

    def affiliation_names_idx(self, country=None):
        return " ".join(self.affiliation_names(country))


class Affiliation(models.Model):
    publication = models.ForeignKey(Publication, on_delete=models.CASCADE, related_name="affiliations")
    author = models.ForeignKey(Person, on_delete=models.CASCADE, related_name="affiliations")
    organizations = models.ManyToManyField(Organization)

    class Meta:
        verbose_name_plural = __("Affiliations")

    def __str__(self):
        return ("Affiliation: (%s) (%s) (%s) " % [str(x) for x in [self.author, self.publication] + self.organizations])


class ExternalId(models.Model):
    organization = models.ForeignKey(Organization, on_delete=models.CASCADE)
    identifier = models.CharField(max_length=20)
    type = models.CharField(max_length=20)
    resource = models.ForeignKey(Resource, on_delete=models.CASCADE, related_name="external_ids")

    def __str__(self):
        return self.organization.name + ":" + self.identifier


class Structure(Resource):
    TYPE = Resource.RESOURCE_TYPES.STRUCTURE

    pdbClass = models.CharField(max_length=50, null=True)
    deposit_date = models.DateField(null=True)
    method = models.CharField(max_length=50, null=True)
    org_list = models.TextField(null=True)

    class Meta:
        verbose_name_plural = __("Structures")

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
        *[(i, x, _(x)) for i, x in enumerate([
            "complete genome", 'chromosome', "scaffold", "contig",
        ])]
    )
    ASSEMBLY_TYPES = Choices(
        *[(i, x, _(x)) for i, x in enumerate([
            "haploid", 'diploid', "other",
        ])]
    )

    TYPE = Resource.RESOURCE_TYPES.ASSEMBLY

    intraspecific_name = models.CharField(max_length=250, null=True)
    species_name = models.CharField(max_length=200, null=True)
    level = models.PositiveIntegerField(null=True, choices=ASSEMBLY_LEVEL)

    ncbi_org = models.CharField(max_length=200, null=True)
    release_date = models.DateField(null=True)
    update_date = models.DateField(null=True)
    assembly_type = models.PositiveIntegerField(null=True, choices=ASSEMBLY_TYPES)

    class Meta:
        verbose_name_plural = __("Assemblies")

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
    TYPE = Resource.RESOURCE_TYPES.EXPRESSION

    pdat = models.DateField(null=True)
    gdstype = models.CharField(max_length=250, null=True)
    submitters = models.TextField(null=True)

    class Meta:
        verbose_name_plural = __("Expression")

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
    TYPE = Resource.RESOURCE_TYPES.BARCODE

    country = models.CharField(max_length=100)
    subdivision = models.CharField(max_length=150)
    marker = models.CharField(max_length=50, null=True)
    image_url = models.URLField(null=True)
    bold_org = models.CharField(max_length=255, null=True)
    collectors = models.CharField(max_length=255, null=True)

    class Meta:
        verbose_name_plural = __("Barcodes")


class ResourceProperty(models.Model):
    organization = models.ForeignKey(Organization, on_delete=models.DO_NOTHING)
    resource = models.ForeignKey(Resource, on_delete=models.CASCADE, related_name="properties")
    term = models.ForeignKey(Term, on_delete=models.DO_NOTHING)


class ResourcePropertyValue(models.Model):
    property = models.OneToOneField(ResourceProperty, on_delete=models.CASCADE, related_name="value")
    value = models.CharField(max_length=200)


class Sample(Resource):
    TYPE = Resource.RESOURCE_TYPES.SAMPLE

    origin_props = ['Origin (developed or donated from)', 'geo_loc_name', 'geographic location (country and/or sea)',
                    'birth_location', 'country_of_birth', 'geo-loc-name']

    country = models.CharField(max_length=100)
    subdivision = models.CharField(max_length=150)
    collection_date = models.DateField(null=True)
    publication_date = models.DateField(null=True)
    update_date = models.DateField(null=True)

    class Meta:
        verbose_name_plural = __("Samples")

    def origin_dict(self):
        props = self.properties.prefetch_related("value").filter(
            Q(term__ontology__name="NCBI sample") & Q(term__identifier__in=Sample.origin_props))
        return {prop.term.name: prop.value.value for prop in props}


class ReadsArchive(Resource):
    TYPE = Resource.RESOURCE_TYPES.READS

    release_date = models.DateField(null=True)
    update_date = models.DateField(null=True)

    class Meta:
        verbose_name_plural = __("Reads Archive")

class Tool(Resource):
    TYPE = Resource.RESOURCE_TYPES.TOOL

    TYPES = Choices(
        *[(i, x, _(x)) for i, x in enumerate([
            'app', 'database', 'library', 'plugin', 'program', 'webserver',
        ])]
    )

    tool_type = models.PositiveSmallIntegerField(choices=TYPES)
    url = models.URLField(null=True)

    class Meta:
        verbose_name_plural = __("Tools")