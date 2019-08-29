# -*- coding: utf-8 -*-
from __future__ import unicode_literals
from model_utils import Choices
from django.utils.translation import gettext_lazy as _, ngettext as __
from django.db import models
from django.urls import reverse

from .Person import Person
from .Organization import Organization
from .RKeyword import RKeyword
from ..managers.BioResourceManager import BioResourceManager

from bioseq.models.Taxon import Taxon
from polymorphic.models import PolymorphicModel


def get_class(kls):
    parts = kls.split('.')
    module = ".".join(parts[:-1])
    m = _import_(module)
    for comp in parts[1:]:
        m = getattr(m, comp)
    return m


class Resource(PolymorphicModel):
    RESOURCE_TYPES = Choices(
        *([(i, x, _(x)) for i, x in enumerate([
            "PUBLICATION", "BIOPROJECT", "SEQUENCE", "ASSEMBLY", "GENOME", "READS",
            "STRUCTURE", "EXPRESSION", "BARCODE", "SAMPLE", "TOOL",
        ])] + [(50, "UNPROCESSED", _("UNPROCESSED")), (40, "PROTEIN", _("PROTEIN")),
               (30, "ORGANIZATION", _("ORGANIZATION")), (20, "PERSON", _("PERSON"))])
    )

    name2code = {
        n: idx for idx, n in RESOURCE_TYPES
    }

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

    creators = models.ManyToManyField(Organization, related_name="created_resources", blank=True)
    publishers = models.ManyToManyField(Organization, related_name="published_resources", blank=True)
    keywords = models.ManyToManyField(RKeyword, related_name="associated_resources")

    ncbi_tax = models.ForeignKey(Taxon, db_column="ncbi_tax", to_field="ncbi_taxon_id",blank=True,
                                 on_delete=models.SET_NULL, null=True, related_name="bioresources")

    deprecated = models.BooleanField(default=False)
    index_updated = models.BooleanField(default=False)

    created_at = models.DateTimeField(auto_now_add=True)
    updated_at = models.DateTimeField(auto_now=True)

    objects = BioResourceManager()

    class Meta:
        unique_together = ('type', 'name',)

    def _str_(self):
        return self.name

    def save(self, *args, **kwargs):
        self.type = self.__class__.TYPE
        super(Resource, self).save(*args, **kwargs)

    def get_absolute_url(self):
        return reverse('bioresources:%s_view' % Resource.RESOURCE_TYPES[self.type].lower(), args=[str(self.id)])

    def type_name(self):
        return Resource.RESOURCE_TYPES[self.type].lower()

    def compile(self):
        templ = get_template("resources/xoai_resource.xml")
        return templ.render({"r": self})

    def ncbi_tax_keywords(self):
        return self.ncbi_tax.keywords.text if self.ncbi_tax else None

    def taxon_name(self):
        return self.ncbi_tax.scientific_name() if self.ncbi_tax else None

    def oai_public(self):
        return True

    def oai_collections(self):
        return ["sndg." + Resource.RESOURCE_TYPES[self.type].lower()]

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


class Collaboration(models.Model):
    COLLABORATION_TYPES = Choices(
        (1, "owner", _("owner")),
        (2, "only_producer", _("only_producer")),
        (3, "only_use", _("only_use")),
        (4, "other", _("other")),
    )
    rev_types = {x[1]: x[0] for x in COLLABORATION_TYPES}

    resource = models.ForeignKey(Resource, related_name="collaborations", on_delete=models.PROTECT)
    person = models.ForeignKey(Person, related_name="collaborations", on_delete=models.PROTECT)
    type = models.PositiveIntegerField(choices=COLLABORATION_TYPES)
    info = models.TextField(null=True, blank=True)


class UnprocessedResource(Resource):
    TYPE = Resource.RESOURCE_TYPES.UNPROCESSED
    future_type = models.PositiveIntegerField(choices=Resource.RESOURCE_TYPES)
