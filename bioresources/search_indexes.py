import datetime
from haystack import indexes
from django.db.models import Prefetch
from .models import Publication, Assembly, BioProject, Expression, Structure

"""
python manage.py build_solr_schema > /opt/solr-7.3.1/server/solr/eze_core/schema.xml
python manage.py rebuild_index
python manage.py update_index
"""


class PublicationIndex(indexes.SearchIndex, indexes.Indexable):
    text = indexes.CharField(document=True, use_template=True)

    name = indexes.CharField(model_attr='name')
    description = indexes.CharField(model_attr='description')

    pub_date = indexes.DateTimeField(model_attr='date_of_publication', null=True)

    type = indexes.CharField(model_attr='type', faceted=True)

    authors = indexes.MultiValueField(model_attr='author_names', faceted=True)
    affiliations = indexes.MultiValueField(model_attr='affiliation_names', faceted=True)

    doi = indexes.CharField(model_attr='doi', null=True)
    electronic_id = indexes.CharField(model_attr='electronic_id', null=True)
    scopus_id = indexes.CharField(model_attr='scopus_id', null=True)
    issn = indexes.CharField(model_attr='issn', null=True)
    pubmed_id = indexes.CharField(model_attr='pubmed_id', null=True)

    def get_model(self):
        return Publication

    def index_queryset(self, using=None):
        """Used when the entire index for model is updated."""
        return self.get_model().objects.prefetch_related(
            "affiliations", "affiliations__organizations", "affiliations__author"
        ).filter(deprecated=False, index_updated=False)


class StructureIndex(indexes.SearchIndex, indexes.Indexable):
    text = indexes.CharField(document=True, use_template=True)
    name = indexes.CharField(model_attr='name')
    description = indexes.CharField(model_attr='description')
    type = indexes.CharField(model_attr='type', faceted=True)

    pdbClass = indexes.CharField(model_attr='pdbClass', faceted=True, null=True)
    deposit_date = indexes.DateTimeField(model_attr='deposit_date', null=True)
    method = indexes.CharField(model_attr='method', faceted=True, null=True)
    org_list = indexes.CharField(model_attr='org_list', null=True)

    authors = indexes.MultiValueField(model_attr='related_author_names', faceted=True)
    affiliations = indexes.MultiValueField(model_attr='related_org_names', faceted=True)
    taxon = indexes.MultiValueField(model_attr='taxon_name', faceted=True)

    def get_model(self):
        return Structure

    def index_queryset(self, using=None):
        """Used when the entire index for model is updated."""
        qs = (Publication.objects.prefetch_related(
            "affiliations__organizations", "affiliations__author")
            .filter(affiliations__organizations__country="Argentina"))
        return (self.get_model().objects.prefetch_related("ncbi_tax__names")
            .prefetch_related(Prefetch("targets__source",
                                       queryset=qs))
            .filter(deprecated=False, index_updated=False))


class AssemblyIndex(indexes.SearchIndex, indexes.Indexable):
    text = indexes.CharField(document=True, use_template=True)
    name = indexes.CharField(model_attr='name')
    description = indexes.CharField(model_attr='description')
    type = indexes.CharField(model_attr='type', faceted=True)

    intraspecific_name = indexes.CharField(model_attr='intraspecific_name', null=True)
    species_name = indexes.CharField(model_attr='species_name', faceted=True, null=True)
    level = indexes.CharField(model_attr='level', faceted=True, null=True)
    ncbi_org = indexes.CharField(model_attr='ncbi_org', null=True)
    release_date = indexes.DateTimeField(model_attr='release_date', null=True)
    update_date = indexes.DateTimeField(model_attr='update_date', null=True)
    assembly_type = indexes.CharField(model_attr='assembly_type', faceted=True, null=True)

    authors = indexes.MultiValueField(model_attr='related_author_names', faceted=True)
    affiliations = indexes.MultiValueField(model_attr='related_org_names', faceted=True)
    taxon = indexes.MultiValueField(model_attr='taxon_name', faceted=True)

    def get_model(self):
        return Assembly

    def index_queryset(self, using=None):
        qs = (Publication.objects.prefetch_related(
            "affiliations__organizations", "affiliations__author")
            .filter(affiliations__organizations__country="Argentina"))

        return (self.get_model().objects.prefetch_related("ncbi_tax__names")
            .prefetch_related(Prefetch("targets__source",
                                       queryset=qs))
            .filter(deprecated=False, index_updated=False))


class ExpressionIndex(indexes.SearchIndex, indexes.Indexable):
    text = indexes.CharField(document=True, use_template=True)
    name = indexes.CharField(model_attr='name')
    description = indexes.CharField(model_attr='description')
    type = indexes.CharField(model_attr='type', faceted=True)

    pdat = indexes.CharField(model_attr='pdat', null=True)
    gdstype = indexes.CharField(model_attr='gdstype', faceted=True, null=True)
    submitters = indexes.CharField(model_attr='submitters', null=True)

    authors = indexes.MultiValueField(model_attr='related_author_names', faceted=True)
    affiliations = indexes.MultiValueField(model_attr='related_org_names', faceted=True)
    taxon = indexes.MultiValueField(model_attr='taxon_name', faceted=True)

    def get_model(self):
        return Expression

    def index_queryset(self, using=None):

        qs = (Publication.objects.prefetch_related(
            "affiliations__organizations", "affiliations__author")
            .filter(affiliations__organizations__country="Argentina"))
        return (self.get_model().objects.prefetch_related("ncbi_tax__names")
            .prefetch_related(Prefetch("targets__source", queryset=qs))
            .filter(deprecated=False, index_updated=False))


class BioProjectIndex(indexes.SearchIndex, indexes.Indexable):
    text = indexes.CharField(document=True, use_template=True)
    name = indexes.CharField(model_attr='name')
    description = indexes.CharField(model_attr='description')
    type = indexes.CharField(model_attr='type', faceted=True)

    sample_scope = indexes.CharField(model_attr='sample_scope', faceted=True, null=True)
    material = indexes.CharField(model_attr='material', faceted=True, null=True)
    capture = indexes.CharField(model_attr='capture', faceted=True, null=True)

    target = indexes.CharField(model_attr='target', faceted=True, null=True)
    submitters = indexes.CharField(model_attr='submitters', null=True)

    method = indexes.CharField(model_attr='method', faceted=True, null=True)

    authors = indexes.MultiValueField(model_attr='related_author_names', faceted=True)
    affiliations = indexes.MultiValueField(model_attr='related_org_names', faceted=True)
    taxon = indexes.MultiValueField(model_attr='taxon_name', faceted=True)

    def get_model(self):
        return BioProject

    def index_queryset(self, using=None):

        qs = (Publication.objects.prefetch_related(
            "affiliations__organizations", "affiliations__author")
            .filter(affiliations__organizations__country="Argentina"))

        return (self.get_model().objects.prefetch_related("ncbi_tax__names")
            .prefetch_related(Prefetch("targets__source",
                                       queryset=qs))
            .filter(deprecated=False, index_updated=False))
