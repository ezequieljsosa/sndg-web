# https://django-best-practices.readthedocs.io/en/latest/applications.html#managers

import datetime
from django.db import models
from django.db.models.query import QuerySet
from django.db.models import Prefetch


class ResourceQuerySet(QuerySet):

    def oai_compliant(self):
        return self.prefetch_related(
            "creators", "publishers"
        ).filter(deprecated=False, index_updated=False).exclude(creators=None, publishers=None)

    def publication_related(self,country):
        from .models import Publication
        qs = (Publication.objects.prefetch_related(
            "affiliations__organizations", "affiliations__author")
            .filter(affiliations__organizations__country=country))

        return (self.prefetch_related("ncbi_tax__names")
            .prefetch_related(Prefetch("targets__source",
                                       queryset=qs))
            .filter(deprecated=False, index_updated=False))


class BioResourceManager(models.Manager):

    def get_query_set(self):
        return ResourceQuerySet(self.model)

    def __getattr__(self, attr, *args):
        # see https://code.djangoproject.com/ticket/15062 for details
        if attr.startswith("_"):
            raise AttributeError
        return getattr(self.get_query_set(), attr, *args)
