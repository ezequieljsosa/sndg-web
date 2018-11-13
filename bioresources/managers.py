# https://django-best-practices.readthedocs.io/en/latest/applications.html#managers

import datetime
from django.db import models
from django.db.models.query import QuerySet

class ResourceQuerySet(QuerySet):
    def live(self):
        """Filter out posts that aren't ready to be published"""
        now = datetime.datetime.now()
        return self.filter(date_published__lte=now, status="published")

class BioResourceManager(models.Manager):
    def get_query_set(self):
        return ResourceQuerySet(self.model)
    def __getattr__(self, attr, *args):
        # see https://code.djangoproject.com/ticket/15062 for details
        if attr.startswith("_"):
            raise AttributeError
        return getattr(self.get_query_set(), attr, *args)