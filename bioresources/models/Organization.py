# -*- coding: utf-8 -*-
from __future__ import unicode_literals
from model_utils import Choices
from django.utils.translation import gettext_lazy as _, ngettext as __
from django.db import models
from django.urls import reverse

class Organization(models.Model):
    name = models.CharField(max_length=200)
    url = models.URLField(null=True)
    country = models.CharField(max_length=200, null=True)
    city = models.CharField(max_length=200, null=True)
    scopus_id = models.CharField(max_length=200, null=True)
    scopus_names = models.TextField(null=True)

    deprecated = models.BooleanField(default=False)
    index_updated = models.BooleanField(default=False)
    created_at = models.DateTimeField(auto_now_add=True)
    updated_at = models.DateTimeField(auto_now=True)

    class Meta:
        verbose_name_plural = _("Organizations")

    def rtype(self):
        return "org"

    def get_absolute_url(self):
        return reverse('bioresources:organization_view', args=[str(self.id)])

    def __str__(self):
        return " ".join([x for x in [self.name, "|", self.country, self.city] if x])