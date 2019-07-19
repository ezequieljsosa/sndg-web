# -*- coding: utf-8 -*-
from __future__ import unicode_literals
from model_utils import Choices
from django.utils.translation import gettext_lazy as _, ngettext as __
from django.db import models

class Affiliation(models.Model):
    publication = models.ForeignKey(Publication, on_delete=models.CASCADE, related_name="affiliations")
    author = models.ForeignKey(Person, on_delete=models.CASCADE, related_name="affiliations")
    organizations = models.ManyToManyField(Organization, related_name="affiliations")

    class Meta:
        verbose_name_plural = __("Affiliations")

    def __str__(self):
        return ("Affiliation: (%s) (%s) (%s) " % [str(x) for x in [self.author, self.publication] + [x.name for x in
                                                                                                     self.organizations.all()]])