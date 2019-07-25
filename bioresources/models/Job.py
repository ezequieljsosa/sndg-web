# -*- coding: utf-8 -*-
from __future__ import unicode_literals
from model_utils import Choices
from django.utils.translation import gettext_lazy as _, ngettext as __
from django.db import models

from sndg.users.models import User
from bioresources.tasks import execute_job


class Job(models.Model):
    STATUS = Choices(
        *[(i, x, _(x)) for i, x in enumerate([
            "NEW", "QUEUED", "ERROR", "RETRYING", "FINISHED"
        ])]
    )

    command = models.TextField()
    start = models.DateTimeField(null=True)
    end = models.DateTimeField(null=True)
    result = models.TextField(null=True)
    status = models.PositiveIntegerField(default=STATUS.NEW, choices=STATUS)
    retry = models.PositiveIntegerField(default=0)
    user = models.ForeignKey(User, on_delete=models.CASCADE, related_name="jobs", null=True)
    dev_error = models.TextField(null=True)

    class Meta:
        verbose_name_plural = _("Job")

    def execute(self):
        execute_job.delay(self.id)
