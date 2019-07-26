# -*- coding: utf-8 -*-
from __future__ import unicode_literals
from django.utils.translation import gettext_lazy as __
from django.shortcuts import redirect, reverse
from django.shortcuts import render

from bioresources.models.Organization import Organization

def organization(request, pk):
    org = Organization.objects.get(id=pk)

    return render(request, 'resources/organization.html', {
    "organization": org, "sidebarleft": 1, })
