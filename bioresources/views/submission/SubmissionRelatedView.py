# -*- coding: utf-8 -*-
from __future__ import unicode_literals
from django.utils.translation import gettext_lazy as __
from django.shortcuts import redirect, reverse
from django.shortcuts import render
from django.db import transaction

from bioresources.models.ExternalId import ExternalId
from bioresources.models.Organization import Organization
from bioresources.models.Resource import Resource

from bioresources.models.Person import Person
from bioresources.models.BioProject import BioProject
from bioresources.models.Sample import Sample
from bioresources.models.Expression import Expression
from bioresources.models.Publication import Publication


class LineResult():
    def __init__(self, line):
        self.line = line
        self.resources = []
        self.message = ""


class ResourceResult():
    def __init__(self, msg, resource):
        self.msg = msg
        self.resource = resource


resource_class = {
    Resource.RESOURCE_TYPES.PERSON: Person,
    Resource.RESOURCE_TYPES.BIOPROJECT: BioProject,
    Resource.RESOURCE_TYPES.SAMPLE: Sample,
    Resource.RESOURCE_TYPES.EXPRESSION: Expression,
    Resource.RESOURCE_TYPES.PUBLICATION: Publication,
}


def SubmissionRelatedView(request):
    current_user = request.user
    if request.method == 'POST':
        Resource.objects.filter(relate)
        resource_type = request.POST["resource_type"]
        resource_id = request.POST["resource_id"]

        process_search = getattr(sys.modules[__name__], Resource.RESOURCE_TYPES[resource_type].lower() + "_process" )
        process_search(resource_class[resource_type].objects.get(id=resource_id))

        return render(request, 'submission/submission_related.html', {
            'resource': "assembly", "search": search, "resources": results})
    else:
        search = request.GET.get("search", "")
        return render(request, 'submission/submission_related.html', {
            'resource': "assembly", "search": search})



