# -*- coding: utf-8 -*-
from __future__ import unicode_literals
from django.utils.translation import gettext_lazy as __
from django.shortcuts import redirect, reverse
from django.shortcuts import render
from django.db import transaction

from bioresources.models.ExternalId import ExternalId
from bioresources.models.Organization import Organization
from bioresources.models.Resource import Resource, Collaboration

from bioresources.models.Person import Person
from bioresources.models.BioProject import BioProject
from bioresources.models.Sample import Sample
from bioresources.models.Expression import Expression
from bioresources.models.Publication import Publication
from bioresources.models.ResourceRelation import ResourceRelation

from bioresources.graph import connect_nodes

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

from django.contrib.auth.decorators import login_required


@login_required
def mark_to_relate(request, resource_id):
    request.session["relate_with"] = resource_id
    return redirect(request.POST["next"])


@login_required
def claim_identity(request, person_id):
    p = Person.objects.get(id=person_id)
    if request.method == 'POST':
        request.user.person = p
        request.user.save()
        return redirect(reverse('bioresources:user_resources'))
    else:
        return render(request, 'submission/claim_identity.html', {'person': p})


@login_required
def claim_resource(request, resource_id):

    r = Resource.objects.get(id=resource_id)


    if request.method == 'POST':
        person = request.user.person
        c = Collaboration(person=person, resource=r,
                          type=request.POST["relation"])
        c.save()
        person.type = Resource.RESOURCE_TYPES.PERSON
        x = connect_nodes(person,r,reltype=str(Collaboration.COLLABORATION_TYPES[int(c.type)]))
        return redirect(r.get_absolute_url())
    else:
        return render(request, 'submission/claim_resource.html', {
            "collaboration_types": {x[0]: str(x[1]) for x in
                                    Collaboration.COLLABORATION_TYPES},
            'resource': r})


@login_required
def SubmissionRelatedView(request, src_id, dst_id):
    r1 = Resource.objects.get(id=src_id)
    r2 = Resource.objects.get(id=dst_id)

    if request.method == 'POST':
        rr = ResourceRelation(source=r1,target=r2,role="uses")
        rr.save()
        connect_nodes(r1,r2)

        return redirect(r2.get_absolute_url())
    else:

        return render(request, 'submission/submission_related.html', {
            'resource1': r1, "resource2": r2})
