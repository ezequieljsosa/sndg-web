# -*- coding: utf-8 -*-
from __future__ import unicode_literals
from django.utils.translation import gettext_lazy as __
from django.shortcuts import redirect, reverse
from django.shortcuts import render

from bioresources.models.BioProject import BioProject
from bioresources.views import labelize

def bioproject(request, pk):
    bioproject = BioProject.objects.get(id=pk)
    gr = GraphRepo()
    graph = {"nodes": [{"id": bioproject.name, "label": labelize(bioproject.name), "color": "orange"}], "edges": []}
    external_orgs = []
    # publications_from_resource_graph(bioproject, graph, external_orgs)
    return render(request, 'resources/bioproject.html', {
        "bioproject": bioproject, "graph": {}, "external_orgs": external_orgs,
        "sidebarleft": 1, })