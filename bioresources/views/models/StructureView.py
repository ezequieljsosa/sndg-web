# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.utils.translation import gettext_lazy as __
from django.shortcuts import redirect, reverse
from django.shortcuts import render

from bioresources.models.Structure import Structure
from bioresources.io.GraphRepo import GraphRepo


def structure(request, pk):
    pdb = Structure.objects.get(id=pk)
    graph, related_resources = GraphRepo.get_neighborhood(pk, "Structure", level=1)
    external_url = "https://www.rcsb.org/structure/" + pdb.name
    return render(request, 'resources/structure.html', {"external_url":external_url,
        "graph": graph, "related_resources": related_resources, "pk": pk, "rtype_src": "Structure",
        "pdb": pdb, "sidebarleft": 1, "level": 1})
