# -*- coding: utf-8 -*-
from __future__ import unicode_literals
from django.utils.translation import gettext_lazy as __
from django.shortcuts import redirect, reverse
from django.shortcuts import render

from bioresources.models.ReadsArchive import ReadsArchive
from bioresources.io.GraphRepo import GraphRepo


def reads(request, pk):
    sra = ReadsArchive.objects.get(id=pk)
    graph, related_resources = GraphRepo.get_neighborhood(pk, "Reads",level=1)

    return render(request, 'resources/readsarchive.html',
                  {"readsarchive": sra, "sidebarleft": 1,"level":1,
                   "graph": graph, "related_resources": related_resources, "pk": pk, "rtype_src": "Reads"
                   })
