# -*- coding: utf-8 -*-
from __future__ import unicode_literals
from django.utils.translation import gettext_lazy as __
from django.shortcuts import redirect, reverse
from django.shortcuts import render

from bioresources.models.Tool import Tool

from bioresources.io.GraphRepo import GraphRepo


def tool(request, pk):
    tool = Tool.objects.get(id=pk)
    graph, related_resources = GraphRepo.get_neighborhood(pk, "Tool", level=1)

    return render(request, 'resources/tool.html', {"external_url":tool.url,
        "graph": graph, "related_resources": related_resources, "pk": pk, "rtype_src": "Tool",
        "tool": tool, "sidebarleft": 1, })
