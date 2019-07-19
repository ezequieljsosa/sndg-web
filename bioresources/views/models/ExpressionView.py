# -*- coding: utf-8 -*-
from __future__ import unicode_literals
from django.utils.translation import gettext_lazy as __
from django.shortcuts import redirect, reverse
from django.shortcuts import render

def expression(request, pk):
    expression = Expression.objects.get(id=pk)
    graph = {"nodes": [{"id": expression.name, "label": labelize(expression.name), "color": "orange"}], "edges": []}
    external_orgs = []
    publications_from_resource_graph(expression, graph, external_orgs)
    return render(request, 'resources/expression.html', {
        "expression": expression, "graph": graph, "external_orgs": external_orgs,
        "sidebarleft": 1, })