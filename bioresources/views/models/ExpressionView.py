# -*- coding: utf-8 -*-
from __future__ import unicode_literals
from django.utils.translation import gettext_lazy as __
from django.shortcuts import redirect, reverse
from django.shortcuts import render

from bioseq.views import labelize

from bioresources.models.Expression import Expression
from bioresources.io.GraphRepo import GraphRepo


def expression(request, pk):
    expression = Expression.objects.get(id=pk)
    graph, related_resources = GraphRepo.get_neighborhood(pk, "Expression",level=2)


    return render(request, 'resources/expression.html', {
        "expression": expression, "sidebarleft": 1,
        "graph": graph, "related_resources": related_resources, "pk": pk, "rtype_src": "Expression"
    })