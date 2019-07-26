# -*- coding: utf-8 -*-
from __future__ import unicode_literals
from django.utils.translation import gettext_lazy as __
from django.shortcuts import redirect, reverse
from django.shortcuts import render

from bioresources.models.Person import Person


def person(request, pk):
    person = Person.objects.get(id=pk)

    return render(request, 'resources/person.html', {
        "person": person, "sidebarleft": 1, })
