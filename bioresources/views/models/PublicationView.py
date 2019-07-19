# -*- coding: utf-8 -*-
from __future__ import unicode_literals
from django.utils.translation import gettext_lazy as __
from django.shortcuts import redirect, reverse
from django.shortcuts import render


def publication(request, pk):
    publication = Publication.objects.prefetch_related("targets__target").get(id=pk)
    orgs = publication.affiliation_names()
    authors = []
    for aff in publication.affiliations.all():
        author = aff.author.complete_name() + " " + " ".join(
            ["(" + str(orgs.index(org.name) + 1) + ")" for org in aff.organizations.all()])
        if author not in authors:
            authors.append(author)

    return render(request, 'resources/publication.html', {
        "publication": publication,
        "sidebarleft": 1})