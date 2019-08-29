# -*- coding: utf-8 -*-
from __future__ import unicode_literals
from django.utils.translation import gettext_lazy as __
from django.shortcuts import redirect, reverse
from django.shortcuts import render

from bioseq.models.Bioentry import Bioentry
from bioseq.models.Seqfeature import Seqfeature
from bioseq.models.Biosequence import Biosequence
from bioresources.models.Assembly import Assembly

def NucleotideView(request,pk):
    be = (Bioentry.objects.select_related("biodatabase").select_related("taxon")
        .prefetch_related("dbxrefs__dbxref", "qualifiers__term", "seq")).get(bioentry_id=pk)

    offset = request.GET.get("offset",0)

    if be.features.count() < 100:
        be = be.prefetch_related("dbxrefs__dbxref", "qualifiers__term", "seq",
                                  "features__locations", "features__source_term",
                                  "features__type_term", "features__qualifiers")
        feature_list = be.features.all()
    else:
        feature_list = be.features.all()[offset:offset + request.GET.get("pageSize",50)]

    assembly = Assembly.objects.get(name=be.biodatabase.name)


    return render(request, 'resources/nucleotide_view.html', {
        "object": be,"assembly":assembly,"feature_list":feature_list,
        "sidebarleft": 0})  # "sidebarrigth": {"news": [{"title": "n1", "text": "lalala"}]