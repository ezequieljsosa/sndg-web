# -*- coding: utf-8 -*-
from __future__ import unicode_literals
from django.utils.translation import gettext_lazy as __
from django.shortcuts import redirect, reverse
from django.shortcuts import render

from bioseq.models.Biodatabase import Biodatabase
from bioseq.models.Biosequence import Biosequence

from bioresources.models.Assembly import Assembly

# def assembly(request, pk):
#     try:
#         pk = int(pk)
#         assembly = Assembly.objects.prefetch_related("external_ids").get(id=pk)
#         sqs = assembly.external_ids.filter(type="accession")
#         accession = sqs.first().identifier if sqs.exists() else assembly.name
#     except ValueError:
#         accession = pk
#         # assembly = (Assembly.objects.prefetch_related("external_ids").filter(external_ids__type="accession",
#         #                                                                      external_ids__identifier=pk))
#
#     pk2 = Biodatabase.objects.get(name=accession).biodatabase_id
#
#     return redirect(reverse("bioseq:assembly_view", kwargs={"pk": pk2}))


def assembly_view(request, pk):
    assembly = Assembly.objects.get(id=pk)

    be = Biodatabase.objects.get(name=assembly.name)
    if "_prots" in be.name:
        be2 = (Biodatabase.objects.prefetch_related("entries__features__type_term"))
        be = be2.get(name=be.name.replace("_prots", ""))  # TODO meter en otro lado esta regla

    lengths = {}
    for x in be.entries.all():
        seq = Biosequence.objects.raw("""
        SELECT bioentry_id, version , length , alphabet 
        FROM biosequence WHERE bioentry_id = %i ;
        """ % (x.bioentry_id))[0]
        lengths[x.accession] = seq.length
    assembly.assembly_type = Assembly.ASSEMBLY_TYPES[assembly.assembly_type]
    assembly.level = Assembly.ASSEMBLY_LEVEL[assembly.level]



    # graph = assembly_graph(assembly)
    return render(request, 'resources/assembly_detail.html', {"lengths": lengths,
                                                           "object": assembly, "graph": {},
                                                           "contigs": be.entries.all(),
                                                           "sidebarleft": {}})