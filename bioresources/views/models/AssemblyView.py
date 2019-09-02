# -*- coding: utf-8 -*-
from __future__ import unicode_literals
from django.utils.translation import gettext_lazy as __
from django.shortcuts import redirect, reverse
from django.shortcuts import render

from bioseq.models.Biodatabase import Biodatabase
from bioseq.models.Biosequence import Biosequence

from bioresources.models.Assembly import Assembly
from bioresources.io.GraphRepo import GraphRepo
from bioresources.io.NCBISearch import NCBISearch
from bioresources.models.jobs.LoadGenomeJob import LoadGenomeJob
from bioresources.models.Job import Job

from bioseq.io.Pagination import Page


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
    graph, related_resources = GraphRepo.get_neighborhood(pk, "Barcode", level=1)

    loaded = True
    job = LoadGenomeJob.objects.filter(assembly=assembly)
    if job.exists():
        job = job.order_by("-id")[0]
        if job.status != Job.STATUS.FINISHED:
            loaded = False

    contigs = None
    lengths = None
    bdb = Biodatabase.objects.filter(name=assembly.name)

    external_ids = [x.identifier for x in assembly.external_ids.all() if x.type == "accession"]
    external_url = ""
    if external_ids:
        external_url = ("https://www.ncbi.nlm.nih.gov/" + NCBISearch.rtype2ncbb[Assembly.TYPE] + "/" + external_ids[
            0])

    if bdb.count():
        bdb = bdb.get()
        beqs = bdb.entries.all()

        lengths = {}
        # for x in be.entries.all():
        # SELECT s.bioentry_id, s.version , s.length , s.alphabet
        seqs = Biosequence.objects.prefetch_related("bioentry").raw("""
            SELECT s.bioentry_id, s.length
            FROM biosequence s,bioentry b WHERE b.biodatabase_id = %i AND  b.bioentry_id = s.bioentry_id ;
            """ % (bdb.biodatabase_id))
        for seq in seqs:
            lengths[seq.bioentry.accession] = seq.length
        # assembly.assembly_type = str(Assembly.ASSEMBLY_TYPES[assembly.assembly_type])
        # assembly.level = str(Assembly.ASSEMBLY_LEVEL[assembly.level])

        page = Page.from_request(request, count=beqs.count())

        contigs = beqs[page.offset(): page.end()]

    params = {"query": "", "page_obj": page,
              "lengths": lengths, "external_url": external_url, "loaded": loaded,
              "level": {k: str(v) for k, v in Assembly.ASSEMBLY_LEVEL},
              "object": assembly, "graph": graph, "atypes": {k: str(v) for k, v in Assembly.ASSEMBLY_TYPES},
              "related_resources": related_resources,
              "contigs": contigs, "sidebarleft": {}}

    return render(request, 'resources/assembly_detail.html', params)
