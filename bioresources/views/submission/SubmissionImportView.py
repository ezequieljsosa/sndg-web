# -*- coding: utf-8 -*-
from __future__ import unicode_literals
from django.utils.translation import gettext_lazy as __
from django.shortcuts import redirect, reverse
from django.shortcuts import render
from django.db import transaction

from bioresources.models.ExternalId import ExternalId
from bioresources.models.Organization import Organization


class LineResult():
    def __init__(self, line):
        self.line = line
        self.resources = []
        self.message = ""


class ResourceResult():
    def __init__(self, msg, resource):
        self.msg = msg
        self.resource = resource


def SubmissionImportView(request):
    current_user = request.user
    if request.method == 'POST':

        search = request.POST["search"]

        results = process_search_lines(search)

        # ncbi_search = NCBISearch()
        # page = Page.from_request(request)
        # records, total = ncbi_search.search(search, retmax=page.size, retstart=page.offset)
        # existing = [x for x in Assembly.objects.filter(name__in=[x.name for x in records])]
        #
        # if current_user:
        #     collaborated = [x.name for x in existing if
        #                     hasattr(current_user, "person") and current_user.person in x.collaborators.all()]
        # else:
        #     collaborated = []
        # collaborate = [x.name for x in existing if x.name not in collaborated]
        #
        # page.set_count(total)
        # return render(request, 'submission/import_resource.html',
        #               {'resource': "assembly", "search": search, "collaborate": collaborate,
        #                "page_obj": page, "results": records, "collaborated": collaborated})
        return render(request, 'submission/submission_import.html', {
            'resource': "assembly", "search": search, "results": results})
    else:
        search = request.GET.get("search", "")
        return render(request, 'submission/submission_import.html', {
            'resource': "assembly", "search": search})


def process_search_lines(search):
    results = []
    for l in search.split("\n"):
        if l.strip():
            results.append(process_line(l))
    return results


def process_line(line):
    r = LineResult(line)
    if len(line.split()) > 1:
        r.message = __("lines must not have any spaces")

    if any([x in line for x in
            ["www.ncbi.nlm.nih.gov", "www.ebi.ac.uk", "www.uniprot.org", "rnacentral.org", "ensembl.org"]]):
        return process_url(line)
    elif "http://" in line or "https://" in line:
        r.message = __("url not recognized")

    else:
        r = process_id(line.strip())
    return r


def process_id(rid):
    line_result = LineResult(rid)
    from bioresources.io.NCBISearch import NCBISearch
    ncbi_search = NCBISearch()
    counts = ncbi_search.search_all(rid)
    for db, count in counts.items():
        if count:
            if (count == 1) and (db not in ["taxonomy", "genome", "biosample"]):
                r = ExternalId.objects.prefetch_related("resource").filter(resource__type= NCBISearch.db_type[db] ,
                    organization=Organization.objects.get(name="NCBI"), identifier=rid)
                if r.count():
                    line_result.resources.append(ResourceResult(__("Already in DB"), r.get().resource))
                else:
                    with transaction.atomic():
                        r = ncbi_search.search_database(db, rid)[0][0]
                        line_result.resources.append(ResourceResult("", r))
                            #     line_result.message = __("Too many results with this identifier. Is there other bioproject or assembly that groups this resources?")

    return line_result
