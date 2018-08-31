# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from haystack.generic_views import SearchView
from haystack.query import SearchQuerySet

from django.shortcuts import render

from django_tables2 import RequestConfig

from .models import Publication, Structure, Assembly, Expression, Person, Organization, Barcode, BioProject, \
    Affiliation, Resource
from .tables import PublicationTable
from .filters import PublicationFilter
from .forms import BioSearchForm

from collections import defaultdict

from django.shortcuts import redirect, reverse

from biosql.models import Tool, Biodatabase
from biosql.views import publications_from_resource_graph, labelize
import urllib
from haystack.inputs import Exact

resources = [
    {"name": "Genomas", "count": 99, "icon": "circle", "type": "genome"},
    {"name": "Proteinas", "count": 99, "icon": "puzzle-piece", "type": "protein"},
    {"name": "Herramientas", "count": 99, "icon": "wrench", "type": "tool"},
    {"name": "Estructuras", "count": 99, "icon": "sitemap", "type": "structure"},
    {"name": "Barcodes", "count": 99, "icon": "barcode", "type": "barcode"},
    {"name": "Ensamblados", "count": 99, "icon": "adn", "type": "assembly"},
    {"name": "Bioproyectos", "count": 99, "icon": "briefcase", "type": "bioproject"},
    {"name": "Publicaciones", "count": 99, "icon": "book", "type": "pmc"},
    {"name": "Organizaciones", "count": 99, "icon": "building", "type": "org"},
    {"name": "Investigadores", "count": 99, "icon": "user", "type": "person"},
    {"name": "Expresión", "count": 99, "icon": "sliders-h", "type": "gds"},

]


def index(request):
    db = request.GET.get("db", "all")
    search = request.GET.get("search", "").strip()

    selected = {x: Exact(request.GET[x]) for x in ["authors", "affiliations", "taxon"]
                if x in request.GET}
    if search:
        sqs = SearchQuerySet().filter(content=search, **selected).facet("type")
    else:
        sqs = SearchQuerySet().filter(**selected).facet("type")
        search = "*"

    params = dict(request.GET)
    params["q"] = [search]

    for ft in ["authors", "affiliations", "taxon"]:
        if ft not in selected:
            sqs = sqs.facet(ft, limit=5)

    facets = sqs.facet_counts()
    rdata = defaultdict(lambda: 0, {k: v for k, v in facets["fields"]["type"]})
    count = 0
    for r in resources:
        r["count"] = rdata[r["type"]]
        count += rdata[r["type"]]

    suggestions = []
    if count == 0:
        suggestions = SearchQuerySet().auto_query(search).spelling_suggestion()
        if suggestions:
            suggestions = [x.strip() for x in suggestions.replace("(", " ").split(")") if x.strip()]

    del facets["fields"]["type"]
    return render(request, 'index.html', {
        "stats": resources, "search": search, "selected": selected,
        "db": db, "suggestions": suggestions, "querystring": params,
        "sidebarleft": facets["fields"], "sidebarrigth": {"news": [{"title": "n1", "text": "lalala"}]}})


def structure(request, pk):
    return redirect(reverse("pdbdb:structure_view", kwargs={"pdbid": Structure.objects.get(id=pk).name}))


def assembly(request, pk):
    assembly = Assembly.objects.prefetch_related("external_ids").get(id=pk)
    accession = assembly.external_ids.filter(type="accession").first().identifier
    pk2 = Biodatabase.objects.get(name=accession).biodatabase_id
    return redirect(reverse("biosql:assembly_view", kwargs={"pk": pk2}))


def publication(request, pk):
    publication = Publication.objects.prefetch_related("targets__target").get(id=pk)
    orgs = publication.affiliation_names(all=True)
    authors = []
    for aff in publication.affiliations.all():
        author = aff.author.complete_name() + " " + " ".join(
            ["(" + str(orgs.index(org.name) + 1) + ")" for org in aff.organizations.all()])
        if author not in authors:
            authors.append(author)

    target_resources = []
    for x in publication.targets.all():
        x.target.url_ = reverse("bioresources:assembly_view", kwargs={"pk": x.target.id})
        target_resources.append(x.target)

    return render(request, 'resources/publication.html', {
        "publication": publication, "orgs": orgs, "authors": authors,
        "sidebarleft": 1, "target_resources": target_resources})


def barcode(request, pk):
    barcode = Barcode.objects.get(id=pk)
    return render(request, 'resources/barcode.html', {
        "barcode": barcode,
        "sidebarleft": 1, })


def organization(request, pk):
    organization = Organization.objects.get(id=pk)

    sqs = SearchQuerySet().filter(affiliations=Exact(organization.name)).facet("type")

    facets = sqs.facet_counts()
    facets = {k: v for k, v in facets["fields"]["type"] if v > 0}

    search_url = (reverse('bioresources:search_view') +
                  "?q=*&affiliations=" + organization.name + "&type=")
    return render(request, 'resources/organization.html', {
        "organization": organization, "facets": facets, "search_url": search_url,
        "sidebarleft": 1, "nameMap": {"person": "Investigadores", "pmc": "Publicaciones"}})


def person(request, pk):
    person = Person.objects.get(id=pk)

    sqs = SearchQuerySet().filter(authors=Exact(person.complete_name())).facet("type")

    facets = sqs.facet_counts()
    facets = {k: v for k, v in facets["fields"]["type"] if v > 0}

    search_url = (reverse('bioresources:search_view') +
                  "?q=*&authors=" + person.complete_name() + "&type=")

    return render(request, 'resources/person.html', {
        "person": person, "facets": facets, "search_url": search_url,
        "nameMap": {"person": "Investigadores", "pmc": "Publicaciones"},
        "sidebarleft": 1, })


def expression(request, pk):
    expression = Expression.objects.get(id=pk)
    graph = {"nodes": [{"id": expression.name, "label": labelize(expression.name), "color": "orange"}], "edges": []}
    external_orgs = []
    publications_from_resource_graph(expression, graph, external_orgs)
    return render(request, 'resources/expression.html', {
        "expression": expression, "graph": graph, "external_orgs": external_orgs,
        "sidebarleft": 1, })


def bioproject(request, pk):
    bioproject = BioProject.objects.get(id=pk)
    graph = {"nodes": [{"id": bioproject.name, "label": labelize(bioproject.name), "color": "orange"}], "edges": []}
    external_orgs = []
    publications_from_resource_graph(bioproject, graph, external_orgs)
    return render(request, 'resources/bioproject.html', {
        "bioproject": bioproject, "graph": graph, "external_orgs": external_orgs,
        "sidebarleft": 1, })


def tool(request, pk):
    tool = Tool.objects.get(id=pk)
    return render(request, 'resources/tool.html', {
        "tool": tool,
        "sidebarleft": 1, })


from django_filters.views import FilterView
from django_tables2.views import SingleTableMixin


class FilteredPersonListView(SingleTableMixin, FilterView):
    table_class = PublicationTable
    model = Publication
    template_name = 'publications/list.html'

    filterset_class = PublicationFilter
    paginate = {'per_page': 25}


class BioSearchView(SearchView):
    form_class = BioSearchForm

    def get_queryset(self):
        queryset = super(BioSearchView, self).get_queryset()
        for x in ["authors", "affiliations", "taxon"] + Resource.facet_dict.get(self.request.GET["type"], []):
            if x not in self.request.GET:
                queryset = queryset.facet(x, limit=5)

        return queryset  # .filter(pub_date__gte=date(2015, 1, 1))

    def get_context_data(self, *args, **kwargs):
        context = super(BioSearchView, self).get_context_data(*args, **kwargs)

        context["sidebarleft"] = {k: [(y1, y2) for y1, y2 in v if y2] for k, v in
                                  context["view"].queryset.query._facet_counts["fields"].items()}

        if hasattr(context["view"].queryset.query, "_raw"):
            context["suggestions"] = context["view"].queryset.query._raw.spellcheck["suggestions"]
        # dbtitle = self.request.GET["db"]
        # type
        # authors
        # affiliations
        context["dbtype"] = self.get_form_kwargs()["data"]["type"]

        for x in context["page_obj"].object_list:
            x.sndg_url = reverse("bioresources:" + context["dbtype"] + '_view', kwargs={"pk": x.object.id})

        prange = list(context["page_obj"].paginator.page_range)

        context["page_obj"].paginator.show_pages = (prange[max(
            context["page_obj"].number - 3, 0):context["page_obj"].number + 2])

        context["selected"] = {x: self.get_form_kwargs()["data"][x] for x in
                               ["authors", "affiliations", "taxon"] +
                               Resource.facet_dict.get(self.request.GET["type"], [])
                               if x in self.get_form_kwargs()["data"]}

        context["params"] = dict(self.request.GET)
        suggestion_list = []
        context["suggestion_list"] = suggestion_list
        if context["suggestions"]:
            for x in context["suggestions"][1::2]:
                for y in x["suggestion"]:
                    suggestion_list.append(y)

        # SearchQuerySet().filter(content="genomea") .facet('affiliations',limit=5).facet("type").facet_counts()
        return context

# def publications(request):
#     table = PublicationTable(Publication.objects.all())
#     RequestConfig(request,paginate={'per_page': 25}).configure(table)
#     return render(request, 'publications/list.html',
#                   {"publications": table})
