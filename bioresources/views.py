# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from haystack.generic_views import SearchView
from haystack.query import SearchQuerySet

from django.shortcuts import render

from django_tables2 import RequestConfig

from .models import Publication
from .tables import PublicationTable
from .filters import PublicationFilter
from .forms import BioSearchForm

from collections import defaultdict

from django.shortcuts import redirect, reverse

resources = [
    {"name": "Genomas", "count": 99, "icon": "circle", "type": "genome"},
    {"name": "Proteinas", "count": 99, "icon": "puzzle-piece", "type": "protein"},
    {"name": "Herramientas", "count": 99, "icon": "wrench", "type": "tool"},
    {"name": "Estructuras", "count": 99, "icon": "sitemap", "type": "structure"},
    {"name": "Barcodes", "count": 99, "icon": "barcode", "type": "barcode"},
    {"name": "Ensamblados", "count": 99, "icon": "adn", "type": "seq"},
    {"name": "Bioproyectos", "count": 99, "icon": "briefcase", "type": "bioproject"},
    {"name": "Publicaciones", "count": 99, "icon": "book", "type": "pmc"},

]


def index(request):
    db = request.GET.get("db", "all")
    search = request.GET.get("search", "").strip()

    if db != "all":
        return redirect(reverse("bioresources:search_view") + "?q=" + search)

    if search:
        sqs = SearchQuerySet().filter(content=request.GET["search"]). \
            facet("authors", limit=5).facet('affiliations', limit=5). \
            facet("type").facet("taxon", limit=5)
    else:
        sqs = SearchQuerySet().facet("type")
    facets = sqs.facet_counts()
    rdata = defaultdict(lambda: 0, {k: v for k, v in facets["fields"]["type"]})
    count = 0
    for r in resources:
        r["count"] = rdata[r["type"]]
        count += rdata[r["type"]]

    suggestions = []
    if count == 0:
        suggestions = SearchQuerySet().auto_query(request.GET["search"]).spelling_suggestion()
        suggestions = [x.strip() for x in suggestions.replace("(", " ").split(")") if x.strip()]

    del facets["fields"]["type"]
    return render(request, 'index.html', {
        "stats": resources, "search": search,
        "db": db, "suggestions": suggestions,
        "sidebarleft": facets["fields"], "sidebarrigth": {"news": [{"title": "n1", "text": "lalala"}]}})


def publication(request, publication_id):
    publication = Publication.objects.get(id=publication_id)
    orgs = publication.affiliation_names()
    authors = []
    for aff in publication.affiliations.all():
        author = aff.author.complete_name() + " " + " ".join(
            ["(" + str(orgs.index(org.name) + 1) + ")" for org in aff.organizations.all()])
        authors.append(author)

    for x in publication.sources.all():
        print(x)
    for x in publication.targets.all():
        print(x)

    return render(request, 'resources/publication.html', {
        "publication": publication, "orgs": orgs, "authors": authors,
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
        for x in ["authors","affiliations","taxon"]:
            if x not in self.request.GET:
                queryset = queryset.facet(x, limit=5)

        return queryset  # .filter(pub_date__gte=date(2015, 1, 1))

    def get_context_data(self, *args, **kwargs):
        context = super(BioSearchView, self).get_context_data(*args, **kwargs)



        context["sidebarleft"] = context["view"].queryset.query._facet_counts["fields"]

        if hasattr(context["view"].queryset.query, "_raw"):
            context["suggestions"] = context["view"].queryset.query._raw.spellcheck["suggestions"]
        # dbtitle = self.request.GET["db"]
        # type
        # authors
        # affiliations
        context["dbtype"] = self.get_form_kwargs()["data"]["type"]
        context["selected"] = {x : self.get_form_kwargs()["data"][x]  for x in ["authors","affiliations","taxon"]
                if x  in self.get_form_kwargs()["data"] }
        # SearchQuerySet().filter(content="genomea") .facet('affiliations',limit=5).facet("type").facet_counts()
        return context

# def publications(request):
#     table = PublicationTable(Publication.objects.all())
#     RequestConfig(request,paginate={'per_page': 25}).configure(table)
#     return render(request, 'publications/list.html',
#                   {"publications": table})
