# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.core.exceptions import ImproperlyConfigured
from django.core.files.storage import FileSystemStorage
from django.forms import Form
from django.http import HttpResponseRedirect, HttpResponse
from haystack.generic_views import SearchView
from haystack.query import SearchQuerySet

from django.shortcuts import render

from django_tables2 import RequestConfig

from .models import Publication, Structure, Assembly, Expression, Person, Organization, Barcode, BioProject, \
    Affiliation, Resource, Sample
from .tables import PublicationTable
# from .filters import PublicationFilter
from .forms import BioSearchForm, AssemblyForm

from collections import defaultdict

from django.shortcuts import redirect, reverse

from biosql.models import Tool, Biodatabase
from biosql.views import labelize
import urllib
from haystack.inputs import Exact

from biosql.io.Pagination import Page

from django.contrib.auth.decorators import login_required
from django.utils.translation import gettext as _
from django.utils.translation import gettext_lazy as __

# noinspection PyProtectedMember
resources = [
    # {"name": "Genomas", "count": 99, "icon": "circle", "type": "genome"},
    {"name": __("Proteins"), "count": 99,
     "icon": "puzzle-piece",
     "type": Resource.RESOURCE_TYPES.PROTEIN},
    {"name": __("Tools"), "count": 99, "icon": "wrench",
     "type": Resource.RESOURCE_TYPES.TOOL},
    {"name": Structure._meta.verbose_name_plural, "count": 99,
     "icon": "sitemap",
     "type": Resource.RESOURCE_TYPES.STRUCTURE},
    {"name": Barcode._meta.verbose_name_plural, "count": 99,
     "icon": "barcode",
     "type": Resource.RESOURCE_TYPES.BARCODE},
    {"name": Assembly._meta.verbose_name_plural, "count": 99,
     "icon": "adn",
     "type": Resource.RESOURCE_TYPES.ASSEMBLY},
    {"name": BioProject._meta.verbose_name_plural, "count": 99,
     "icon": "briefcase",
     "type": Resource.RESOURCE_TYPES.BIOPROJECT},
    {"name": Publication._meta.verbose_name_plural, "count": 99,
     "icon": "book",
     "type": Resource.RESOURCE_TYPES.PUBLICATION},
    {"name": __("Organizations"), "count": 99, "icon": "building", "type": "org"},
    {"name": __("Persons"), "count": 99, "icon": "user", "type": "person"},
    {"name": Expression._meta.verbose_name_plural, "count": 99,
     "icon": "sliders-h",
     "type": Resource.RESOURCE_TYPES.EXPRESSION},

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
    if "fields" in facets:
        rdata = defaultdict(lambda: 0, {k: v for k, v in facets["fields"]["type"]})
    else:
        rdata = defaultdict(lambda: 0)
    count = 0
    for r in resources:
        r["count"] = rdata[str(r["type"])]
        count += rdata[r["type"]]

    suggestions = []
    if count == 0:
        suggestions = SearchQuerySet().auto_query(search).spelling_suggestion()
        if suggestions:
            suggestions = [x.strip() for x in suggestions.replace("(", " ").split(")") if x.strip()]

    if "fields" in facets:
        del facets["fields"]["type"]
    else:
        facets["fields"] = {}

    return render(request, 'index.html', {
        "stats": resources, "search": search, "selected": selected,
        "db": db, "suggestions": suggestions, "querystring": params,
        "sidebarleft": facets["fields"], "sidebarrigth": {"news": [{"title": "n1", "text": "lalala"}]}})


def structure(request, pk):
    return redirect(reverse("pdbdb:structure_view", kwargs={"pdbid": Structure.objects.get(id=pk).name}))


def assembly(request, pk):
    from xml.sax.saxutils import escape
    try:
        pk = int(pk)
        assembly = Assembly.objects.prefetch_related("external_ids").get(id=pk)
        sqs = assembly.external_ids.filter(type="accession")
        accession = sqs.first().identifier if sqs.exists() else assembly.name
    except ValueError:
        accession = pk
        # assembly = (Assembly.objects.prefetch_related("external_ids").filter(external_ids__type="accession",
        #                                                                      external_ids__identifier=pk))

    pk2 = Biodatabase.objects.get(name=accession).biodatabase_id

    return redirect(reverse("biosql:assembly_view", kwargs={"pk": pk2}))


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
        "sidebarleft": 1, "nameMap":
            {"person": Person._meta.verbose_name_plural,
             Resource.RESOURCE_TYPES.PUBLICATION: Publication._meta.verbose_name_plural}})


def person(request, pk):
    person = Person.objects.get(id=pk)

    sqs = SearchQuerySet().filter(authors=Exact(person.complete_name())).facet("type")

    facets = sqs.facet_counts()
    facets = {k: v for k, v in facets["fields"]["type"] if v > 0}

    search_url = (reverse('bioresources:search_view') +
                  "?q=*&authors=" + person.complete_name() + "&type=")

    return render(request, 'resources/person.html', {
        "person": person, "facets": facets, "search_url": search_url,
        "nameMap": {"person": Person._meta.verbose_name_plural,
                    Resource.RESOURCE_TYPES.PUBLICATION: Publication._meta.verbose_name_plural},
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


def reads(request, pk):
    tool = Tool.objects.get(id=pk)
    return render(request, 'resources/reads.html', {
        "tool": tool,
        "sidebarleft": 1, })


def submission(request):
    return render(request, "forms/submission.html")


# from django_filters.views import FilterView
# from django_tables2.views import SingleTableMixin


# class FilteredPersonListView(SingleTableMixin, FilterView):
#     table_class = PublicationTable
#     model = Publication
#     template_name = 'publications/list.html'
#
#     filterset_class = PublicationFilter
#     paginate = {'per_page': 25}


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
        if "suggestions" in context and context["suggestions"]:
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

# @login_required
def assembly_new(request):
    if request.method == 'POST':
        # create a form instance and populate it with data from the request:
        form = AssemblyForm(request.POST)

        if form.is_valid():
            assembly = form.save()

            return HttpResponseRedirect('/assembly/' + str(assembly.id))


    else:
        form = AssemblyForm()

    return render(request, 'forms/assembly_new.html', {'form': form})


def import_resource(request):
    current_user = request.user
    if request.method == 'POST':
        from .io.NCBISearch import NCBISearch
        search = request.POST["search"]
        ncbi_search = NCBISearch()
        page = Page.from_request(request)
        records, total = ncbi_search.search(search, retmax=page.size, retstart=page.offset)
        existing = [x for x in Assembly.objects.filter(name__in=[x.name for x in records])]

        if current_user:
            collaborated = [x.name for x in existing if
                           hasattr(current_user, "person") and current_user.person in x.collaborators.all()]
        else:
            collaborated = []
        collaborate = [x.name for x in existing if x.name not in collaborated]

        page.set_count(total)
        return render(request, 'forms/import_resource.html',
                      {'resource': "assembly", "search": search, "collaborate": collaborate,
                       "page_obj": page, "results": records, "collaborated": collaborated})
    else:
        search = request.GET.get("search", "")
        return render(request, 'forms/import_resource.html', {
            'resource': "assembly", "search": search})


# from django_fine_uploader.views import FineUploaderView
# class NotConcurrentUploaderView(FineUploaderView):
#     """Example of a chunked, but NOT concurrent upload.
#     Disabling concurrent uploads per view.
#     Remember, you can turn off concurrent uploads on your settings, with:
#     FU_CONCURRENT_UPLOADS = False
#     """
#     @property
#     def concurrent(self):
#         return False
#
#     def form_valid(self, form):
#         self.process_upload(form)
#         return self.make_response({'success': True})
#
from django.views import generic


class UploadView(generic.TemplateView):
    template_name = 'forms/resource_upload.html'


from django.conf import settings

from resumable.fields import ResumableFileField


class ResumableForm(Form):
    file = ResumableFileField(
        # allowed_mimes=("audio/ogg",),
        upload_url=lambda: reverse('upload_api'),
        chunks_dir=getattr(settings, 'FILE_UPLOAD_TEMP_DIR')
    )


from resumable.files import ResumableFile


def upload_view(request):
    form = ResumableForm()
    return render(request, 'forms/resource_upload.html', {'form': form})


class ResumableUploadView(generic.TemplateView):
    template_name = 'forms/resource_upload.html'

    def get(self, request, *args, **kwargs):
        if "resumableChunkNumber" in request.GET:
            r = ResumableFile(self.storage, self.request.GET)
            if not (r.chunk_exists or r.is_complete):
                return HttpResponse('chunk not found', status=404)
            return HttpResponse('chunk already exists')
        else:
            form = ResumableForm()
            return render(request, self.template_name, {'form': form})

    def post(self, *args, **kwargs):
        """Saves chunks then checks if the file is complete.
        """
        chunk = self.request.FILES.get('file')
        r = ResumableFile(self.storage, self.request.POST)
        if r.chunk_exists:
            return HttpResponse('chunk already exists')
        r.process_chunk(chunk)
        if r.is_complete:
            self.process_file(r.filename, r)
            r.delete_chunks()
        return HttpResponse()

    def process_file(self, filename, file):
        """Process the complete file.
        """
        self.storage.save(filename, file)

    @property
    def chunks_dir(self):
        chunks_dir = getattr(settings, 'FILE_UPLOAD_TEMP_DIR', None)
        if not chunks_dir:
            raise ImproperlyConfigured(
                'You must set settings.FILE_UPLOAD_TEMP_DIR')
        return chunks_dir

    @property
    def storage(self):
        return FileSystemStorage(location=self.chunks_dir)


from django.contrib.auth.decorators import login_required


def sample_view(request, pk):
    sample = Sample.objects.get(id=pk)
    graph = {"nodes": [{"id": expression.name, "label": labelize(expression.name), "color": "orange"}], "edges": []}
    external_orgs = []
    publications_from_resource_graph(sample, graph, external_orgs)
    return render(request, 'resources/expression.html', {
        "expression": expression, "graph": graph, "external_orgs": external_orgs,
        "sidebarleft": 1, })


def assembly_view(request, pk):
    be = Biodatabase.objects.get(biodatabase_id=pk)
    if "_prots" in be.name:
        be2 = (Biodatabase.objects.prefetch_related("entries__features__type_term"))
        be = be2.get(name=be.name.replace("_prots", ""))  # TODO meter en otro lado esta regla
    else:
        be = (Biodatabase.objects.prefetch_related("entries__features__type_term"))
        be = be.get(biodatabase_id=pk)

    assembly = Assembly.objects.get(external_ids__identifier=be.name)
    lengths = {}
    for x in be.entries.all():
        seq = Biosequence.objects.raw("""
        SELECT bioentry_id, version , length , alphabet ,"" seq
        FROM biosequence WHERE bioentry_id = %i ;
        """ % (x.bioentry_id))[0]
        lengths[x.accession] = seq.length

    graph = assembly_graph(assembly)
    return render(request, 'biosql/assembly_detail.html', {"lengths": lengths,
                                                           "object": assembly, "graph": graph,
                                                           "contigs": be.entries.all(),
                                                           "sidebarleft": {}})


def sequence_view(request, pk):
    be = (Bioentry.objects.select_related("biodatabase").select_related("taxon")
        .prefetch_related("dbxrefs__dbxref", "qualifiers__term", "seq",
                          "features__locations", "features__source_term",
                          "features__type_term", "features__qualifiers"))
    be = be.get(bioentry_id=pk)

    if be.biodatabase.name.endswith("prots"):
        beg = Biodatabase.objects.get(name=be.biodatabase.name.replace("_prots", ""))
        taxon = beg.entries.first().taxon

        feature = Seqfeature.objects.seqfeature_from_locus_tag(beg.biodatabase_id, be.accession)
        feature = list(feature)[0]

        locations = list(feature.locations.all())
        start = locations[0].start_pos
        end = locations[-1].end_pos

        seq = Biosequence.objects.raw("""
        SELECT bioentry_id, version , length , alphabet ,SUBSTRING( seq,%i,%i ) seq
        FROM biosequence WHERE bioentry_id = %i ;
        """ % (start, end - start, feature.bioentry_id))[0]
        functions = {"biological_process": [], "molecular_function": [], "cellular_component": []}
        for qual in be.qualifiers.filter(term__dbxrefs__dbxref__dbname="go",
                                         term__dbxrefs__dbxref__accession="goslim_generic"):
            for dbxref in qual.term.dbxrefs.all():
                if dbxref.dbxref.accession in functions:
                    functions[dbxref.dbxref.accession].append(qual.term)

    else:
        beg = be.biodatabase

    graph = entry_graph(be, beg)

    if be.biodatabase.name.endswith("prots"):
        return render(request, 'biosql/protein_detail.html', {
            "functions": functions, "graph": graph, "accession": be.biodatabase.name.replace("_prots", ""),
            "object": be, "feature": feature, "taxon": taxon, "seq": seq, "start": start, "end": end,
            "sidebarleft": 1})
    else:
        return render(request, 'biosql/biosequence_detail.html', {
            "object": be, "graph": graph,
            "sidebarleft": 0})  # "sidebarrigth": {"news": [{"title": "n1", "text": "lalala"}]
