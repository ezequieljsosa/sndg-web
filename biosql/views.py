import math
from collections import defaultdict

from django.db.models import Prefetch
from django.shortcuts import render
from django.views.generic import TemplateView

from bioresources.models import Assembly, Publication, Sample, Resource
from .models import Taxon, Biosequence, Bioentry, Seqfeature, Biodatabase, BioentryQualifierValue, Ontology
from SNDG import Struct


class Page:
    def __init__(self):
        self.count = None
        self.size = None
        self.offset = None
        self.page = None
        self.number = None
        self.previous_page_number = None
        self.next_page_number = None
        self.has_previous = None
        self.has_next = None
        self.num_pages = None
        self.paginator = {}


class TableAttribute:
    def __init__(self, name, entity, render, query, url=None, value=None):
        self.name = name
        self.entity = entity
        self.render = render
        self.query = query
        self.url = url
        self.value = value


class AboutView(TemplateView):
    template_name = "about.html"

    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)
        return context


class TaxView(TemplateView):
    model = Taxon
    template_name = "biosql/taxon_detail.html"

    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)
        context['object'] = Taxon.objects.prefetch_related("names").get(ncbi_taxon_id=self.kwargs["pk"])
        return context


def labelize(long_string, size=4):
    idx = 0
    lines = []
    line = ""
    for x in long_string.split():
        idx += 1
        line += x + " "
        if idx == size:
            lines.append(line)
            line = " "
            idx = 0
    if line:
        lines.append(line)
    return "\\n".join(lines)


def sequence_lt_view(request, locus_tag):
    return sequence_view(request, Bioentry.objects.get(accession=locus_tag).bioentry_id)


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

    if be.biodatabase.name.endswith("prots"):
        return render(request, 'biosql/protein_detail.html', {
            "functions": functions, "accession": be.biodatabase.name.replace("_prots", ""),
            "object": be, "feature": feature, "taxon": taxon, "seq": seq, "start": start, "end": end,
            "sidebarleft": 1})
    else:
        return render(request, 'biosql/biosequence_detail.html', {
            "object": be,
            "sidebarleft": 0})  # "sidebarrigth": {"news": [{"title": "n1", "text": "lalala"}]


def filter_from_req(request):
    Prefetch("targets__source",
             queryset=qs)

    filter_map = {
        "tax": {
            "2query": lambda v: Q(bioentry__taxon__ncbi_taxon_id=v)
        }, "assembly": {
            "2query": lambda v: Q(bioentry__biodatabase_id=v)
        }, "ftype": {
            "2query": lambda v: Q(name=v) if v else (Q(name="CDS") | Q(name="tRNA") | Q(name="repeat_region"))
        }, "gene": {

        }}
    filters = []
    for fname in filter_map:
        if fname in request.GET:
            f = filter_map[fname]
            f["value"] = request.GET[fname]
            f["display"] = fname
            filters.append(f)

    return filters


def page_from_request(count, request):
    page_size = int(request.GET.get("page_size", "20"))
    if page_size > 500:
        page_size = 500
    page_num = int(request.GET.get("page", "1"))
    page_offset = (page_num - 1) * page_size

    prange = range(1, math.ceil(count / page_size))
    number = math.floor(page_offset / page_size) + 1
    page = Page()
    page.count = count
    page.size = page_size
    page.offset = page_offset
    page.page = page_num
    page.number = number
    page.previous_page_number = number - 1
    page.next_page_number = number + 1
    page.has_previous = page_num > 1
    page.has_next = page_num < math.ceil(count / page_size)
    page.num_pages = math.ceil(count / page_size)
    page.paginator = {"show_pages": (prange[max(
        number - 3, 0):number + 2]),
                      "num_pages": math.ceil(count / page_size)}

    return page


def gene_product_list_view(query_set, request):
    pass


def protein_list_view(request, assembly_id):
    # filters = filter_from_req(request)

    # TableAttribute(name,entity,render,query,url,value)
    be = Biodatabase.objects.get(biodatabase_id=assembly_id)
    if "_prots" not in be.name:
        be = Biodatabase.objects.get(name=be.name + "_prots")
        assembly_id = be.biodatabase_id
    # taxon = beg.entries.first().taxon

    filters = [TableAttribute("assembly", "bioentry", lambda e: e.biodatabase.name,
                              lambda param: Q(biodatabase_id=param), url=lambda e: e.biodatabase.get_absolute_url(),
                              value=Biodatabase.objects.get(biodatabase_id=assembly_id))]
    scorers = []

    def description(entity):
        beg = Biodatabase.objects.get(name=entity.biodatabase.name.replace("_prots", ""))
        f = Seqfeature.objects.seqfeature_from_locus_tag(beg.biodatabase_id, entity.accession)
        f = list(f)[0]
        return f.qualifiers.filter(term__name__in=["product","description","note"]).first().value

    def location(entity):
        beg = Biodatabase.objects.get(name=entity.biodatabase.name.replace("_prots", ""))
        f = Seqfeature.objects.seqfeature_from_locus_tag(beg.biodatabase_id, entity.accession)
        f = list(f)[0]
        locations = list(f.locations.all())
        start = locations[0].start_pos
        end = locations[-1].end_pos
        return {"start": start, "end": end, "ref": f.bioentry.name}

    def genes(entity):
        beg = Biodatabase.objects.get(name=entity.biodatabase.name.replace("_prots", ""))
        f = Seqfeature.objects.seqfeature_from_locus_tag(beg.biodatabase_id, entity.accession)
        f = list(f)[0]
        return ", ".join([f.qualifiers.get(term__name=x).value
                          for x in ["gene_symbol", "old_locus_tag", "protein_id", "Alias", "gene"]
                          if f.qualifiers.filter(term__name=x).exists()]),

    columns = [TableAttribute("Description", "feature", lambda e: description(e),
                              query=lambda param: Q(value__icontains=param)),
               TableAttribute("Accession", "entry", lambda e: e.accession,
                              query=lambda param: Q(accession__contains=param),
                              url=lambda e: e.get_absolute_url()),
               TableAttribute("Location", "feature", lambda e: location(e),
                              query=None),
               TableAttribute("Aliases", "feature", lambda e: ", ".join(genes(e)),
                              query=lambda param: Q(value__icontains=param)),

               ]

    queryset = Bioentry.objects.prefetch_related("dbxrefs__dbxref", "qualifiers__term").filter(
        biodatabase_id=assembly_id)

    page = page_from_request(queryset.count(), request)
    queryset = queryset[page.offset:page.offset + page.size]

    return render(request, 'biosql/protein_list.html', {
        "objects": enumerate(queryset, page.offset + 1), "filters": filters, "columns": columns, "page_obj": page,
        "query": "-",
        "params": request.GET, "sidebarleft": 1})


from django.db.models import Q


def assembly_view(request, pk):
    be = Biodatabase.objects.get(biodatabase_id=pk)
    if "_prots" in be.name:
        be2 = (Biodatabase.objects.prefetch_related("entries__features__type_term"))
        be = be2.get(name=be.name.replace("_prots", ""))  # TODO meter en otro lado esta regla
    else:
        be = Biodatabase.objects.get(biodatabase_id=pk)

    page = None
    if be.entries.count() < 500:
        # ontology = Ontology.objects.get(name="Stats")
        # qvs = BioentryQualifierValue.objects.prefetch_related("bioentry","term").filter(
        #     bioentry__biodatabase=be, term__ontology=ontology)
        # qualifiers = defaultdict(lambda: {})
        # for qv in qvs:
        #     qualifiers[qv.bioentry.accession][qv.term.identifier.split("_")[-1]] = qv.value
        # # x.accession: {y.term.identifier.split("_")[-1]: y.value for y in x.qualifiers.all() }
        be = Biodatabase.objects.prefetch_related("entries__qualifiers__terms").get(biodatabase_id=pk)
        queryset = be.entries
    else:
        queryset = be.entries.prefetch_related("qualifiers__term")
        page = page_from_request(queryset.count(), request)
        queryset = queryset[page.offset:page.offset + page.size]

    tax = be.entries.first().taxon
    assembly = {
        "name": be.name, "description": be.description, "ncbi_tax": tax, "biodatabase_id": be.biodatabase_id
    }
    lengths = {}
    for x in queryset.all():
        seq = Biosequence.objects.raw("""
        SELECT bioentry_id, version , length , alphabet ,"" seq
        FROM biosequence WHERE bioentry_id = %i ;
        """ % (x.bioentry_id))[0]
        lengths[x.accession] = seq.length

    return render(request, 'biosql/assembly_detail.html', {"lengths": lengths,
                                                           "object": assembly,
                                                           "contigs": queryset, "page_obj": page, "query": None,
                                                           "sidebarleft": {}})
