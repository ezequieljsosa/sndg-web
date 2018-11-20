from django.shortcuts import render
from django.views.generic import TemplateView, DetailView

from django.db.models import Prefetch

from .models import Taxon, Biosequence, Bioentry, Seqfeature, Biodatabase, Term, Ontology,SeqfeatureQualifierValue
from bioresources.models import Assembly, Publication,Sample
from django.db.models import Q


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

    ont_sfk = Ontology.objects.get(name="SeqFeature Keys")
    ont_at = Ontology.objects.get(name="Annotation Tags")

    term_cds = Term.objects.get(ontology=ont_sfk, name="CDS")
    term_lt = Term.objects.get(ontology=ont_at, name="locus_tag")

    if be.biodatabase.name.endswith("prots"):
        beg = Biodatabase.objects.get(name=be.biodatabase.name.replace("_prots", ""))
        taxon = beg.entries.first().taxon
        feature = Seqfeature.objects.raw("""
        SELECT sf.seqfeature_id, sf.bioentry_id, sf.type_term_id, sf.source_term_id , sf.display_name,  sf.rank
        FROM biodatabase bdb,bioentry be, seqfeature sf, seqfeature_qualifier_value sfqv
        WHERE bdb.biodatabase_id = %i AND be.biodatabase_id = bdb.biodatabase_id AND 
          sf.bioentry_id = be.bioentry_id AND sf.type_term_id = %i AND 
          sf.seqfeature_id = sfqv.seqfeature_id AND sfqv.term_id = %i AND 
          sfqv.value = "%s"
        """ % (beg.biodatabase_id, term_cds.term_id, term_lt.term_id, be.accession))
        feature = list(feature)[0]

        locations = list(feature.locations.all())
        start = locations[0].start_pos
        end = locations[-1].end_pos

        seq = Biosequence.objects.raw("""
        SELECT bioentry_id, version , length , alphabet ,SUBSTRING( seq,%i,%i ) seq
        FROM biosequence WHERE bioentry_id = %i ;
        """ % (start, end - start, feature.bioentry_id))[0]
        functions = {"biological_process": [], "molecular_function": [], "cellular_component": []}
        for qual in be.qualifiers.filter(term__dbxrefs__dbxref__dbname="go",term__dbxrefs__dbxref__accession="goslim_generic"):
            for dbxref in qual.term.dbxrefs.all():
                if dbxref.dbxref.accession in functions:
                    functions[dbxref.dbxref.accession].append(qual.term)

    else:
        beg = be.biodatabase

    graph = entry_graph(be, beg)

    if be.biodatabase.name.endswith("prots"):
        return render(request, 'biosql/protein_detail.html', {
            "functions": functions, "graph": graph,"accession":be.biodatabase.name.replace("_prots", ""),
            "object": be, "feature": feature, "taxon": taxon, "seq": seq, "start": start, "end": end,
            "sidebarleft": 1})
    else:
        return render(request, 'biosql/biosequence_detail.html', {
            "object": be, "graph": graph,
            "sidebarleft": 0})  # "sidebarrigth": {"news": [{"title": "n1", "text": "lalala"}]


def publications_from_resource_graph(resource, graph, external_orgs, resource_identifier=None):
    if not resource_identifier:
        resource_identifier = resource.name

    for x in resource.sources.filter(source__type="pmc").all():

        graph["nodes"].append({"id": x.source.name, "label": labelize(x.source.name), "color": "SlateBlue"})
        graph["edges"].append({"from": resource_identifier, "to": x.source.name})

        p = Publication.objects.get(id=x.source.id)
        for i, affiliation in enumerate(p.affiliations.filter(organizations__country="Argentina").all()):
            if (affiliation.author.complete_name not in
                                                       [x["id"] for x in graph["nodes"]]):
                graph["nodes"].append({"id": affiliation.author.complete_name,
                                       "label": affiliation.author.complete_name() + " #" + str(i),
                                       "color": "cyan" })

                for org in affiliation.organizations.all():
                    if org.name not in [y["id"] for y in graph["nodes"]]:
                        graph["nodes"].append({"id": org.name,"type":"organization",
                                               "label": labelize(org.name),
                                               "color": "green" if org.country == "Argentina" else "purple"})
                    graph["edges"].append({"from": x.source.name, "to": org.name})
                    graph["edges"].append({"from": org.name, "to": affiliation.author.complete_name})
        for i, affiliation in enumerate(p.affiliations.filter(~Q(organizations__country="Argentina")).all()):
            for org in affiliation.organizations.all():
                if org not in external_orgs:
                    external_orgs.append(org)


from django.db.models import Q


def assembly_graph(assembly):
    external_orgs = []
    assembly_name = assembly.external_ids.filter(type="accession").first().identifier
    graph = {"nodes": [{"id": assembly_name, "label": assembly_name,"type":"assembly", "color": "orange"}],
             "edges": []}

    publications_from_resource_graph(assembly, graph, external_orgs, assembly_name)
    for rel in assembly.sources.filter(source__type__in = ["publication","sample"]).all(): #.filter(~Q(source__type="publication"))
        if rel.source.name not in graph["nodes"]:

            graph["edges"].append({"from": assembly_name, "to": rel.source.name})
            if rel.source.type == "publication":
                graph["nodes"].append({"id": rel.source.name, "type":"publication","label": labelize(rel.source.name), "color": "SlateBlue"})
                publications_from_resource_graph(rel.source, graph, external_orgs)

            else:
                graph["nodes"].append({"id": rel.source.name,"type":"sample", "label": labelize(rel.source.name), "color": "Grey"})
                sample = Sample.objects.get(id=rel.source.id)
                sample_origin = sample.origin_dict()
                for k,v in sample_origin.items():
                    name = k + ": " + v
                    if name not in graph["nodes"]:
                        graph["nodes"].append({"id":name , "label": labelize(name), "color": "AliceBlue"})
                        graph["edges"].append({"from": name, "to": sample.name})

    nodes = []
    added = []
    for x in graph["nodes"]:
        if x["id"] not in added:
            added.append(x["id"])
            nodes.append(x)
    graph["nodes"] = nodes
    graph["external_orgs"] = external_orgs
    return graph


def entry_graph(be, beg):
    assembly = Assembly.objects.get(external_ids__identifier=beg.name)
    graph = assembly_graph(assembly)
    graph["nodes"].append({"id": be.accession, "label": be.accession,"type":"seq", "color": "red"})
    graph["edges"].append({"from": be.accession, "to": beg.name})
    return graph





def assembly_view(request, pk):
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
