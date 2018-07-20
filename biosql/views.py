from django.shortcuts import render
from django.views.generic import TemplateView, DetailView

from django.db.models import Prefetch

from .models import Taxon, Biosequence, Bioentry, Seqfeature, Biodatabase
from bioresources.models import Assembly, Publication


class AboutView(TemplateView):
    template_name = "about.html"

    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)
        context['jojo'] = "sdalmkghdaf"
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


def sequence_view(request, pk):
    # pf_features = Prefetch()
    be = (Bioentry.objects.select_related("biodatabase").select_related("taxon")
        .prefetch_related("dbxrefs__dbxref", "qualifiers__term", "seq",
                          "features__locations", "features__source_term",
                          "features__type_term", "features__qualifiers"))
    be = be.get(bioentry_id=pk)

    if be.biodatabase.name.endswith("prots"):
        beg = Biodatabase.objects.get(name=be.biodatabase.name.replace("_prots", ""))
        taxon = beg.entries.first().taxon
        feature = Seqfeature.objects.raw("""
        SELECT sf.seqfeature_id, sf.bioentry_id, sf.type_term_id, sf.source_term_id , sf.display_name,  sf.rank
        FROM biodatabase bdb,bioentry be, seqfeature sf, seqfeature_qualifier_value sfqv
        WHERE bdb.biodatabase_id = %i AND be.biodatabase_id = bdb.biodatabase_id AND 
          sf.bioentry_id = be.bioentry_id AND sf.type_term_id = 20 AND 
          sf.seqfeature_id = sfqv.seqfeature_id AND sfqv.term_id = 18 AND 
          sfqv.value = "%s"
        """ % (beg.biodatabase_id, be.accession))
        feature = list(feature)[0]

        locations = list(feature.locations.all())
        start = locations[0].start_pos
        end = locations[-1].end_pos

        seq = Biosequence.objects.raw("""
        SELECT bioentry_id, version , length , alphabet ,SUBSTRING( seq,%i,%i ) seq
        FROM biosequence WHERE bioentry_id = %i ;
        """ % (start, end - start, feature.bioentry_id))[0]
        functions = {"biological_process": [], "molecular_function": [], "cellular_component": []}
        for qual in be.qualifiers.all():
            for dbxref in qual.term.dbxrefs.all():
                if dbxref.dbxref.accession in functions:
                    functions[dbxref.dbxref.accession].append(qual.term)

    else:
        beg = be.biodatabase

    graph = entry_graph(be, beg)

    if be.biodatabase.name.endswith("prots"):
        return render(request, 'biosql/protein_detail.html', {
            "functions": functions, "graph": graph,
            "object": be, "feature": feature, "taxon": taxon, "seq": seq, "start": start, "end": end,
            "sidebarleft": 1})
    else:
        return render(request, 'biosql/biosequence_detail.html', {
            "object": be,"graph": graph,
            "sidebarleft": 0 }) # "sidebarrigth": {"news": [{"title": "n1", "text": "lalala"}]


def assembly_graph(assembly):
    external_orgs = []
    assembly_name = assembly.external_ids.filter(type="accession").first().identifier
    graph = {"nodes": [ {"id": assembly_name, "label": assembly_name, "color": "orange"}],
             "edges": []}
    for rel in assembly.sources.all():
        if rel.source.name not in graph["nodes"]:
            graph["nodes"].append({"id": rel.source.name, "label": labelize(rel.source.name), "color": "SlateBlue"})
            graph["edges"].append({"from": assembly_name, "to": rel.source.name})
            for x in rel.source.sources.all():
                if x.source.type == "pmc":
                    p = Publication.objects.get(id=x.source.id)
                    for i, affiliation in enumerate(p.affiliations.all()):
                        if affiliation.author.arg_affiliation:
                            graph["nodes"].append({"id": affiliation.author.complete_name,
                                                   "label": affiliation.author.complete_name() + " #" + str(i),
                                                   "color": "cyan" if affiliation.author.arg_affiliation else "grey"})

                            for org in affiliation.organizations.all():
                                if org.name not in [y["id"] for y in graph["nodes"]]:
                                    graph["nodes"].append({"id": org.name,
                                                           "label": labelize(org.name),
                                                           "color": "green" if org.country == "Argentina" else "purple"})
                                graph["edges"].append({"from": rel.source.name, "to": org.name})
                                graph["edges"].append({"from": org.name, "to": affiliation.author.complete_name})
                        else:
                            for org in affiliation.organizations.all():
                                if org not in external_orgs:
                                    external_orgs.append(org)
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
    graph["nodes"].append({"id": be.accession, "label": be.accession, "color": "red"})
    graph["edges"].append({"from": be.accession, "to": beg.name})
    return graph

def assembly_view(request, pk):
    # pf_features = Prefetch()
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

    return render(request, 'biosql/assembly_detail.html', { "lengths": lengths,
            "object": assembly,"graph": graph, "contigs":be.entries.all(),
            "sidebarleft": {} })
