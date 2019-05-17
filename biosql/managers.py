"""

"""
import datetime
from django.db import models
from django.db.models.query import QuerySet
from django.db.models import Prefetch, Q


class SeqfeatureQuerySet(QuerySet):

    def geneproducts(self, index_updated=False):
        from .models import Seqfeature
        return Seqfeature.objects.prefetch_related("bioentry", "dbxrefs__dbxref", "qualifiers__term",
                                                   "locations").filter(index_updated=index_updated,
                                                                       bioentry__biodatabase_id=1,
                                                                       bioentry_id = 3724603,
                                                                       type_term__name__in=Seqfeature.ENTRY_TYPES)

        # def iter_qs(qs):
        #     for idx, f in enumerate(qs, page_offset + 1):
        #         lts = f.qualifiers.get(term__name="locus_tag").value
        #         product = f.qualifiers.filter(term__name="product")
        #         product = product.first().value if product.exists() else "-"
        #
        #         # feature = Seqfeature.objects.seqfeature_from_locus_tag(biodbid, lts)
        #         # feature = list(feature)[0]
        #
        #         yield {
        #             "idx": idx,
        #             "ftype": f.type_term.name,
        #             "product": product,
        #             "lt": lts,
        #             "genes": ", ".join([f.qualifiers.get(term__name=x).value
        #                                 for x in ["gene_symbol", "old_locus_tag", "protein_id", "Alias", "gene"]
        #                                 if f.qualifiers.filter(term__name=x).exists()]),
        #             "location": f.bioentry.accession + ":" + ";".join(
        #                 [str(l) for l in f.locations.all()])
        #         }

    def seqfeature_from_locus_tag(self, biodatabase_id, accession):
        from .models import Ontology, Term, Seqfeature
        ont_sfk = Ontology.objects.get(name="SeqFeature Keys")
        ont_at = Ontology.objects.get(name="Annotation Tags")

        term_cds = Term.objects.get(ontology=ont_sfk, name="mRNA")
        term_lt = Term.objects.get(ontology=ont_at, name="locus_tag")

        return Seqfeature.objects.raw("""
        SELECT sf.seqfeature_id, sf.bioentry_id, sf.type_term_id, sf.source_term_id , sf.display_name,  sf.rank
        FROM biodatabase bdb,bioentry be, seqfeature sf, seqfeature_qualifier_value sfqv
        WHERE bdb.biodatabase_id = %i AND be.biodatabase_id = bdb.biodatabase_id AND 
          sf.bioentry_id = be.bioentry_id AND sf.type_term_id = %i AND 
          sf.seqfeature_id = sfqv.seqfeature_id AND sfqv.term_id = %i AND 
          sfqv.value = "%s"
        """ % (biodatabase_id, term_cds.term_id, term_lt.term_id, accession))


class SeqfeatureManager(models.Manager):

    def get_query_set(self):
        return SeqfeatureQuerySet(self.model)

    def __getattr__(self, attr, *args):
        # see https://code.djangoproject.com/ticket/15062 for details
        if attr.startswith("_"):
            raise AttributeError
        return getattr(self.get_query_set(), attr, *args)


class BioentryQuerySet(QuerySet):

    def proteins(self, index_updated=False):
        from .models import Bioentry
        return Bioentry.objects.prefetch_related("dbxrefs__dbxref", "qualifiers__term__dbxrefs__dbxref", "seq").filter(
            index_updated=index_updated,
            bioentry_id = 3724603,
            biodatabase_id=388)


class BioentryManager(models.Manager):

    def get_query_set(self):
        return BioentryQuerySet(self.model)

    def __getattr__(self, attr, *args):
        # see https://code.djangoproject.com/ticket/15062 for details
        if attr.startswith("_"):
            raise AttributeError
        return getattr(self.get_query_set(), attr, *args)
