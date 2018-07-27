from django.core.management.base import CommandError
from django.db import transaction
from django.db.utils import IntegrityError

from django_tqdm import BaseCommand

import traceback

import datetime
import json
import os
from goatools.obo_parser import GODag
from biosql.models import Taxon, TaxonName, TaxIdx, Term, TermRelationship, Ontology, TermIdx


class Command(BaseCommand):
    help = 'Loads the obo files to the database. Prepared for GO and SO'

    def __init__(self, stdout=None, stderr=None, no_color=False):
        super().__init__(stdout=stdout, stderr=stderr, no_color=no_color)
        self.is_a = Term.objects.get(identifier="is_a")

    def travel_tax_child(self, parent, iter_adv, txt_acc, buffer):

        txt = [x.name for x in parent.names.all()] + [str(parent.ncbi_taxon_id)]
        txt += txt_acc
        for t in Taxon.objects.prefetch_related("names", "children").filter(parent_taxon=parent):
            if t.ncbi_taxon_id != 1:
                self.travel_tax_child(t, iter_adv, txt,buffer)

        if parent.ncbi_taxon_id != 1:
            buffer["data"].append(TaxIdx(tax=parent, text=txt))
        if len(buffer["data"]) > 5000:
            TaxIdx.objects.bulk_create(buffer["data"])
            buffer["data"] = []

        iter_adv.update(1)

    def travel_go_child(self, cache, parent, iter_adv, text_acc, buffer):
        if parent.term_id in cache:
            return
        txt = " ".join([parent.name, parent.identifier, parent.definition] +
                       [x.term.name for x in parent.dbxrefs.all()])  # [x.synonym for x in parent.synonyms.all()]
        txt += text_acc
        for c in TermRelationship.objects.prefetch_related("object_term__dbxrefs__term",
                                                           ).filter(subject_term=parent,
                                                                    predicate_term=self.is_a):
            if c.object_term.identifier not in cache:
                self.travel_go_child(cache, c.object_term, iter_adv, txt,buffer)

        if parent.term_id in cache:
            return
        cache[parent.term_id] = 1
        buffer["data"].append(TermIdx(term=parent, text=txt))
        iter_adv.update(1)
        if len(buffer["data"]) == 5000:
            with transaction.atomic():
                TermIdx.objects.bulk_create(buffer["data"])
            buffer["data"] = []

    def add_arguments(self, parser):
        pass

    def handle(self, *args, **options):
        # root = Taxon.objects.prefetch_related("names").get(ncbi_taxon_id=1)
        # count = Taxon.objects.count()
        # iter_adv = self.tqdm(total=count)
        #
        # buffer = {"data":[]}
        # for c in root.children.all():
        #     self.travel_tax_child(c, iter_adv, "",buffer)
        # bulk_create(buffer["data"])


        buffer = {"data":[]}
        cache = {x.term_id: 1 for x in self.tqdm(TermIdx.objects.select_related("term").all())}
        go = Ontology.objects.get(name="Gene Ontology")
        count = Term.objects.filter(ontology=go).count()

        # Term.objects.prefetch_related("synonyms").filter(ontology=go).first() "synonyms",
        iter_adv = self.tqdm(total=count)
        iter_adv.update( len(cache))
        gos = list(Term.objects.filter(ontology=go))
        for c in gos:
            iter_adv.set_description(c.identifier)
            self.travel_go_child(cache, c, iter_adv, "",buffer)
