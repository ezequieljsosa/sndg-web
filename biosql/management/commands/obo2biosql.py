from django.core.management.base import BaseCommand, CommandError
from django.db import transaction
from django.db.utils import IntegrityError

import traceback

import datetime
import json
import os
from tqdm import tqdm

from goatools.obo_parser import GODag
from biosql.models import Ontology, Term, TermRelationship, Dbxref, TermDbxref, TermSynonym


class Command(BaseCommand):
    help = 'Loads the obo files to the database. Prepared for GO and SO'

    def __init__(self, stdout=None, stderr=None, no_color=False):
        super().__init__(stdout=stdout, stderr=stderr, no_color=no_color)
        self.ontology = None
        self.dbmap = {}
        self.is_a = None
        self.relmap = {}
        self.ontology_name = ""
        self.obo_path = None

    def add_arguments(self, parser):
        parser.add_argument('--obo_path', required=True)
        parser.add_argument('--ontology', required=True, choices=["go", "so"])

    def create_base_terms(self):
        sfs_ontology = Ontology.objects.get_or_create(name="SeqFeature Sources")[0]
        Term.objects.get_or_create(identifier="manual", name="manual", version=1, ontology=sfs_ontology,
                                   definition="added or corrected by a person")
        Term.objects.get_or_create(identifier="bibliography", name="bibliography", version=1, ontology=sfs_ontology,
                                   definition="found in bibliography")
        Term.objects.get_or_create(identifier="experimental", name="experimental", version=1, ontology=sfs_ontology,
                                   definition="the annotation was obtained experimentally")
        Term.objects.get_or_create(identifier="calculated", name="calculated", version=1, ontology=sfs_ontology,
                                   definition="the annotation was obtained using software")
        Term.objects.get_or_create(identifier="other", name="other", version=1, ontology=sfs_ontology,
                                   definition="")

        self.is_a = Term.objects.get(identifier="is_a", name="is_a", version=1, ontology=graph_ontology, definition="")

        self.ontology = Ontology.objects.get_or_create(name=self.ontology_name)

        graph_ontology = Ontology.objects.get_or_create(name="Graph", definition="")[0]
        self.dbmap = {
            "biological_process": Dbxref.objects.get_or_create(dbname="go", accession="biological_process", version=1)[
                0],
            "molecular_function": Dbxref.objects.get_or_create(dbname="go", accession="molecular_function", version=1)[
                0],
            "cellular_component": Dbxref.objects.get_or_create(dbname="go", accession="cellular_component", version=1)[
                0],
            "biosapiens": Dbxref.objects.get_or_create(dbname="so", accession="biosapiens", version=1)[0],
            "DBVAR": Dbxref.objects.get_or_create(dbname="so", accession="DBVAR", version=1)[0],
            "SOFA": Dbxref.objects.get_or_create(dbname="so", accession="SOFA", version=1)[0],
        }

        for x in """goslim_agr "AGR slim"
        goslim_aspergillus "Aspergillus GO slim"
        goslim_candida "Candida GO slim"
        goslim_chembl "ChEMBL protein targets summary"
        goslim_generic "Generic GO slim"
        goslim_goa "GOA and proteome slim"
        goslim_metagenomics "Metagenomics GO slim"
        goslim_mouse "Mouse GO slim"
        goslim_pir "PIR GO slim"
        goslim_plant "Plant GO slim"
        goslim_pombe "Fission yeast GO slim"
        goslim_synapse "synapse GO slim"
        goslim_virus "Viral GO slim"
        goslim_yeast "Yeast GO slim"
        gosubset_prok "Prokaryotic GO subset"
        virus_checked "Viral overhaul terms" """.split("\n"):
            k = x.strip().split(' ')[0]
            self.dbmap[k] = Dbxref.objects.get_or_create(dbname="go", accession=k, version=1)[0]

        part_of = Term.objects.get(identifier="part_of", )
        regulates = Term.objects.get(identifier="regulates")
        negatively_regulates = Term.objects.get(identifier="negatively_regulates")
        positively_regulates = Term.objects.get(identifier="positively_regulates")

        has_quality = Term.objects.get_or_create(identifier="has_quality", name="has_quality", version=1, ontology=
        graph_ontology)[0]
        derives_from = Term.objects.get_or_create(identifier="derives_from", name="derives_from", version=1, ontology=
        graph_ontology)[0]
        has_origin = Term.objects.get_or_create(identifier="has_origin", name="has_origin", version=1, ontology=
        graph_ontology)[0]
        has_part = Term.objects.get_or_create(identifier="has_part", name="has_part", version=1, ontology=
        graph_ontology)[0]
        transcribed_to = \
            Term.objects.get_or_create(identifier="transcribed_to", name="transcribed_to", version=1, ontology=
            graph_ontology)[0]

        variant_of = Term.objects.get_or_create(identifier="variant_of", name="variant_of", version=1, ontology=
        graph_ontology)[0]

        transcribed_from = \
            Term.objects.get_or_create(identifier="transcribed_from", name="transcribed_from", version=1, ontology=
            graph_ontology)[0]

        adjacent_to = Term.objects.get_or_create(identifier="adjacent_to", name="adjacent_to", version=1, ontology=
        graph_ontology)[0]

        member_of = Term.objects.get_or_create(identifier="member_of", name="member_of", version=1, ontology=
        graph_ontology)[0]

        contains = Term.objects.get_or_create(identifier="contains", name="contains", version=1, ontology=
        graph_ontology)[0]
        non_functional_homolog_of = Term.objects.get_or_create(identifier="non_functional_homolog_of",
                                                               name="non_functional_homolog_of", version=1, ontology=
                                                               graph_ontology)[0]

        overlaps = Term.objects.get_or_create(identifier="overlaps", name="overlaps", version=1, ontology=
        graph_ontology)[0]

        guided_by = Term.objects.get_or_create(identifier="guided_by", name="guided_by", version=1, ontology=
        graph_ontology)[0]

        self.relmap = {
            'negatively_regulates': negatively_regulates, 'regulates': regulates,
            'positively_regulates': positively_regulates, 'part_of': part_of,
            "has_quality": has_quality, 'derives_from': derives_from, 'has_origin': has_origin,
            "has_part": has_part, "transcribed_to": transcribed_to, "variant_of": variant_of,
            "transcribed_from": transcribed_from, "adjacent_to": adjacent_to,
            "member_of": member_of, "contains": contains, "non_functional_homolog_of": non_functional_homolog_of,
            "overlaps": overlaps, "guided_by": guided_by
        }

    def create_terms(self):
        ids = []
        with open(self.obo_path) as h:  # filter the ids because GODag iterates over alt_ids too
            for l in h.readlines():
                if l.startswith("id: GO:"):
                    ids.append(l.split(" ")[1].strip())

        go_dag = GODag(self.obo_path, load_obsolete=True,
                       optional_attrs=["def", "synonym", "subset", "alt_id", "dbxref"])

        for go in tqdm(ids):
            term = go_dag[go]
            with transaction.atomic():
                dbTerm = Term(name=term.name,
                              definition=term.defn,
                              identifier=go,
                              is_obsolete="T" if term.is_obsolete else "F",
                              ontology=ontology)
                dbTerm.save()
                # intersetc = (set([ x.dbxref.accession for x in  dbTerm.dbxrefs.all()]) &
                #              set(["biological_process","molecular_function","cellular_component"]))
                # if not intersetc:
                #
                #     termdbref = TermDbxref(term=t, dbxref=self.dbmap[term.namespace], rank=1)
                #     termdbref.save()

                termdbref = TermDbxref(term=dbTerm, dbxref=self.dbmap[term.namespace], rank=1)
                termdbref.save()

                for subset in term.subset:
                    if subset in dbmap:
                        termdbref = TermDbxref(term=dbTerm, dbxref=self.dbmap[subset], rank=1)
                        termdbref.save()
                for synonym in term.synonym:
                    syn = TermSynonym(term=dbTerm, synonym=synonym[0][:255])
                    syn.save()
                for synonym in term.alt_ids:
                    syn = TermSynonym(term=dbTerm, synonym=synonym[0][:255])
                    syn.save()

                self.cache[go] = dbTerm

    def create_relationships(self):
        go_dag = GODag(self.obo_path, optional_attrs=['relationship'], load_obsolete=False)

        terms = Term.objects.filter(ontology=self.ontology)

        for dbTerm in tqdm(terms.all(), total=terms.count()):
            go = dbTerm.identifier
            term = go_dag[dbTerm.identifier]
            with transaction.atomic():
                for child in term.children:
                    if go in cache and child.id in cache:
                        r = TermRelationship(
                            subject_term=cache[go],  # parent
                            predicate_term=is_a,
                            object_term=cache[child.id],  # child
                            ontology=ontology
                        )
                        r.save()

                for rel, terms in term.relationship.items():
                    for child in terms:
                        if go in cache and child.id in cache:
                            r = TermRelationship(
                                subject_term=cache[go],  # parent
                                predicate_term=relmap[rel],
                                object_term=cache[child.id],  # child
                                ontology=ontology
                            )
                            r.save()

    def handle(self, *args, **options):

        if options["ontology"] == "go":
            self.ontology_name = "Gene Ontology"
        elif options["ontology"] == "so":
            self.ontology_name = "Sequence Ontology"

        self.obo_path = options["obo_path"]

        self.create_base_terms()

        self.create_terms()
        self.create_relationships()
