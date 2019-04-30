import datetime
import json
import os
import sys
from django.db.utils import IntegrityError
from django.core.management.base import BaseCommand, CommandError
from tqdm import tqdm
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import Bio.SeqIO as bpio
import pandas as pd
from chembl_model.models import TargetDictionary, ZincCompound,ZincProperty
from itertools import islice

class Command(BaseCommand):
    help = 'Loads gene ontology terms'

    def add_arguments(self, parser):
        parser.add_argument('command')

    def handle(self, *args, **options):
        if options["command"] == "fasta":
            """SELECT t.chembl_id AS target_chembl_id,
            t.pref_name        AS target_name,
            t.target_type,
            c.accession        AS protein_accession,
            c.sequence         AS protein_sequence
            FROM target_dictionary t
              JOIN target_type tt ON t.target_type = tt.target_type
              JOIN target_components tc ON t.tid = tc.tid
              JOIN component_sequences c ON tc.component_id = c.component_id
            AND tt.parent_type  = 'PROTEIN';
            """
            accessions = {}
            with sys.stdout as h:
                qs = TargetDictionary.objects.using("chembl_25").prefetch_related("components__component").filter(
                    target_type__parent_type="PROTEIN")
                total = qs.count()
                for target_dictionary in tqdm(qs, total=total, file=sys.stderr):
                    for target_component in target_dictionary.components.all():
                        rid = target_component.component.accession
                        if rid not in accessions:
                            accessions[rid] = 1
                            r = SeqRecord(id=rid, name="", description=target_dictionary.pref_name,
                                          seq=Seq(target_component.component.sequence))
                            bpio.write(r, h, "fasta")
        if options["command"] == "load_zinc":
            ZincCompound.objects.all().delete()
            ZincProperty.objects.all().delete()
            df = pd.read_table("/data/databases/zinc/6_prop.xls")
            bulk = []
            for i, r in tqdm(df.iterrows()):
                d = r.to_dict()
                d["code"] = d["ZINC_ID"]
                del d["ZINC_ID"]
                bulk.append(Zinc(**d))

                if i and (i % 1000 == 0):
                    Zinc.objects.bulk_create(bulk)
                    bulk = []
