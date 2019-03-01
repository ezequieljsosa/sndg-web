from django.core.management.base import CommandError
from django.db import transaction
from django.db.utils import IntegrityError

from django_tqdm import BaseCommand
from django.db.models import Q

import traceback

import datetime
import json
import os
import pandas as pd

from vardb.models import Protocol, AntibioticResistance,ReportedAllele,Allele,Effect
from biosql.models import Term, Ontology,Biodatabase,Bioentry

from tqdm import tqdm


class Command(BaseCommand):
    help = 'Loads the obo files to the database. Prepared for GO and SO'

    def __init__(self, stdout=None, stderr=None, no_color=False):
        super().__init__(stdout=stdout, stderr=stderr, no_color=no_color)
        self.resist = None

    def add_arguments(self, parser):
        pass

    def handle(self, *args, **options):
        """https://www.ncbi.nlm.nih.gov/biosample/docs/antibiogram-myco/
        https://www.ncbi.nlm.nih.gov/biosample/docs/antibiogram/
         https://www.ncbi.nlm.nih.gov/biosample/docs/beta-lactamase/
        """
        # Ontology Assay Result
        if not Ontology.objects.filter(name="Assay Results").exists():
            aro = Ontology(name="Assay Results")
            aro.save()
            Term(name="Positive", identifier="Positive", ontology=aro).save()
            Term(name="Negative", identifier="Negative", ontology=aro).save()
            Term(name="Inconclusive", identifier="Inconclusive", ontology=aro).save()

        # Ontology Genotype Support Status
        if not Ontology.objects.filter(name="Genotype Support Status").exists():
            gsso = Ontology(name="Genotype Support Status")
            gsso.save()
            Term(name="Conclusive", identifier="Conclusive", ontology=gsso).save()
            Term(name="Possible", identifier="Possible", ontology=gsso).save()
            Term(name="Hint", identifier="Hint", ontology=gsso).save()

        # Ontology Antibiotic
        # Protocol : unknown...
        if not Ontology.objects.filter(name="Antibiotics").exists():
            ao = Ontology(name="Antibiotics")
            ao.save()

            for l in """streptomycin\tSTR
                        isoniazid\tINH
                        rifampicin\tRIF
                        ethambutol\tEMB
                        kanamycin\tKAN
                        amikacin\tAMK
                        capreomycin\tCAP
                        ethionamide\tETH
                        fluoroquinolone\tFLQ
                        aminoglycoside\tAMI
                        viomycin\tVIO
                        linezolid\tZD                  
                        para-aminosalicylic acid\tPAS""".split("\n"):
                name, ident = [x.strip() for x in l.split("\t")][:2]

                t = Term(name=name, identifier=ident, ontology=ao)
                t.save()
                p = AntibioticResistance(name="%s resistant" % name, antibiotic=t).save()
                Protocol(name="%s antibiogram" % name, phenotype=p).save()

        path_db = "/home/eze/Downloads/andytb - andytb.csv"
        self.resist = pd.read_csv(path_db)
        self.resist["AApos"] = [int(x) if x != "-" else "" for x in self.resist.AApos]
        self.resist["LocusTag"] = [x.split("_")[0] for x in self.resist.GeneID]
        self.resist["NucleotidePosH37"] = [ int(str(x).split("/")[0]) if str(x) not in ["-","nan"] else "" for x in
                                           self.resist.NucleotidePosH37]


        # {'AMINOGLYCOSIDES', 'ETHAMBUTOL', 'CAPREOMYCIN', 'STREPTOMYCIN', 'ISONIAZID', 'PYRAZINAMIDE', 'FLUOROQUINOLONES', 'FOSFOMYCIN', 'ETHIONAMIDE', 'KANAMYCIN', 'AMIKACIN', 'LINEZOLID', 'PARA-AMINOSALISYLIC_ACID', 'CLOFAZIMINE', 'BEDAQUILINE', 'VIOMYCIN', 'RIFAMPICIN', '-'}

        bdb = Biodatabase.objects.get(name="GCF_000195955.2")
        seq = Bioentry.objects.get(biodatabase=bdb, identifier="NC_000962.3")
        # for _,r in tqdm(self.resist.iterrows(),total=len(self.resist)):
        #     pheno = AntibioticResistance.objects.get(name=map_tb_anti[""])
        #
        #     ra = ReportedAllele( phenotype = pheno,   reported_in = r.Source)
        #     ra.save()
        #
        #     if r.NucleotidePosH37:
        #         allele = Allele.objects.filter(variant_fk__pos=r.NucleotidePosH37,variant_fk__contig=seq,
        #                                     alt=r.ALT)
        #         if allele.exists():
        #             allele = Allele.objects.get(variant_fk__pos=r.NucleotidePosH37,variant_fk__contig=seq,
        #                                            alt=r.ALT)
        #         else:
        #             allele = Allele.objects.get(variant_fk__pos=r.NucleotidePosH37,variant_fk__contig=seq,
        #                                         alt=r.ALT)
        #
        #     effect =
        #     variant_query = Variant.objects.filter(
        #         pos=var.POS, contig=seq, ref=var.REF)
        #
        #     if not variant_query.exists():
        #         new_variant = Variant(contig=seq, pos=var.POS, ref=var.REF)
        #         new_variant.save()
        #     else:
        #         new_variant = variant_query.get()



