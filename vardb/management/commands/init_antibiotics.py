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

from vardb.models import Protocol, AntibioticResistance,ReportedAllele,Allele,Effect,GenotypeSupport
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
        if not Ontology.objects.filter(name=GenotypeSupport.STATUS_ONTOLOGY).exists():
            gsso = Ontology(name=GenotypeSupport.STATUS_ONTOLOGY)
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
                        clofazimine\tCLO
                        bedaquilin\t
                        viomycin\tVIO
                        linezolid\tZD                  
                        para-aminosalicylic_acid\tPAS""".split("\n"):
                name, ident = [x.strip() for x in l.split("\t")][:2]

                t = Term(name=name, identifier=ident, ontology=ao)
                t.save()
                p = AntibioticResistance(name="%s resistant" % name, antibiotic=t).save()
                Protocol(name="%s antibiogram" % name, phenotype=p).save()








