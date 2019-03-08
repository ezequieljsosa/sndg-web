from django.core.management.base import CommandError
from django.db import transaction
from django.db.utils import IntegrityError

from django_tqdm import BaseCommand
from django.db.models import Q

import traceback

import datetime
import json
import os
import logging
from biosql.models import Bioentry, Biodatabase, Seqfeature
from vardb.models import Variant, Allele, Variantannotation, Variantassignment, Variantcollection, Effect,AlleleEffect, AntibioticResistance
import vcf
import subprocess
from tqdm import tqdm
import hgvs.parser

_log = logging.getLogger(__name__)

# log = logging.getLogger('django.db.backends')
# log.setLevel(logging.DEBUG)
# log.addHandler(logging.StreamHandler())

class Command(BaseCommand):
    help = 'Updates a variant collection genotypes'

    def __init__(self, stdout=None, stderr=None, no_color=False):
        super().__init__(stdout=stdout, stderr=stderr, no_color=no_color)




    def add_arguments(self, parser):
        parser.add_argument('-vc','--variant_collection', help="variant collection sample name", required=True)
        # parser.add_argument('--reference', required=True)

    def handle(self, *args, **options):

        vc = Variantcollection.objects.get(sample=options["variant_collection"])
        pbar = tqdm(list(AntibioticResistance.objects.all()))
        for ar in pbar:
            pbar.set_description(str(ar))
            genotype, support = ar.process_variant_collection(vc)

