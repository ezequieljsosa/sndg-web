import json
import os
import subprocess as sp

from bioseq.models import Taxon
from django.core.management.base import BaseCommand
from django.db import transaction
from tqdm import tqdm

from bioresources.models import Barcode, Organization


def execute(cmd,**kwargs):
    sp.call(cmd.format(**kwargs),shell=True)

def download_file(complete_url, target, ovewrite=False, retries=3):
    if not target.strip():
        target = "./"
    if not os.path.exists(os.path.dirname(os.path.abspath(target))):
        raise Exception("%s does not exists" % os.path.dirname(target))
    if os.path.exists(target) and not ovewrite:
        raise OvewriteFileException("%s already exists" % target)

    execute(' wget  --timeout=20 --tries={retries} -O {target} "{url}"',
            url=complete_url, retries=retries, target=target)


class Command(BaseCommand):
    help = 'Download and loads all barcodes of Bold of a given country'

    def __init__(self, stdout=None, stderr=None, no_color=False):
        super().__init__(stdout=stdout, stderr=stderr, no_color=no_color)
        self.ontology = None
        self.dbmap = {}
        self.is_a = None
        self.relmap = {}
        self.ontology_name = ""
        self.obo_path = None

    def add_arguments(self, parser):

        parser.add_argument('--country', required=True)
        parser.add_argument('--bold_url',
                            default="http://www.boldsystems.org/index.php/API_Public/combined?format=json&geo=")
        parser.add_argument('--json', default="/tmp/bold.json")

    def handle(self, *args, **options):

        if not os.path.exists(options["json"]):
            download_file(options["bold_url"] + options["country"], options["json"])

        data = json.load(open(options["json"]))["bold_records"]["records"].values()
        bcodes = []
        org_bold = Organization.objects.get_or_create(name="BOLD")[0]

        for d in tqdm(data):
            if "sequences" in d:
                for x in ["species", "genus", "subfamily", "family", "order", "class", "phylum"]:
                    if x in d["taxonomy"]:
                        tax = d["taxonomy"][x]["taxon"]["taxID"]
                        break

                d["description"] = (d["taxonomy"][x]["taxon"]["name"] + " " +
                                    d["sequences"]["sequence"][0]["markercode"] + " "
                                    + d["specimen_identifiers"]["institution_storing"])

                d["tax"] = int(tax)

                bc = Barcode(
                    name=d["processid"],
                    description=d["description"],
                    country=d["collection_event"]["country"],
                    collectors=d["collection_event"]["collectors"],
                    type="barcode",
                    marker=d["sequences"]["sequence"][0]["markercode"],

                    bold_org=d["specimen_identifiers"]["institution_storing"]
                )
                try:
                    bc.ncbi_tax = Taxon.objects.get(ncbi_taxon_id=d["tax"])
                except Taxon.DoesNotExist:
                    pass
                try:
                    bc.image_url = d["specimen_imagery"]["media"][0]["image_file"]
                except KeyError:
                    pass
                try:
                    bc.subdivision = d["collection_event"]["province_state"]
                except KeyError:
                    pass

                bcodes.append(bc)
        bcodes_col = []
        for i,bc in enumerate(tqdm(bcodes)):
            bcodes_col.append(bc)
            if i == 5000:
                with transaction.atomic():
                    for bc2 in bcodes_col:
                        bc2.save()
                        bc2.publishers.add(org_bold)
                    bcodes_col = []

        with transaction.atomic():
            for bc2 in bcodes_col:
                bc2.save()
                bc2.publishers.add(org_bold)

