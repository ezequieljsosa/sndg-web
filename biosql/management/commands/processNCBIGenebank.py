import gzip
import os
import subprocess as sp
from glob import glob

import Bio.SeqIO as bpio
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from BioSQL import BioSeqDatabase
from django.conf import settings
from django.core.management.base import BaseCommand

from bioresources.models import Assembly
from biosql.models import Biodatabase, Term, Bioentry, Ontology, BioentryQualifierValue, Dbxref, BioentryDbxref
from tqdm import tqdm
import math
from django.db import transaction
import pandas as pd


class Command(BaseCommand):
    help = 'Loads and '

    def __init__(self, stdout=None, stderr=None, no_color=False):
        super().__init__(stdout=stdout, stderr=stderr, no_color=no_color)
        self.ontology = None
        self.dbmap = {}
        self.is_a = None
        self.relmap = {}
        self.ontology_name = ""
        self.obo_path = None

    def add_arguments(self, parser):
        parser.add_argument('--accession', required=True)
        parser.add_argument('--driverDB', default="MySQLdb")
        parser.add_argument('--skipAnn', action='store_false')
        parser.add_argument('--workDir', default='/tmp')

    def download(self, assembly, options):
        arr = assembly.name.split("_")
        acc = "_".join(arr[:2])
        name = "_".join(arr[2:])
        asspath = "/".join([acc[0:3], acc[4:7], acc[7:10], acc[10:13],
                            acc + "_" + name.replace(" ", "_").replace("#", "_")])

        cmd = 'rsync --recursive --include="*_genomic.gbff.gz" --exclude="*"   rsync://ftp.ncbi.nlm.nih.gov/genomes/all/' + asspath + '/ "' + \
              options["workDir"] + '"'
        sp.call(cmd, shell=True)  # , stdout=FNULL

    def connect_to_server(self, options):
        if "user" in options:
            user = options["user"]
        else:
            data = {x.split("=")[0]: x.split("=")[1].strip() for x in
                    open(os.path.expanduser(settings.DATABASES["default"]["OPTIONS"]["read_default_file"])).readlines()[
                    1:]}
            user = data["user"]

        if "pass" in options:
            password = options["pass"]
        else:
            data = {x.split("=")[0]: x.split("=")[1].strip() for x in
                    open(os.path.expanduser(settings.DATABASES["default"]["OPTIONS"]["read_default_file"])).readlines()[
                    1:]}
            password = data["password"]

        dbname = settings.DATABASES["default"]["NAME"]
        return BioSeqDatabase.open_database(driver=options["driverDB"], user=user,
                                            passwd=password, host="localhost", db=dbname)

    def get_or_create_biodb(self, accession, description, server):
        if accession in server:
            return server[accession]
        db = server.new_database(accession, description=description)
        server.commit()
        return db

    def handle(self, *args, **options):
        options["workDir"] = os.path.abspath(options["workDir"])
        assert os.path.exists(options["workDir"]), "%s does not exists" % options["workDir"]
        assemblyqs = Assembly.objects.filter(external_ids__identifier=options["accession"])
        if not assemblyqs.exists():
            self.stderr.write(options["accession"] + " not found")
            return 1
        assembly = assemblyqs.first()

        self.download(assembly, options)
        filename = glob(options["workDir"] + "/" + options["accession"] + "*.gbff.gz")[0]
        server = self.connect_to_server(options)
        if not Biodatabase.objects.filter(name=options["accession"]).exists():
            db = self.get_or_create_biodb(options["accession"], assembly.description, server)
            with tqdm(bpio.parse(gzip.open(filename, "rt"), "genbank")) as pbar:
                pbar.set_description("Processing contigs: ")
                for sequence in pbar:
                    features = []
                    for f in sequence.features:
                        try:
                            int(f.location.start)
                            int(f.location.end)
                            features.append(f)
                        except ValueError:
                            pass
                    sequence.features = features
                    db.load([sequence])
                    server.commit()
        if not Biodatabase.objects.filter(name=options["accession"] + "_prots").exists():
            db_seqs = server[options["accession"]]
            db_prot = self.get_or_create_biodb(options["accession"] + "_prots", "", server)

            with tqdm(db_seqs) as pbar:
                pbar.set_description("Processing contigs")
                for seqid in pbar:
                    sequence = db_seqs[seqid]
                    prots = []
                    for f in sequence.features:
                        if "translation" in f.qualifiers and "locus_tag" in f.qualifiers:
                            r = SeqRecord(id=f.qualifiers["locus_tag"][0], name="", description="",
                                          seq=Seq(f.qualifiers["translation"][0]))
                            prots.append(r)
                    with tqdm(prots) as pbar2:
                        pbar2.set_description("Saving contig proteins")
                        db_prot.load(pbar2)
                        server.commit()

        faa_path = options["workDir"] + "/" + options["accession"] + ".faa"
        if not os.path.exists(faa_path):
            db_prot = self.get_or_create_biodb(options["accession"] + "_prots", "", server)
            with tqdm(db_prot) as pbar:
                def iter():
                    for x in pbar:
                        yield db_prot[x]

                bpio.write(iter(), faa_path, "fasta")
        self.annotate(options["accession"], faa_path)

        self.jbrowse(options)
        self.krona(options)

    def krona(self):
        pass

    def jbrowse(self):
        pass

    def annotate(self, accession, faa_path):

        ann_eggnog = faa_path + ".emapper.annotations"
        if not os.path.exists(ann_eggnog):
            "run emmaper"

        cols = "query_name,seed_eggNOG_ortholog,seed_ortholog_evalue,seed_ortholog_score,predicted_gene_name,GO_terms,KEGG_KOs,BiGG_reactions,Annotation_tax_scope,OGs,bestOG|evalue|score,COG cat,eggNOG annot"
        df_egg = pd.read_table(ann_eggnog, comment="#", names=cols.split(","))
        GO_ROOT_TERMS = ["GO:0008150", "GO:0005575", "GO:0003674"]
        go_ontology = Ontology.objects.filter(name="Gene Ontology").get()
        df_egg = df_egg.fillna(value=False)
        for _, row in tqdm(df_egg.iterrows(), total=len(df_egg)):
            with transaction.atomic():
                entry = Bioentry.objects.filter(biodatabase__name=accession + "_prots",
                                                accession=row["query_name"]).get()
                entry.description = row['eggNOG annot']

                dbxref = Dbxref.objects.get_or_create(dbname="eggnog", accession=row["seed_eggNOG_ortholog"])[0]
                BioentryDbxref.objects.get_or_create(bioentry=entry, dbxref=dbxref)

                if row["predicted_gene_name"]:
                    entry.name = row["predicted_gene_name"]
                if row["GO_terms"]:
                    for go in row["GO_terms"].split(","):
                        if go not in GO_ROOT_TERMS:
                            gotermqs = Term.objects.filter(ontology=go_ontology, identifier=go)
                            if gotermqs.exists():
                                go_term = gotermqs.get()
                                BioentryQualifierValue(bioentry=entry, term=go_term).save()
                            else:
                                self.stderr.write("%s term not found" % go)

                if row["KEGG_KOs"]:
                    for ko in row["KEGG_KOs"].split(","):
                        dbxref = Dbxref.objects.get_or_create(dbname="kegg_ko", accession=ko)[0]
                        BioentryDbxref.objects.get_or_create(bioentry=entry, dbxref=dbxref)
