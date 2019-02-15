import gzip
import os
import subprocess as sp
from glob import glob

import Bio.SeqIO as bpio
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from BioSQL import BioSeqDatabase
from bioresources.models import Assembly
from biosql.models import Biodatabase, Term, Bioentry, Ontology, BioentryQualifierValue, Dbxref, BioentryDbxref
from tqdm import tqdm
import math
from django.db import transaction
import pandas as pd


class NCBI2SQL():

    def __init__(self):
        self.server = None

    def connect_to_server(self, dbuser, dbpass, dbname, dbdriver="MySQLdb", dbhost="localhost"):

        self.server = BioSeqDatabase.open_database(driver=dbdriver, user=dbuser,
                                            passwd=dbpass, host=dbhost, db=dbname)
        return self.server

    def get_or_create_biodb(self, accession, description=""):
        if accession in self.server:
            return self.server[accession]
        db = self.server.new_database(accession, description=description)
        self.server.commit()

        return db

    def download(self, assembly_name, workdir):
        arr = assembly_name.split("_")
        acc = "_".join(arr[:2])
        name = "_".join(arr[2:])
        asspath = "/".join([acc[0:3], acc[4:7], acc[7:10], acc[10:13],
                            acc + "_" + name.replace(" ", "_").replace("#", "_")])

        cmd = 'rsync --recursive --include="*_genomic.gbff.gz" --exclude="*"   rsync://ftp.ncbi.nlm.nih.gov/genomes/all/' + asspath + '/ "' + \
              workdir + '"'
        with open(os.devnull, 'w') as FNULL:
            sp.call(cmd, shell=True,stdout=FNULL)

    def create_contigs(self, accession, genebank,description=""):
        if not Biodatabase.objects.filter(name=accession).exists():
            db = self.get_or_create_biodb(accession, description)
            with tqdm(bpio.parse(gzip.open(genebank, "rt"), "genbank")) as pbar:
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
                    self.server.commit()
        else:
            print("%s already loaded in the database" % accession)

    def create_proteins(self, accession):
        if not Biodatabase.objects.filter(name=accession + "_prots").exists():
            db_seqs = self.server[accession]
            db_prot = self.get_or_create_biodb(accession + "_prots")

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
                        self.server.commit()
        else:
            print("%s proteome already exists" % accession)

    def protein_fasta(self, accession, faa_path, ovewrite=False):
        if ovewrite or not os.path.exists(faa_path):
            db_prot = self.get_or_create_biodb(accession + "_prots")
            with tqdm(db_prot) as pbar:
                def iter():
                    for x in pbar:
                        yield db_prot[x]

                bpio.write(iter(), faa_path, "fasta")
        else:
            print("%s already exists" % faa_path)
