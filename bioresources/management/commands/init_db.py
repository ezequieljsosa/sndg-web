import os
import subprocess as sp
from tqdm import tqdm
import requests
from datetime import datetime

from SNDG.WebServices import download_file

from django.core.management.base import BaseCommand

from bioseq.models.Ontology import Ontology
from bioseq.models.Term import Term, TermRelationship
from bioseq.models.Biodatabase import Biodatabase
from bioseq.models.Bioentry import Bioentry, BioentryQualifierValue, BioentryDbxref
from bioseq.models.Biosequence import Biosequence
from bioseq.models.Seqfeature import Seqfeature, SeqfeatureDbxref, SeqfeatureQualifierValue
from bioseq.models.Location import Location
from bioseq.models.Dbxref import Dbxref, DbxrefQualifierValue

from bioresources.models.Resource import Resource
from bioresources.models.Person import Person
from bioresources.models.Organization import Organization
from bioresources.models.Affiliation import Affiliation
from bioresources.models.Identity import Identity
from bioresources.models.ExternalId import ExternalId
from bioresources.models.Publication import Publication
from bioresources.models.Tool import Tool
from bioresources.models.Assembly import Assembly
from bioresources.models.ReadsAchive import ReadsArchive
from bioresources.models.BioProject import BioProject
from bioresources.models.Sample import Sample
from bioresources.models.ResourceRelation import ResourceRelation
from bioresources.models.ResourceProperty import ResourceProperty, ResourcePropertyValue


# from pdbdb.models import PDB


class InitDB(BaseCommand):
    """

    """

    help = 'Tinny DB for documentation and testing purposes'

    def add_arguments(self, parser):
        parser.add_argument('--remove', action='store_true')

    def handle(self, *args, **options):
        assert os.path.exists(options["dst_dir"])
        Ontology.load_ann_terms()
        Ontology.load_go_base(Ontology.GO)
        self.load_ontology()

    def load_ontology(self):
        go = Ontology.objects.get(name=Ontology.GO)

        """
        [Inferred is_a relation] GO:0008150 biological_process
            [Inferred is_a relation] GO:0050896 response to stimulus            
                [Inferred is_a relation] GO:0006950 response to stress
        """

        bp = Term.objects.get_or_create(identifier="GO:0008150", name="biological_process", version=1, ontology=go)[0]
        r2s = Term.objects.get_or_create(identifier="GO:0050896", name="response to stimulus ", version=1, ontology=go)[
            0]
        r2ss = Term.objects.get_or_create(identifier="GO:0006950", name="response to stress", version=1, ontology=go)[0]
        TermRelationship(
            subject_term=bp,  # parent
            predicate_term=Ontology.relmap["is_a"],
            object_term=r2s,  # child
            ontology=go).save()

        TermRelationship(
            subject_term=r2s,  # parent
            predicate_term=Ontology.relmap["is_a"],
            object_term=r2ss,  # child
            ontology=go).save()

    def load_bacteria(self):
        sfk_ontology = Ontology.objects.get(name=Ontology.SFK)
        sfs_ontology = Ontology.objects.get(name=Ontology.SFS)
        ann_ontology = Ontology.objects.get(name=Ontology.ANNTAGS)

        cds = Term.objects.get(identifier="CDS", ontology=sfk_ontology)
        trna = Term.objects.get(identifier="tRNA", ontology=sfk_ontology)
        calculated = Term.objects.get(identifier="calculated", ontology=sfs_ontology)
        gene = Term.objects.get(identifier="gene", ontology=ann_ontology)
        locus_tag = Term.objects.get(identifier="locus_tag", ontology=ann_ontology)
        product = Term.objects.get(identifier="product", ontology=ann_ontology)

        genome = Biodatabase.objects.create(name="TestBacteria", description="Bacteria for testing purposes")
        contig = Bioentry.objects.create(biodatabase_id=genome, description="bcontig1", name="bcontig1")
        seq = """TTGACCGATGACCCCGGTTCAGGCTTCACCACAGTGTGGAACGCGGTCGTCTCCGAACTTAACGGCGACC
                 CTAAGGTTGACGACGGACCCAGCAGTGATGCTAATCTCAGCGCTCCGCTGACCCCTCAGCAAAGGGCTTG
                 GCTCAATCTCGTCCAGCCATTGACCATCGTCGAGGGGTTTGCTCTGTTATCCGTGCCGAGCAGCTTTGTC"""
        Biosequence.objects.create(bioentry=contig, seq=seq.replace("\n", "").replace(" ", ""))
        tRNA = Seqfeature.objects.create(bioentry=contig, type_term=trna, source_term=calculated)
        SeqfeatureQualifierValue.objects.create(seqfeature=tRNA, term=gene, value="alaT")
        SeqfeatureQualifierValue.objects.create(seqfeature=tRNA, term=locus_tag, value="Rvnt02")
        SeqfeatureQualifierValue.objects.create(seqfeature=tRNA, term=product, value="tRNA-Ile")
        Location.objects.create(seqfeature=tRNA, start_pos=10, end_pos=30)
        prot = Seqfeature.objects.create(bioentry=contig, type_term=cds, source_term=calculated)
        SeqfeatureQualifierValue.objects.create(seqfeature=prot, term=locus_tag, value="Rv0009")
        Location.objects.create(seqfeature=prot, start_pos=35, end_pos=55)

        proteome = Biodatabase.objects.create(name="TestBacteria_prots", description="Bacteria for testing purposes")
        prot = Bioentry.objects.create(biodatabase_id=proteome, description="ppiA", name="Rv0009")
        Biosequence.objects.create(bioentry=prot, seq="MADCDSVTNSPLATATA")
        unip = Dbxref(dbname="Uniprot", accession="P9WHW3")
        BioentryDbxref.objects.create(bioentry=prot, dbxref=unip)
        BioentryQualifierValue.objects.create(bioentry=prot, term=product,
                                              value="iron-regulated peptidyl-prolyl cis-transisomerase PpiA")
        BioentryQualifierValue.objects.create(bioentry=prot, term=gene, value="ppiA")

    def load_eukaryote(self):
        pass

    def load_structure(self):
        tmp_dir = "/tmp/2PZI.pdb"
        download_file("https://files.rcsb.org/view/2PZI.pdb", target=tmp_dir, ovewrite=True)

    def load_resources(self):
        mixs_ontology = Ontology.objects.get(name="MIxS")

        abs1 = """The World Health Organization (WHO) estimates that 40% of tuberculosis (TB) cases are not diagnosed and treated correctly. Even though there are several diagnostic tests available in the market, rapid, easy, inexpensive detection, and drug susceptibility testing (DST) of Mycobacterium tuberculosis is still of critical importance specially in low and middle-income countrie"""
        publication_fp = Publication.objects.create(
            name="Fluoromycobacteriophages can detect viable Mycobacterium tuberculosis and determine phenotypic rifampicin resistance in 3-5 days from sputum collection",
            description=abs1,
            doi="10.3389/fmicb.2018.01471", date_of_publication=datetime(2018, 7, 5), pubmed_id="PMC6041418")
        publication_tp = Publication.objects.create(
            name="Target-Pathogen: A structural bioinformatic approach to prioritize drug targets in pathogens",
            doi="10.1093/nar/gkx1015", date_of_publication=datetime(2018, 7, 5))

        eze = Person.objects.create(surname="Sosa", name="Ezequiel Jorge", scopus_id="57057186700")
        adri = Person.objects.create(surname="Turjanski", name="Gustavo Adrián", scopus_id="6602847199")
        spo = Person.objects.create(surname="Poggi", name="Susana")

        ncbi = Organization.objects.create(name="NCBI")
        ic = Organization.objects.create(name="University of Buenos Aires Instituto de Cálculo",
                                         url="http://www.ic.fcen.uba.ar/",
                                         country="Argentina")
        uba = Organization.objects.create(
            name="Facultad de Ciencias Exactas y Naturales de la Universidad en Buenos Aires",
            url="http://www.ic.fcen.uba.ar/", country="Argentina")
        pbi = Organization.objects.create(
            name="Department of Biological Sciences, Pittsburgh Bacteriophage Institute, University of Pittsburgh",
            url="http://www.pitt.edu/~duda/PBImemberlist.html", country="United States")

        Affiliation.objects.create(publication=publication_fp, author=eze, organizations=[ic, uba])
        Affiliation.objects.create(publication=publication_fp, author=spo, organizations=[pbi, uba])
        Affiliation.objects.create(publication=publication_tp, author=eze, organizations=[uba])
        Affiliation.objects.create(publication=publication_tp, author=adri, organizations=[uba])

        bioproject = BioProject.objects.create(
            name="XDR Mycobacterium tuberculosis strains from sputum. Argentine, Muniz/Malbran Hospital",
            description="XDR Mycobacterium tuberculosis strains from sputum. Argentine, Muniz/Malbran Hospital")
        ExternalId.objects.create(organization=ncbi, identifier="PRJNA438689", type="acc", resource=bioproject)
        ExternalId.objects.create(organization=ncbi, identifier="438689", type="id", resource=bioproject)

        sample = Sample.create(name="SAMN08724043",
                               description="Pathogen: clinical or host-associated sample from Mycobacterium tuberculosis")
        ExternalId.objects.create(organization=ncbi, identifier="SAMN08724043", type="acc", resource=sample)
        ExternalId.objects.create(organization=ncbi, identifier="8724043", type="id", resource=sample)
        host = Term.objects.get(identifier="host", ontology=mixs_ontology)
        geo_loc_name = Term.objects.get(identifier="geo_loc_name", ontology=mixs_ontology)
        biovar = Term.objects.get(identifier="biovar", ontology=mixs_ontology)
        p = ResourceProperty.objects.create(organization=ncbi, resource=sample, term=host)
        ResourcePropertyValue.objects.create(property=p, value="human")
        p = ResourceProperty.objects.create(organization=ncbi, resource=sample, term=geo_loc_name)
        ResourcePropertyValue.objects.create(property=p, value="Argentina: Buenos Aires")
        p = ResourceProperty.objects.create(organization=ncbi, resource=sample, term=biovar)
        ResourcePropertyValue.objects.create(property=p, value="XDR-4f")

        ra = ReadsArchive.objects.create(name="SRX3803595")

        ResourceRelation.objects.create(source=publication_fp, target=bioproject,
                                        role=Resource.RESOURCE_TYPES[
                                                 Resource.RESOURCE_TYPES.PUBLICATION] +
                                             "_" + Resource.RESOURCE_TYPES.BIOPROJECT)
        ResourceRelation.objects.create(source=bioproject, target=sample,
                                        role=Resource.RESOURCE_TYPES[
                                                 Resource.RESOURCE_TYPES.BIOPROJECT] +
                                             "_" + Resource.RESOURCE_TYPES.SAMPLE)

        ResourceRelation.objects.create(source=sample, target=ra,
                                        role=Resource.RESOURCE_TYPES[
                                                 Resource.RESOURCE_TYPES.SAMPLE] +
                                             "_" + Resource.RESOURCE_TYPES.READS)

        tp = Tool.objects.create(name="TargetPathogen")
        ResourceRelation.objects.create(source=publication_tp, target=tp,
                                        role=Resource.RESOURCE_TYPES[
                                                 Resource.RESOURCE_TYPES.PUBLICATION] +
                                             "_" + Resource.RESOURCE_TYPES.TOOL)

        # '2', 'PUBLICATION_EXPRESSION', '0', '8', '4990'
        # '78', 'PUBLICATION_ASSEMBLY', '0', '185', '5066'
        # '290', 'PUBLICATION_STRUCTURE', '0', '661', '5278'
        #
