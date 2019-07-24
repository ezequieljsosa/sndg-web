# -*- coding: utf-8 -*-
from __future__ import unicode_literals
import os

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import Bio.SeqIO as bpio

from bioseq.models.Bioentry import Bioentry
from bioseq.models.Seqfeature import Seqfeature
from collections import defaultdict


class Model2JBrowse():

    def __init__(self):
        self.name = ""
        self.bioentry_iterator = None
        self.fasta_path = "/tmp/jbrowse.fasta"
        self.gff_path = "/tmp/jbrowse.gff"
        self.ovewrite = False
        self.source = "."
        self.unique_ids = defaultdict(lambda: 0)

    def get_next_id(self, feature_type: str):
        idx = self.unique_ids[feature_type]
        self.unique_ids[feature_type] += 1
        return feature_type + "_" + str(idx)

    def build(self):
        assert os.path.exists(os.path.dirname(self.fasta_path))
        assert os.path.exists(os.path.dirname(self.gff_path))
        assert self.bioentry_iterator, "No bioentries selected..."
        if not os.path.exists(os.path.dirname(self.fasta_path)) or self.ovewrite:
            with open(self.fasta_path, "w") as h_fasta, open(self.gff_path, "w") as h_gff:
                for bioentry in self.bioentry_iterator:
                    self.process_bioentry(bioentry, h_fasta, h_gff)

    def process_bioentry(self, bioentry: Bioentry, h_fasta, h_gff):
        r = SeqRecord(id=bioentry.accession, name="", description="", seq=Seq(bioentry.seq.seq))
        bpio.write(r, h_fasta, "fasta")

        for feature in bioentry.features.all():
            self.process_seqfeature(bioentry.accession, feature, h_gff)

    def process_seqfeature(self, accession: str, feature: Seqfeature, h_gff):
        loc = feature.locations.all()[0]
        fid = self.get_next_id(feature.type_term.identifier)
        attributes = ["ID=" + fid]
        qualifiers = {x.term.identifier: x.value for x in feature.qualifiers.all()}
        if "gene" in qualifiers:
            attributes.append("Name=" + qualifiers["gene"])
        if "locus_tag" in qualifiers:
            attributes.append("locus_tag=" + qualifiers["locus_tag"])

        attributes_str = ",".join(attributes)
        r = "\t".join(
            [str(x) for x in [accession, self.source, feature.type_term.identifier, loc.start_pos, loc.end_pos,
                              ".", loc.strand, 0, attributes_str]])
        h_gff.write(r + "\n")

        for sf in feature.subfeatures():
            self.process_subfeature(accession, sf, h_gff, fid)

    def process_subfeature(self, accession: str, feature: Seqfeature, h_gff, parent):
        loc = feature.locations.all()[0]
        sfs = feature.subfeatures()
        attributes = ["Parent=" + parent]
        if "locus_tag" in qualifiers:
            attributes.append("locus_tag=" + qualifiers["locus_tag"])
        if sfs:
            fid = self.get_next_id(feature.type_term.identifier)
            attributes.append("ID=" + fid)
        attributes_str = ",".join(attributes)
        r = "\t".join(
            [str(x) for x in [bioentry.accession, self.source, feature.type_term.identifier, loc.start_pos, loc.end_pos,
                              ".", loc.strand, 0, attributes_str]])
        h_gff.write(r + "\n")
        if sfs:
            for sf in sfs:
                self.process_subfeature(accession, sf, h_gff, fid)

    """
    docker run -v /data/xomeq/jbrowse/data/:/jbrowse/data/ -v /tmp/jbrowse.fasta:/tmp/jbrowse.fasta  jbrowse/jbrowse-1.12.0 bin/prepare-refseqs.pl --fasta /tmp/jbrowse.fasta --out data/TestBacteria --key "Sequence"

    """