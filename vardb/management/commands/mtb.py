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

from vardb.models import Protocol, AntibioticResistance, ReportedAllele, Allele, Effect, Variant, AlleleEffect
from biosql.models import Term, Ontology, Biodatabase, Bioentry

from tqdm import tqdm
import Bio.SeqIO as bpio
from Bio.Seq import Seq
from Bio.SeqUtils import seq3


class Command(BaseCommand):
    help = 'Loads the obo files to the database. Prepared for GO and SO'

    def __init__(self, stdout=None, stderr=None, no_color=False):
        super().__init__(stdout=stdout, stderr=stderr, no_color=no_color)
        self.resist = None

    def add_arguments(self, parser):
        pass

    def handle(self, *args, **options):
        from Bio.Data import CodonTable
        codontable = CodonTable.ambiguous_dna_by_id[1]
        aa2codon = codontable.back_table
        # aa2aa = defaultdict(lambda : defaultdict(lambda:[]))
        # for aa in table.protein_alphabet.letters:
        #     for i1 in ["A","C","G","T"]:
        #         for i2 in ["A","C","G","T"]:
        #             for i3 in ["A","C","G","T"]:
        #                 aa2aa[aa]

        ann = "/mnt/data2/data/projects/mtbxdr/external/GCF_000195955.2_ASM19595v2_genomic.gbff"
        ann = bpio.read(ann, "gb")

        def search_cds(locus_tag):
            return [x for x in ann.features if x.type == "gene" and locus_tag in x.qualifiers["locus_tag"]][0]

        path_db = "/home/eze/Downloads/andytb - andytb (1).csv"
        self.resist = pd.read_csv(path_db)
        self.resist["AApos"] = [int(x) if x != "-" else "" for x in self.resist.AApos]
        self.resist["LocusTag"] = [x.split("_")[0] for x in self.resist.GeneID]
        self.resist["NucleotidePosH37"] = [int(str(x).split("/")[0]) if str(x) not in ["-", "nan"] else "" for x in
                                           self.resist.NucleotidePosH37]

        # {'AMINOGLYCOSIDES', 'ETHAMBUTOL', 'CAPREOMYCIN', 'STREPTOMYCIN', 'ISONIAZID', 'PYRAZINAMIDE', 'FLUOROQUINOLONES', 'FOSFOMYCIN', 'ETHIONAMIDE', 'KANAMYCIN', 'AMIKACIN', 'LINEZOLID', 'PARA-AMINOSALISYLIC_ACID', 'CLOFAZIMINE', 'BEDAQUILINE', 'VIOMYCIN', 'RIFAMPICIN', '-'}

        errores = []
        idx = 0
        from Bio.Data import CodonTable
        codontable = CodonTable.standard_dna_table
        codon2aa = codontable.forward_table

        # for ridx, r in tqdm(self.resist.iterrows(), total=len(self.resist)):
        #     if r.NucleotidePosH37:
        #         assert r.REF
        #         assert r.ALT
        #     else:
        #         idx += 1
        #         if r.Effect in ["S531", "+nt420:GG", "+nt349:CACTG"]:
        #             continue
        #         cds = search_cds(r.LocusTag)
        #         record = cds.extract(ann)
        #         if r.AAref != str(record.seq.translate())[r.AApos - 1]:
        #             errores.append(r)
        #         else:
        #             muts = []
        #             if cds.strand == -1:
        #
        #                 subseq = ann.seq[cds.location.end - (r.AApos * 3): cds.location.end - ((r.AApos) * 3) + 3]
        #                 assert r.AAref == codon2aa[subseq.reverse_complement()]
        #                 for i in range(3):
        #                     for nu in ["A", "C", "G", "T"]:
        #                         replace = subseq[:i] + nu + subseq[i + 1:]
        #                         if (replace.reverse_complement() in codon2aa and codon2aa[
        #                             replace.reverse_complement()] == r.AAalt):
        #                             muts.append([subseq[i], i + cds.location.end - (r.AApos * 3), str(replace[i]),replace])
        #                         elif replace in codontable.stop_codons and r.AAalt == "STOP":
        #                             muts.append([subseq[i], i + cds.location.end - (r.AApos * 3), str(replace[i]),replace,"muahahha"])
        #
        #                 if muts:
        #                     for mut in muts:
        #                         seq2 = ann.seq[:mut[1]] + mut[2] + ann.seq[mut[1] + 1:]
        #                         seq3 = cds.extract(seq2)
        #                         assert seq3.translate()[r.AApos - 1] == r.AAalt
        #
        #
        #             else:
        #
        #                 start = cds.location.start + (r.AApos - 1) * 3
        #                 subseq = ann.seq[start: start + 3]
        #                 assert r.AAref == codon2aa[subseq]
        #                 for i in range(3):
        #                     for nu in ["A", "C", "G", "T"]:
        #                         replace = str(subseq[:i] + Seq(nu) + subseq[i + 1:])
        #                         if replace in codon2aa and codon2aa[replace] == r.AAalt:
        #                             muts.append([subseq[i], i + ((r.AApos - 1) * 3) + cds.location.start, replace[i]])
        #                         elif replace in codontable.stop_codons and r.AAalt == "STOP":
        #                             muts.append([subseq[i], i + ((r.AApos - 1) * 3) + cds.location.start, replace[i]])
        #                 if muts:
        #                     for mut in muts:
        #                         seq2 = ann.seq[:mut[1]] + mut[2] + ann.seq[mut[1] + 1:]
        #                         seq3 = cds.extract(seq2)
        #                         assert seq3.translate()[r.AApos - 1] == r.AAalt.replace("STOP", "*"), [
        #                             seq3.translate()[r.AApos - 1], r.AAalt]
        #
        #             # if record.seq.translate()[r.AApos ] != r.AAref:
        #             # print ()
        #             # print ([record.seq.translate()[r.AApos ] , r.AAref,r.AApos,r.LocusTag,   r])
        #             # break
        #
        #         assert r.AApos, r
        #         assert r.AAref
        #         assert r.AAalt
        # print([len(errores), idx])

        bdb = Biodatabase.objects.get(name="GCF_000195955.2")
        seq = Bioentry.objects.get(biodatabase=bdb, identifier="NC_000962.3")
        cache_lt = {}
        def search_gene(locustag):
            if locustag in cache_lt:
                 return cache_lt[locustag]
            else:
                f =  Seqfeature.objects.get(Q(bioentry=seq, type_term__name="gene") &
                                    Q(qualifiers__value=locustag,
                                   qualifiers__term__name="locus_tag"))
                cache_lt[locustag] = f
                return f


        map_antibiotic_name = {x.antibiotic.name: x for x in
                               AntibioticResistance.objects.filter(antibiotic__ontology__name="Antibiotics")}

        for i, r in tqdm(self.resist.iterrows(), total=len(self.resist)):

            if r.Drug.lower() not in ["-"]:

                pheno = (map_antibiotic_name[r.Drug.lower()] if r.Drug.lower() in map_antibiotic_name
                            else map_antibiotic_name[r.Drug.lower()[:-1]])
                ra = ReportedAllele(phenotype=pheno, reported_in=r.Source)
                gene = search_gene(r.LocusTag)

                with transaction.atomic():

                    try:
                        cds = search_cds(r.LocusTag)
                    except:
                        self.stderr.write("locus tag not found: '%s'\n" % r.LocusTag  )
                        continue
                    record = cds.extract(ann)
                    allele = None
                    new_effect = None

                    if r.NucleotidePosH37:

                        variant_query = Variant.objects.filter(
                            pos=r.NucleotidePosH37, contig=seq, ref=r.REF)

                        if not variant_query.exists():
                            gene = Seqfeature.objects.get(Q(bioentry=new_variant.contig, type_term__name="gene") &
                                                   Q(qualifiers__value=r.LocusTag,
                                                     qualifiers__term__name="locus_tag"))

                            new_variant = Variant(contig=seq, pos=r.NucleotidePosH37, ref=r.REF,gene=gene)
                            new_variant.save()
                        else:
                            new_variant = variant_query.get()

                        allele_query = Allele.objects.filter(variant_fk=new_variant, alt=r.ALT)
                        if not allele_query.exists():
                            if cds.strand == -1:
                                genepos = cds.location.end - r.NucleotidePosH37
                            else:
                                genepos = r.NucleotidePosH37 - cds.location.start

                            hgvs_c = str(genepos) + r.REF + ">" + r.ALT
                            allele = Allele(variant_fk=new_variant, alt=r.ALT, hgvs_c=hgvs_c)
                            allele.save()

                            if genepos < 0:
                                new_effect = Effect(transcript=r.LocusTag,ref_organism=bdb,gene=gene,
                                                    variant_type='upstream_gene_variant')
                                new_effect.save()
                                allele.main_effect_fk = new_effect
                                allele.save()
                                AlleleEffect(allele_fk=allele, effect_fk=new_effect).save()

                        else:
                            allele = allele_query.get()

                            ra.allele = allele
                            ra.save()

                    if r.AApos:
                        if r.Effect in ["S531", "+nt420:GG", "+nt349:CACTG"]:
                            continue

                        if r.AAref == str(record.seq.translate())[r.AApos - 1]:

                            effect_query = Effect.objects.filter(transcript=r.LocusTag,
                                                                 aa_pos=r.AApos, aa_ref=r.AAref, aa_alt=r.AAalt)
                            if not effect_query.exists():
                                variant_type = 'frameshift_variant' if "fs" == r.AAalt else (
                                    'stop_gained' if "STOP" == r.AAalt else "missense_variant"
                                    if r.AAref != r.AAalt else 'synonymous_variant')

                                new_effect = Effect(transcript=r.LocusTag,ref_organism=bdb,
                                                    variant_type=variant_type)
                                new_effect.aa_pos = r.AApos
                                new_effect.aa_ref = r.AAref
                                new_effect.aa_alt = "*" if r.AAalt in ["fs", "STOP"] else r.AAalt
                                new_effect.hgvs_p = (
                                        (seq3(new_effect.aa_ref) if new_effect.aa_ref != "*" else "*") +
                                        str(new_effect.aa_pos) +
                                        (seq3(new_effect.aa_alt) if new_effect.aa_alt != "*" else "*"))
                                new_effect.save()
                            else:
                                new_effect = effect_query.get()
                            ra.effect = new_effect
                            ra.save()
                            if r.NucleotidePosH37:
                                allele.main_effect_fk = new_effect
                                allele.save()
                                AlleleEffect(allele_fk=allele, effect_fk=new_effect).save()
