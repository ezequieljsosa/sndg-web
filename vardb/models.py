from django.db import models
from django.db.models import Q
from biosql.models import Biodatabase, Bioentry, Seqfeature, Term, Ontology
from .managers import VariantannotationManager, ReportedAlleleManager, VariantassignmentManager


class Allele(models.Model):
    id = models.AutoField(primary_key=True)
    variant_fk = models.ForeignKey('Variant', models.DO_NOTHING, db_column='variant_fk', related_name="alleles")
    alt = models.CharField(max_length=255)
    main_effect_fk = models.ForeignKey('Effect', models.DO_NOTHING, db_column='main_effect_fk', blank=True, null=True)
    hgvs_c = models.CharField(max_length=255, blank=True, null=True)

    class Meta:
        unique_together = (('variant_fk', 'alt'),)

    def __str__(self):
        return str(self.variant_fk) + "->" + self.alt


class AlleleEffect(models.Model):
    allele_fk = models.ForeignKey(Allele, models.CASCADE, db_column='allele_fk', related_name="effects")
    effect_fk = models.ForeignKey('Effect', models.CASCADE, db_column='effect_fk', related_name="alleles")


class Effect(models.Model):
    id = models.AutoField(primary_key=True)

    transcript = models.CharField(max_length=255, blank=True, null=True, db_index=True)
    gene = models.ForeignKey(Seqfeature, models.CASCADE, related_name="effects", null=True)
    variant_type = models.CharField(max_length=255)

    predicted_impact = models.CharField(max_length=255, blank=True, null=True)
    aa_ref = models.CharField(max_length=255, blank=True, null=True)
    aa_pos = models.IntegerField(blank=True, null=True)
    aa_alt = models.CharField(max_length=255, blank=True, null=True)
    hgvs_p = models.CharField(max_length=255, blank=True, null=True)

    ref_organism = models.ForeignKey(Biodatabase, models.CASCADE, "variant_effects", db_index=True)

    def __str__(self):
        return self.transcript + " " + self.variant_type + " " + self.hgvs_p


class Variant(models.Model):
    id = models.AutoField(primary_key=True)
    pos = models.IntegerField()
    # gene = models.CharField(max_length=255, blank=True, null=True)
    gene = models.ForeignKey(Seqfeature, models.CASCADE, related_name="variants", null=True)
    gene_pos = models.IntegerField(blank=True, null=True)

    description = models.CharField(max_length=255, blank=True, null=True)

    contig = models.ForeignKey(Bioentry, models.CASCADE, related_name="variants")
    ref = models.CharField(max_length=255)

    created_at = models.DateTimeField(auto_now_add=True)
    updated_at = models.DateTimeField(auto_now=True)

    def __str__(self):
        return self.contig.name + " " + self.ref + str(self.pos)

    class Meta:
        unique_together = (('contig', 'pos', 'ref'),)
        index_together = [
            ["contig", "gene"],
        ]


class Variantannotation(models.Model):
    DP = "DP"

    id = models.AutoField(primary_key=True)
    assignment_fk = models.ForeignKey('Variantassignment', models.DO_NOTHING, db_column='assignment_fk',
                                      related_name="annotations")
    source_type = models.CharField(max_length=255)
    source = models.CharField(max_length=255)
    prop = models.CharField(max_length=255)
    value = models.CharField(max_length=255)
    description = models.CharField(max_length=255, blank=True, null=True)

    objects = VariantannotationManager()

    def __str__(self):
        return self.prop + "=" + self.value


class Variantassignment(models.Model):
    id = models.AutoField(primary_key=True)
    variant_collection_fk = models.ForeignKey('Variantcollection', models.CASCADE, related_name="assignments",
                                              db_column='variant_collection_fk')
    variant_fk = models.ForeignKey(Variant, models.CASCADE, db_column='variant_fk', related_name="assignments")
    allele_fk = models.ForeignKey(Allele, models.CASCADE, db_column='allele_fk', related_name="assignments")

    objects = VariantassignmentManager()

    class Meta:
        unique_together = (('variant_collection_fk', 'variant_fk', 'allele_fk'),)


class Variantcollection(models.Model):
    id = models.AutoField(primary_key=True)
    ref_organism = models.ForeignKey(Biodatabase, models.CASCADE, "strains", null=True)
    sample = models.CharField(max_length=255)
    description = models.CharField(max_length=255, blank=True, null=True)
    created_at = models.DateTimeField(auto_now_add=True)
    updated_at = models.DateTimeField(auto_now=True)

    class Meta:
        unique_together = (('ref_organism', 'sample'),)


class Phenotype(models.Model):
    id = models.AutoField(primary_key=True)
    name = models.CharField(max_length=255)
    description = models.TextField()


class AntibioticResistance(Phenotype):
    antibiotic = models.ForeignKey(Term, models.DO_NOTHING)

    def __str__(self):
        return self.antibiotic.identifier

    def __repr__(self):
        return str(self)

    def process_variant_collection(self, variant_collection, depth=30):
        """
        processes a variant collection and creates a GenotypeSupport if needed
        :param vc: variant collection
        :return:
        """
        genotype = Genotype.objects.filter(phenotype=self, variant_collection=variant_collection)
        if genotype.exists():
            genotype = genotype.get()
        else:
            genotype = Genotype(phenotype=self, variant_collection=variant_collection)
        support = []
        status_ontology = Ontology.objects.get(name=GenotypeSupport.STATUS_ONTOLOGY)
        status_conclusive = Term.objects.get(ontology=status_ontology, name="Conclusive")
        status_possible = Term.objects.get(ontology=status_ontology, name="Possible")
        status_hint = Term.objects.get(ontology=status_ontology, name="Hint")

        low_conf_alleles = [x.assignment_fk.allele_fk for x in
                            Variantannotation.objects.dp_filter(variant_collection, depth).prefetch_related(
                                "assignment_fk__allele_fk")]

        exact_reported = list(ReportedAllele.objects.exact_reported(self, variant_collection).exclude(
            allele__in=low_conf_alleles))

        support += self.get_genotype_support(exact_reported, variant_collection, status_conclusive, genotype)

        pos_reported = [x for x in ReportedAllele.objects.pos_reported(self, variant_collection).exclude(
            allele__in=low_conf_alleles) if x not in exact_reported]

        support += self.get_genotype_support(pos_reported, variant_collection, status_possible, genotype)

        reported_genes = self.reported_genes()
        gene_variants = Variantassignment.objects.reported_positions(variant_collection,
                                                                     reported_genes).exclude(
            allele_fk__in=low_conf_alleles)

        for assignment in gene_variants:
            if assignment not in exact_reported and assignment not in pos_reported:

                if (assignment.allele_fk.main_effect_fk and
                        assignment.allele_fk.main_effect_fk.variant_type == "frameshift_variant"):
                    gene = [x for x in reported_genes if x.allele == assignment.allele_fk][0]
                    variant_fs = [x for x in gene.effects.all() if "frameshift_variant" in x.variant_type]

                    if variant_fs:
                        gs = GenotypeSupport(genotype=genotype, status=status_conclusive, assignment=assignment,
                                             reported_genes=variant_fs[0].reported.first())
                    else:
                        gs = GenotypeSupport(genotype=genotype, status=status_possible, assignment=assignment)

                else:
                    gs = GenotypeSupport(genotype=genotype, status=status_hint, assignment=assignment)

                support.append(gs)
        return genotype, support

    def get_genotype_support(self, reported_qs, variant_collection, status_term, genotype):
        support = []
        for reported in reported_qs:
            if reported.allele:

                assignment = Variantassignment.objects.get(variant_collection_fk=variant_collection,
                                                       variant_fk=reported.allele.variant_fk)
            else:
                alleles = [x.allele_fk for x in reported.effect.alleles.all()]
                assignment = Variantassignment.objects.get(variant_collection_fk=variant_collection,
                                                           allele_fk__in=alleles)

            gs = GenotypeSupport(genotype=genotype, status=status_term, reported_allele=reported, assignment=assignment)
            support.append(gs)
        return support

    def reported_genes(self):

        return Seqfeature.objects.prefetch_related("qualifiers__term").filter(
            Q(variants__alleles__main_effect_fk__reported__phenotype=self) |
            Q(variants__alleles__reported__phenotype=self)).distinct()


class Protocol(models.Model):
    id = models.AutoField(primary_key=True)
    phenotype = models.ForeignKey(Phenotype, models.CASCADE, "protocols", null=True)
    name = models.CharField(max_length=255)


class Assay(models.Model):
    id = models.AutoField(primary_key=True)
    protocol = models.ForeignKey(Protocol, models.CASCADE, "assays", null=True)
    variant_collection = models.ForeignKey(Variantcollection, models.CASCADE)
    created_at = models.DateTimeField(auto_now_add=True)
    updated_at = models.DateTimeField(auto_now=True)
    result = models.ForeignKey(Term, models.DO_NOTHING)


class Genotype(models.Model):
    id = models.AutoField(primary_key=True)
    phenotype = models.ForeignKey(Phenotype, models.CASCADE, "genotypes", null=True)
    variant_collection = models.ForeignKey(Variantcollection, models.CASCADE)


class ReportedAllele(models.Model):
    id = models.AutoField(primary_key=True)
    allele = models.ForeignKey(Allele, models.CASCADE, null=True, related_name="reported")
    effect = models.ForeignKey(Effect, models.CASCADE, null=True, related_name="reported")
    phenotype = models.ForeignKey(Phenotype, models.CASCADE, "reported", null=True)
    reported_in = models.CharField(max_length=255)
    created_at = models.DateTimeField(auto_now_add=True)
    updated_at = models.DateTimeField(auto_now=True)

    objects = ReportedAlleleManager()


class GenotypeSupport(models.Model):
    STATUS_ONTOLOGY = "Genotype Support Status"

    id = models.AutoField(primary_key=True)
    genotype = models.ForeignKey(Genotype, models.CASCADE, "supportedby")
    status = models.ForeignKey(Term, models.DO_NOTHING)
    reported_allele = models.ForeignKey(ReportedAllele, models.DO_NOTHING, null=True)
    assignment = models.ForeignKey(Variantassignment, models.DO_NOTHING)
