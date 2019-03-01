from django.db import models

from biosql.models import Biodatabase, Bioentry, Seqfeature,Term


class Allele(models.Model):
    id = models.AutoField(primary_key=True)
    variant_fk = models.ForeignKey('Variant', models.DO_NOTHING, db_column='variant_fk')
    alt = models.CharField(max_length=255)
    main_effect_fk = models.ForeignKey('Effect', models.DO_NOTHING, db_column='main_effect_fk', blank=True, null=True)
    hgvs_c = models.CharField(max_length=255, blank=True, null=True)

    class Meta:
        unique_together = (('variant_fk', 'alt'),)


class Effect(models.Model):
    id = models.AutoField(primary_key=True)
    allele_fk = models.ForeignKey(Allele, models.DO_NOTHING, db_column='allele_fk',null=True)
    transcript = models.CharField(max_length=255, blank=True, null=True)
    variant_type = models.CharField(max_length=255)

    predicted_impact = models.CharField(max_length=255, blank=True, null=True)
    aa_ref = models.CharField(max_length=255, blank=True, null=True)
    aa_pos = models.IntegerField(blank=True, null=True)
    aa_alt = models.CharField(max_length=255, blank=True, null=True)
    hgvs_p = models.CharField(max_length=255, blank=True, null=True)


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

    # class Meta:
    #     unique_together = (( 'contig', 'pos', 'ref'),)


class Variantannotation(models.Model):
    id = models.AutoField(primary_key=True)
    assignment_fk = models.ForeignKey('Variantassignment', models.DO_NOTHING, db_column='assignment_fk')
    source_type = models.CharField(max_length=255)
    source = models.CharField(max_length=255)
    prop = models.CharField(max_length=255)
    value = models.CharField(max_length=255)
    description = models.CharField(max_length=255, blank=True, null=True)


class Variantassignment(models.Model):
    id = models.AutoField(primary_key=True)
    variant_collection_fk = models.ForeignKey('Variantcollection', models.CASCADE, db_column='variant_collection_fk')
    variant_fk = models.ForeignKey(Variant, models.CASCADE, db_column='variant_fk')
    allele_fk = models.ForeignKey(Allele, models.CASCADE, db_column='allele_fk')

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

    def process_allele(self,allele):
        """
        processes an allele from a variant collection and creates a GenotypeSupport if needed
        :param allele:
        :return:
        """

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
    allele = models.ForeignKey(Allele, models.CASCADE)
    effect = models.ForeignKey(Effect, models.CASCADE)
    phenotype = models.ForeignKey(Phenotype, models.CASCADE, "reported", null=True)
    reported_in = models.CharField(max_length=255)
    created_at = models.DateTimeField(auto_now_add=True)
    updated_at = models.DateTimeField(auto_now=True)


class GenotypeSupport(models.Model):
    id = models.AutoField(primary_key=True)
    genotype = models.ForeignKey(Genotype, models.CASCADE, "supportedby")
    status = models.ForeignKey(Term, models.DO_NOTHING)
    reported_allele = models.ForeignKey(ReportedAllele, models.DO_NOTHING, null=True)
    assignment = models.ForeignKey(Variantassignment, models.DO_NOTHING)