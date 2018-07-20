
from django.db import models

class Allele(models.Model):
    id = models.AutoField(primary_key=True)
    variant_fk = models.ForeignKey('Variant', models.DO_NOTHING, db_column='variant_fk')
    alt = models.CharField(max_length=255)
    main_effect_fk = models.ForeignKey('Effect', models.DO_NOTHING, db_column='main_effect_fk', blank=True, null=True)

    class Meta:
        unique_together = (('variant_fk', 'alt'),)

class Effect(models.Model):
    id = models.AutoField(primary_key=True)
    allele_fk = models.ForeignKey(Allele, models.DO_NOTHING, db_column='allele_fk')
    transcript = models.CharField(max_length=255, blank=True, null=True)
    variant_type = models.CharField(max_length=255)
    reported = models.IntegerField()
    predicted_impact = models.CharField(max_length=255, blank=True, null=True)
    aa_ref = models.CharField(max_length=255, blank=True, null=True)
    aa_pos = models.IntegerField(blank=True, null=True)
    aa_alt = models.CharField(max_length=255, blank=True, null=True)



class Variant(models.Model):
    id = models.AutoField(primary_key=True)
    pos = models.IntegerField()
    gene = models.CharField(max_length=255, blank=True, null=True)
    gene_pos = models.IntegerField(blank=True, null=True)
    contig = models.CharField(max_length=255)
    description = models.CharField(max_length=255, blank=True, null=True)
    ref_organism = models.CharField(max_length=255)
    ref = models.CharField(max_length=255)
    modified_date = models.DateTimeField()

    class Meta:
        unique_together = (('ref_organism', 'contig', 'pos', 'ref'),)


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
    ref_organism = models.CharField(max_length=255)
    sample = models.CharField(max_length=255)
    description = models.CharField(max_length=255, blank=True, null=True)
    modified_date = models.DateTimeField()

    class Meta:
        unique_together = (('ref_organism', 'sample'),)
