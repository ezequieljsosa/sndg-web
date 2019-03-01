# This is an auto-generated Django model module.
# You'll have to do the following manually to clean this up:
#   * Rearrange models' order
#   * Make sure each model has one field with primary_key=True
#   * Make sure each ForeignKey has `on_delete` set to the desired behavior.
#   * Remove `managed = True` lines if you wish to allow Django to create, modify, and delete the table
# Feel free to rename the models, but don't rename db_table values or field names.
from django.db import models
from collections import defaultdict

from itertools import groupby
from .managers import SeqfeatureManager,BioentryManager

from django.shortcuts import redirect, reverse



"""
ALTER TABLE taxon_name ADD COLUMN id INT AUTO_INCREMENT PRIMARY KEY;
ALTER TABLE bioentry_qualifier_value  ADD COLUMN id INT AUTO_INCREMENT PRIMARY KEY;
ALTER TABLE term_synonym  ADD COLUMN id INT AUTO_INCREMENT PRIMARY KEY;
ALTER TABLE term ADD version INT UNSIGNED DEFAULT 1;

"""


class Biodatabase(models.Model):
    biodatabase_id = models.AutoField(primary_key=True)
    name = models.CharField(unique=True, max_length=128)
    authority = models.CharField(max_length=128, blank=True, null=True)
    description = models.TextField(blank=True, null=True)

    class Meta:
        managed = True
        db_table = 'biodatabase'

    def get_absolute_url(self):
        return reverse('biosql:assembly_view', args=[str(self.biodatabase_id)])


class Bioentry(models.Model):
    bioentry_id = models.AutoField(primary_key=True)
    biodatabase = models.ForeignKey(Biodatabase, models.CASCADE, "entries")
    taxon = models.ForeignKey('Taxon', models.DO_NOTHING, blank=True, null=True)
    name = models.CharField(max_length=40)
    accession = models.CharField(max_length=128)
    identifier = models.CharField(max_length=40, blank=True, null=True)
    division = models.CharField(max_length=6, blank=True, null=True)
    description = models.TextField(blank=True, null=True)
    version = models.PositiveSmallIntegerField(default=1, null=True)

    index_updated = models.BooleanField(default=False)

    objects = BioentryManager()

    class Meta:
        managed = True
        db_table = 'bioentry'
        unique_together = (('accession', 'biodatabase', 'version'), ('identifier', 'biodatabase'),)

    def get_absolute_url(self):
        return reverse('biosql:seq_view', args=[str(self.bioentry_id)])

    def groupedFeatures(self):
        group = defaultdict(lambda: [])
        for f in self.features.all():
            group[f.type_term.identifier].append(f.locations.first())

        return dict(group)

    def feature_counts(self):
        data = defaultdict(lambda: 0)
        for f in self.features.all():
            data[f.type_term.name] += 1
        return dict(data)

    def genes(self):
        beg = Biodatabase.objects.get(name=self.biodatabase.name.replace("_prots", ""))
        feature = Seqfeature.objects.seqfeature_from_locus_tag(beg.biodatabase_id, self.accession)
        feature = list(feature)[0]
        return [x.value for x in
                feature.qualifiers.filter(
                    term__name__in=["gene_symbol", "old_locus_tag", "protein_id", "Alias", "gene"])]

    def product_description(self):
        qs = self.qualifiers.filter(term__name="product")
        return qs.value if qs.exists() else None

    def molecular_weight(self):
        qs = self.qualifiers.filter(term__name="molecular_weight")
        return qs.value if qs.exists() else None



    def go_terms(self, database):
        return self.qualifiers.filter(term__dbxrefs__dbxref__dbname="go",
                                      term__dbxrefs__dbxref__accession=database)
        # term__dbxrefs__dbxref__accession="goslim_generic",

    def biological_process(self):
        return self.go_terms("biological_process")

    def molecular_function(self):
        return self.go_terms("molecular_function")

    def cellular_component(self):
        return self.go_terms("cellular_component")

    def ftype(self):
        return "protein"




class BioentryDbxref(models.Model):
    bioentry_dbxref_id = models.AutoField(primary_key=True)
    bioentry = models.ForeignKey(Bioentry, models.CASCADE, related_name="dbxrefs")
    dbxref = models.ForeignKey('Dbxref', models.DO_NOTHING)
    rank = models.SmallIntegerField(default=1, null=True)

    class Meta:
        managed = True
        db_table = 'bioentry_dbxref'
        unique_together = (('bioentry', 'dbxref', 'rank'),)


class BioentryPath(models.Model):
    bioentry_path_id = models.AutoField(primary_key=True)
    object_bioentry = models.ForeignKey(Bioentry, models.DO_NOTHING, related_name="object_bioentry_path")
    subject_bioentry = models.ForeignKey(Bioentry, models.DO_NOTHING, related_name="subject_bioentry_path")
    term = models.ForeignKey('Term', models.DO_NOTHING)
    distance = models.PositiveIntegerField(blank=True, null=True)

    class Meta:
        managed = True
        db_table = 'bioentry_path'
        unique_together = (('object_bioentry', 'subject_bioentry', 'term', 'distance'),)


class BioentryQualifierValue(models.Model):
    bioentry_qualifiervalue_id = models.AutoField(primary_key=True)
    bioentry = models.ForeignKey(Bioentry, models.CASCADE, related_name="qualifiers")
    term = models.ForeignKey('Term', models.DO_NOTHING)
    value = models.TextField(blank=True, null=True)
    rank = models.IntegerField(default=1, null=True)

    class Meta:
        managed = True
        db_table = 'bioentry_qualifier_value'
        unique_together = (('bioentry', 'term', 'rank'),)


class BioentryReference(models.Model):
    bioentry_relationship_id = models.AutoField(primary_key=True)
    bioentry = models.ForeignKey(Bioentry, models.CASCADE)
    reference = models.ForeignKey('Reference', models.DO_NOTHING)
    start_pos = models.IntegerField(blank=True, null=True)
    end_pos = models.IntegerField(blank=True, null=True)
    rank = models.SmallIntegerField(default=1, null=True)

    class Meta:
        managed = True
        db_table = 'bioentry_reference'
        unique_together = (('bioentry', 'reference', 'rank'),)


class BioentryRelationship(models.Model):
    bioentry_relationship_id = models.AutoField(primary_key=True)
    object_bioentry = models.ForeignKey(Bioentry, models.DO_NOTHING, related_name="object_bioentry_relationship")
    subject_bioentry = models.ForeignKey(Bioentry, models.DO_NOTHING, related_name="subject_bioentry_relationship")
    term = models.ForeignKey('Term', models.DO_NOTHING)
    rank = models.IntegerField(default=1, null=True)

    class Meta:
        managed = True
        db_table = 'bioentry_relationship'
        unique_together = (('object_bioentry', 'subject_bioentry', 'term'),)


class Biosequence(models.Model):
    bioentry = models.OneToOneField(Bioentry, models.CASCADE, primary_key=True, related_name="seq")
    version = models.SmallIntegerField(blank=True, default=1, null=True)
    length = models.IntegerField(blank=True, null=True)
    alphabet = models.CharField(max_length=10, blank=True, null=True)
    seq = models.TextField(blank=True, null=True)

    class Meta:
        managed = True
        db_table = 'biosequence'


class Comment(models.Model):
    comment_id = models.AutoField(primary_key=True)
    bioentry = models.ForeignKey(Bioentry, models.CASCADE)
    comment_text = models.TextField()
    rank = models.SmallIntegerField(default=1, null=True)

    class Meta:
        managed = True
        db_table = 'comment'
        unique_together = (('bioentry', 'rank'),)


class Dbxref(models.Model):
    dbxref_id = models.AutoField(primary_key=True)
    dbname = models.CharField(max_length=40)
    accession = models.CharField(max_length=128)
    version = models.PositiveSmallIntegerField(default=1, null=True)

    class Meta:
        managed = True
        db_table = 'dbxref'
        unique_together = (('accession', 'dbname', 'version'),)

    def __str__(self):
        return self.dbname + ":" + self.accession + "." + str(self.version)


class DbxrefQualifierValue(models.Model):
    dbxref = models.OneToOneField(Dbxref, models.DO_NOTHING, primary_key=True)
    term = models.ForeignKey('Term', models.DO_NOTHING)
    rank = models.SmallIntegerField(default=1, null=True)
    value = models.TextField(blank=True, null=True)

    class Meta:
        managed = True
        db_table = 'dbxref_qualifier_value'
        unique_together = (('dbxref', 'term', 'rank'),)


class Location(models.Model):
    location_id = models.AutoField(primary_key=True)
    seqfeature = models.ForeignKey('Seqfeature', models.CASCADE, related_name="locations")
    dbxref = models.ForeignKey(Dbxref, models.DO_NOTHING, blank=True, null=True)
    term = models.ForeignKey('Term', models.DO_NOTHING, blank=True, null=True)
    start_pos = models.IntegerField(blank=True, null=True)
    end_pos = models.IntegerField(blank=True, null=True)
    strand = models.IntegerField()
    rank = models.SmallIntegerField(default=1, null=True)

    class Meta:
        managed = True
        db_table = 'location'
        unique_together = (('seqfeature', 'rank'),)

    def __str__(self):
        return "%i-%i(%s)" % (self.start_pos, self.end_pos, "+" if self.strand == 1 else "-")


class LocationQualifierValue(models.Model):
    location = models.OneToOneField(Location, models.DO_NOTHING, primary_key=True)
    term = models.ForeignKey('Term', models.DO_NOTHING)
    value = models.CharField(max_length=255)
    int_value = models.IntegerField(blank=True, null=True)

    class Meta:
        managed = True
        db_table = 'location_qualifier_value'
        unique_together = (('location', 'term'),)


class Ontology(models.Model):
    GO = "Gene Ontology"
    SO = "Sequence Ontology"
    EC = "Enzyme Commission number"
    SFS = "SeqFeature Sources"
    GRAPH = "Graph"

    ontology_id = models.AutoField(primary_key=True)
    name = models.CharField(unique=True, max_length=32)
    definition = models.TextField(blank=True, null=True)

    class Meta:
        managed = True
        db_table = 'ontology'

    def __str__(self):
        return "%s  (%i)" % (self.name, self.ontology_id)


class Reference(models.Model):
    reference_id = models.AutoField(primary_key=True)
    dbxref = models.OneToOneField(Dbxref, models.DO_NOTHING, unique=True, blank=True, null=True)
    location = models.TextField()
    title = models.TextField(blank=True, null=True)
    authors = models.TextField(blank=True, null=True)
    crc = models.CharField(unique=True, max_length=32, blank=True, null=True)

    class Meta:
        managed = True
        db_table = 'reference'


class Tool(models.Model):
    name = models.CharField(max_length=120, blank=False)
    description = models.TextField(default="")
    version = models.CharField(max_length=64, blank=True, default=1, null=True)
    url = models.URLField()

    def rtype(self):
        return "tool"


class ToolRun(models.Model):
    tool = models.ForeignKey(Tool, models.DO_NOTHING, related_name="runs")
    parameters = models.TextField(blank=False)
    result_url = models.CharField(max_length=255, blank=True, null=True)
    result_path = models.FilePathField(path="/data/runs", null=True, recursive=True, allow_folders=True)
    created_at = models.DateTimeField(auto_now_add=True)
    executed_at = models.DateTimeField(auto_now=True)


# class ToolRunResult(models.Model):
#     toolrunresult_id = models.AutoField(primary_key=True)
#     run = models.ForeignKey(ToolRun, models.DO_NOTHING, related_name="results")
#     description = models.CharField(max_length=255, blank=True, null=True)
#
#
# class ToolRunResultAttrs(models.Model):
#     runResult = models.ForeignKey(ToolRunResult, models.DO_NOTHING, related_name="attrs")
#     name = models.CharField(max_length=255, blank=True, null=True)

# class AligmentSimple(ToolRunResult):
#     query = models.TextField(blank=False)
#     query_seq = models.TextField(blank=False)
#     query_start = models.IntegerField()
#     query_end = models.IntegerField()
#     query_strand = models.IntegerField(default=1)


class Seqfeature(models.Model):

    ENTRY_TYPES = ["CDS","rRNA","tRNA","regulatory","ncRNA","mRNA","repeat"]

    seqfeature_id = models.AutoField(primary_key=True)
    bioentry = models.ForeignKey(Bioentry, models.CASCADE, related_name="features")
    type_term = models.ForeignKey('Term', models.DO_NOTHING, related_name="features_of_type")
    source_term = models.ForeignKey('Term', models.DO_NOTHING, related_name="source_of")
    display_name = models.CharField(max_length=64, blank=True, null=True)
    rank = models.PositiveSmallIntegerField(default=1, null=True)

    index_updated = models.BooleanField(default=False)

    def qualifiers_dict(self):
        return {x.term.name: x.value for x in self.qualifiers.all()}

    objects = SeqfeatureManager()

    class Meta:
        managed = True
        db_table = 'seqfeature'
        unique_together = (('bioentry', 'type_term', 'source_term', 'rank'),)

    def locus_tag(self):
        return self.qualifiers.get(term__name='locus_tag').value

    def genes(self):
        return [x.value for x in
                self.qualifiers.filter(term_name__in=["gene_symbol", "old_locus_tag", "protein_id", "Alias", "gene"])]

    def description(self):
        qs = self.qualifiers.filter(term__name='product')
        if qs.exists():
            return qs.get().value
        return self.type_term.name

    def length(self):
        return sum([abs(x.end_pos - x.start_pos) for x in self.locations])


class SeqfeatureDbxref(models.Model):
    seqfeature_dbxref_id = models.AutoField(primary_key=True)
    seqfeature = models.ForeignKey(Seqfeature, models.DO_NOTHING, related_name="dbxrefs")
    dbxref = models.ForeignKey(Dbxref, models.DO_NOTHING)
    rank = models.SmallIntegerField(default=1, null=True)

    class Meta:
        managed = True
        db_table = 'seqfeature_dbxref'
        unique_together = (('seqfeature', 'dbxref'),)


class SeqfeaturePath(models.Model):
    seqfeature_path_id = models.AutoField(primary_key=True)
    object_seqfeature = models.ForeignKey(Seqfeature, models.DO_NOTHING, related_name="object_paths")
    subject_seqfeature = models.ForeignKey(Seqfeature, models.DO_NOTHING, related_name="subject_paths")
    term = models.ForeignKey('Term', models.DO_NOTHING)
    distance = models.PositiveIntegerField(blank=True, null=True)

    class Meta:
        managed = True
        db_table = 'seqfeature_path'
        unique_together = (('object_seqfeature', 'subject_seqfeature', 'term', 'distance'),)


class SeqfeatureQualifierValue(models.Model):
    seqfeature_qualifiervalue_id = models.AutoField(primary_key=True)
    seqfeature = models.ForeignKey(Seqfeature, models.DO_NOTHING, related_name="qualifiers")
    term = models.ForeignKey('Term', models.DO_NOTHING)
    rank = models.SmallIntegerField(default=1, null=True)
    value = models.TextField()

    class Meta:
        managed = True
        db_table = 'seqfeature_qualifier_value'
        unique_together = (('seqfeature', 'term', 'rank'),)
        indexes = [
            models.Index(fields=['term']),
        ]

    def __str__(self):
        return str(self.term.name) + ":" + self.value


class SeqfeatureRelationship(models.Model):
    seqfeature_relationship_id = models.AutoField(primary_key=True)
    object_seqfeature = models.ForeignKey(Seqfeature, models.DO_NOTHING, related_name="object_relationships")
    subject_seqfeature = models.ForeignKey(Seqfeature, models.DO_NOTHING, related_name="subject_relationships")
    term = models.ForeignKey('Term', models.DO_NOTHING)
    rank = models.IntegerField(default=1, null=True)

    class Meta:
        managed = True
        db_table = 'seqfeature_relationship'
        unique_together = (('object_seqfeature', 'subject_seqfeature', 'term'),)


class Taxon(models.Model):
    SCIENTIFIC_NAME = "scientific name"

    taxon_id = models.AutoField(primary_key=True)
    ncbi_taxon_id = models.IntegerField(unique=True, blank=True, null=True)
    parent_taxon = models.ForeignKey("self", models.DO_NOTHING, related_name="children", null=True)
    node_rank = models.CharField(max_length=32, blank=True, null=True)
    genetic_code = models.PositiveIntegerField(blank=True, null=True)
    mito_genetic_code = models.PositiveIntegerField(blank=True, null=True)
    left_value = models.PositiveIntegerField(unique=True, blank=True, null=True)
    right_value = models.PositiveIntegerField(unique=True, blank=True, null=True)

    def scientific_name(self):
        return self.names.filter(name_class=Taxon.SCIENTIFIC_NAME).first().name

    def other_names(self):
        return {x.name_class: x.name for x in self.names.filter(name_class__in=["synonym", "common name"])}

    def parents(self):
        parents = []
        curr = self.parent_taxon
        while curr.taxon_id != 1:
            parents.append(curr)
            curr = curr.parent_taxon

        return parents

    class Meta:
        managed = True
        db_table = 'taxon'

    def __str__(self):
        return "(" + str(self.ncbi_taxon_id) + ")" + self.scientific_name()


class TaxonName(models.Model):
    id = models.AutoField(primary_key=True)
    taxon = models.ForeignKey(Taxon, models.DO_NOTHING, related_name="names")
    name = models.CharField(max_length=255)
    name_class = models.CharField(max_length=32)

    class Meta:
        managed = True
        db_table = 'taxon_name'
        unique_together = (('taxon', 'name', 'name_class'),)


class Term(models.Model):
    term_id = models.AutoField(primary_key=True)
    name = models.TextField(blank=True, null=True)
    definition = models.TextField(blank=True, null=True)
    identifier = models.CharField(unique=True, max_length=255, blank=True, null=True)
    is_obsolete = models.CharField(max_length=1, blank=True, null=True)
    ontology = models.ForeignKey(Ontology, models.DO_NOTHING)
    version = models.PositiveSmallIntegerField(default=1, null=True)

    class Meta:
        managed = True
        db_table = 'term'
        # unique_together = (('name', 'ontology', 'is_obsolete'),)

    def __str__(self):
        return "%s - %s (%i)" % (self.identifier, self.name, self.ontology.ontology_id)


class TermDbxref(models.Model):
    term_dbxref_id = models.AutoField(primary_key=True)
    term = models.ForeignKey(Term, models.CASCADE, related_name="dbxrefs")
    dbxref = models.ForeignKey(Dbxref, models.DO_NOTHING)
    rank = models.SmallIntegerField(default=1, null=True)

    class Meta:
        managed = True
        db_table = 'term_dbxref'
        unique_together = (('term', 'dbxref'),)


class TermPath(models.Model):
    term_path_id = models.AutoField(primary_key=True)
    subject_term = models.ForeignKey(Term, models.DO_NOTHING, related_name="subject_termpaths")
    predicate_term = models.ForeignKey(Term, models.DO_NOTHING, related_name="predicate_termpaths")
    object_term = models.ForeignKey(Term, models.DO_NOTHING, related_name="object_termpaths")
    ontology = models.ForeignKey(Ontology, models.DO_NOTHING)
    distance = models.PositiveIntegerField(blank=True, null=True)

    class Meta:
        managed = True
        db_table = 'term_path'
        unique_together = (('subject_term', 'predicate_term', 'object_term', 'ontology', 'distance'),)


class TermRelationship(models.Model):
    term_relationship_id = models.AutoField(primary_key=True)
    subject_term = models.ForeignKey(Term, models.DO_NOTHING, related_name="subject_termrelationships")  # parent
    predicate_term = models.ForeignKey(Term, models.DO_NOTHING, related_name="predicate_termrelationships")
    object_term = models.ForeignKey(Term, models.DO_NOTHING, related_name="object_termrelationships")  # child
    ontology = models.ForeignKey(Ontology, models.DO_NOTHING)

    class Meta:
        managed = True
        db_table = 'term_relationship'
        unique_together = (('subject_term', 'predicate_term', 'object_term', 'ontology'),)


class TermRelationshipTerm(models.Model):
    term_relationship = models.OneToOneField(TermRelationship, models.DO_NOTHING, primary_key=True)
    term = models.OneToOneField(Term, models.DO_NOTHING, unique=True)

    class Meta:
        managed = True
        db_table = 'term_relationship_term'


class TermSynonym(models.Model):
    term_synonym_id = models.AutoField(primary_key=True)
    synonym = models.CharField(max_length=255)
    term = models.ForeignKey(Term, models.CASCADE, related_name="synonyms")

    class Meta:
        managed = True
        db_table = 'term_synonym'
        unique_together = (('term', 'synonym'),)


class TaxIdx(models.Model):
    tax = models.OneToOneField(Taxon, models.CASCADE, primary_key=True, db_column="tax_id", related_name="keywords")
    text = models.TextField()

    class Meta:
        managed = True
        db_table = 'tax_idx'


class TermIdx(models.Model):
    term = models.OneToOneField(Term, models.CASCADE, primary_key=True, db_column="term_id", related_name="keywords")
    text = models.TextField()

    class Meta:
        managed = True
        db_table = 'term_idx'
