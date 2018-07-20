from django.contrib import admin
from .models import Bioentry, Biosequence, Biodatabase, Ontology, Term, Dbxref, BioentryQualifierValue, Seqfeature, \
    Tool, ToolRun, Location, SeqfeatureDbxref,Taxon

admin.site.register(Biosequence)
admin.site.register(Biodatabase)

admin.site.register(Ontology)
admin.site.register(Term)


admin.site.register(Tool)
admin.site.register(ToolRun)



from easy_select2 import select2_modelform

BioentryForm = select2_modelform(Bioentry, attrs={'width': '250px'})


@admin.register(Bioentry)
class BiosequenceAdmin(admin.ModelAdmin):
    form = BioentryForm
    # fields = ('name','biodatabase','accession','identifier','division','version')
    raw_id_fields = (
        'taxon',
    )
    # readonly_fields = ('taxon_id',)


@admin.register(BioentryQualifierValue)
class BioentryQualifierValueAdmin(admin.ModelAdmin):
    raw_id_fields = (
        'bioentry', "term"
    )


@admin.register(Seqfeature)
class SeqfeatureAdmin(admin.ModelAdmin):
    raw_id_fields = (
        'bioentry', "type_term", "source_term"
    )

@admin.register(Location)
class LocationAdmin(admin.ModelAdmin):
    raw_id_fields = (
        'seqfeature', "dbxref", "term"
    )


@admin.register(Dbxref)
class DbxrefAdmin(admin.ModelAdmin):
    search_fields = ["accession"]

@admin.register(SeqfeatureDbxref)
class SeqfeatureDbxrefAdmin(admin.ModelAdmin):

    autocomplete_fields = ["dbxref"]
    raw_id_fields = (
        'seqfeature',
    )

@admin.register(Taxon)
class TaxonAdmin(admin.ModelAdmin):
    search_fields = ["accession"]
