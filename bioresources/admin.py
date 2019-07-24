# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.contrib import admin
from django.urls import reverse
from django.utils.html import format_html

from .models.Person import Person
from .models.Affiliation import Affiliation
from .models.Organization import Organization
from .models.Resource import Resource
from .models.Publication import Publication
from .models.Identity import Identity
from .models.Expression import Expression
from .models.BioProject import BioProject
from .models.Structure import Structure
from .models.Sample import Sample
from .models.ReadsAchive import ReadsArchive
from .models.Assembly import Assembly
from .models.Barcode import Barcode
from .models.ResourceRelation import ResourceRelation
from .models.RKeyword import RKeyword
from .models.ExternalId import ExternalId
from .models.ResourceProperty import ResourceProperty


admin.site.register(Person)
admin.site.register(Identity)
admin.site.register(RKeyword)
admin.site.register(Affiliation)
admin.site.register(ResourceProperty)
admin.site.register(Sample)
admin.site.register(ReadsArchive)
admin.site.register(ExternalId)




@admin.register(Resource)
class ResourceAdmin(admin.ModelAdmin):
    autocomplete_fields = ["ncbi_tax"]
    search_fields = ["name","description"]
    list_display = ["name", "description", "deprecated", "updated_at"]

    # def link(self, obj):
    #     # https://en.proft.me/2014/10/12/reversing-admin-urls-django/
    #     return format_html('<a href="{url}?pdb_id={{pdb_id}}">Resources</a>',
    #                        pdb_id=obj.id, url=reverse('admin:bioresources_resource_changelist'))


@admin.register(ResourceRelation)
class ResourceRelationAdmin(admin.ModelAdmin):
    autocomplete_fields = ["source", "target"]


@admin.register(Expression)
class ExpressionAdmin(admin.ModelAdmin):
    autocomplete_fields = ["ncbi_tax"]


@admin.register(BioProject)
class BioProjectAdmin(admin.ModelAdmin):
    autocomplete_fields = ["ncbi_tax"]


@admin.register(Assembly)
class AssemblyAdmin(ResourceAdmin):
    autocomplete_fields = ["ncbi_tax"]
    # search_fields = ["name","description"]


@admin.register(Organization)
class OrganizationAdmin(admin.ModelAdmin):
    search_fields = ["name"]


@admin.register(Structure)
class StructureAdmin(admin.ModelAdmin):
    autocomplete_fields = ["ncbi_tax"]


@admin.register(Barcode)
class BarcodeAdmin(admin.ModelAdmin):
    autocomplete_fields = ["ncbi_tax"]




@admin.register(Publication)
class PublicationAdmin(admin.ModelAdmin):
    autocomplete_fields = ["ncbi_tax"]
    list_display = ["name", "description", "links"]
    search_fields = ["name", "description"]

    def links(self, obj):
        # https://en.proft.me/2014/10/12/reversing-admin-urls-django/
        return format_html('<a href="{url}?pdb_id={{pdb_id}}">Resources</a>',
                           pdb_id=obj.id, url=reverse('admin:bioresources_resource_changelist'))
