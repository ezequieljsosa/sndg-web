# -*- coding: utf-8 -*-
from __future__ import unicode_literals
from django.utils.html import format_html
from django.contrib import admin
from django.urls import reverse

from .models import Person, Organization, Resource, Publication, Identity, Expression, BioProject, Structure, Assembly, \
    Barcode, ResourceRelation

admin.site.register(Person)

admin.site.register(Identity)


@admin.register(Resource)
class ResourceAdmin(admin.ModelAdmin):
    autocomplete_fields = ["ncbi_tax"]
    search_fields = ["accession"]


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
class AssemblyAdmin(admin.ModelAdmin):
    autocomplete_fields = ["ncbi_tax"]


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
