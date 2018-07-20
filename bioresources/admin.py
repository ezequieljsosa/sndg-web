# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.contrib import admin

from .models import Person,Organization,Resource,Publication,Identity,Expression,BioProject,Structure,Assembly

admin.site.register(Person)

admin.site.register(Resource)
admin.site.register(Publication)
admin.site.register(Identity)

admin.site.register(Expression)
admin.site.register(BioProject)
admin.site.register(Assembly)


@admin.register(Organization)
class OrganizationAdmin(admin.ModelAdmin):

    search_fields = ["name"]

@admin.register(Structure)
class StructureAdmin(admin.ModelAdmin):

    autocomplete_fields = ["ncbi_tax"]
