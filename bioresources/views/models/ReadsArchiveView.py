# -*- coding: utf-8 -*-
from __future__ import unicode_literals
from django.utils.translation import gettext_lazy as __
from django.shortcuts import redirect, reverse
from django.shortcuts import render

from bioresources.models.ReadsArchive import ReadsArchive

def reads(request, pk):
    sra = ReadsArchive.objects.get(id=pk)
    return render(request, 'resources/readsarchive.html', {
        "readsarchive": sra,
        "sidebarleft": 1, })