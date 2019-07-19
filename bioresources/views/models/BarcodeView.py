# -*- coding: utf-8 -*-
from __future__ import unicode_literals
from django.utils.translation import gettext_lazy as __
from django.shortcuts import redirect, reverse
from django.shortcuts import render


def barcode(request, pk):
    barcode = Barcode.objects.get(id=pk)
    return render(request, 'resources/barcode.html', {
        "barcode": barcode,
        "sidebarleft": 1, })
