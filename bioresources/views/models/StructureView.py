# -*- coding: utf-8 -*-
from __future__ import unicode_literals
from django.utils.translation import gettext_lazy as __
from django.shortcuts import redirect, reverse
from django.shortcuts import render

def structure(request, pk):
    return redirect(reverse("pdbdb:structure_view", kwargs={"pdbid": Structure.objects.get(id=pk).name}))

