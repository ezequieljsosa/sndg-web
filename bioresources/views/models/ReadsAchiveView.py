# -*- coding: utf-8 -*-
from __future__ import unicode_literals
from django.utils.translation import gettext_lazy as __
from django.shortcuts import redirect, reverse
from django.shortcuts import render

def reads(request, pk):
    tool = Tool.objects.get(id=pk)
    return render(request, 'resources/reads.html', {
        "tool": tool,
        "sidebarleft": 1, })