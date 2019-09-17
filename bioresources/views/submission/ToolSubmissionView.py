# -*- coding: utf-8 -*-
from __future__ import unicode_literals
from django.utils.translation import gettext_lazy as __
from datetime import datetime

from django.shortcuts import redirect, reverse
from django.shortcuts import render
from django import forms

from crispy_forms.helper import FormHelper
from crispy_forms.layout import Submit
from django.http import HttpResponseRedirect
from django.contrib.auth.decorators import login_required

from bioresources.models.Tool import Tool
from bioresources.models.Resource import Collaboration
from django.db import transaction



class ToolForm(forms.ModelForm):
    release_date = forms.DateField(required=True, widget=forms.SelectDateWidget(years=range(1990, datetime.now().year)))

    class Meta:
        model = Tool
        fields = ["name", "description", "url", "tool_type"]

    def __init__(self, *args, **kwargs):
        super(ToolForm, self).__init__(*args, **kwargs)
        self.helper = FormHelper()
        self.helper.add_input(Submit('submit', 'Submit'))

    def clean(self):
        cleaned_data = super(ToolForm, self).clean()
        qs = Tool.objects.filter(name=cleaned_data["name"])
        if "pk" in self.data:
            if qs.exclude(id=self.data["pk"]).exists():
                self._errors['name'] = self._errors.get('name', [])
                self._errors['name'].append(__("%s already exists") % cleaned_data["name"])
        else:
            if qs.exists():
                self._errors['name'] = self._errors.get('name', [])
                self._errors['name'].append(__("%s already exists") % cleaned_data["name"])


@login_required
def ToolSubmissionView(request):
    if request.method == 'POST':
        if "pk" in request.GET:
            resource = Tool.objects.get(id=request.GET["pk"])
            form = ToolForm(request.POST,instance=resource)
        else:
            form = ToolForm(request.POST)

        if form.is_valid():
            with transaction.atomic():
                tool = form.save()
                if not Collaboration.objects.filter(resource=tool,person=request.user.person).exists():
                    Collaboration.objects.create(resource=tool,person=request.user.person,type=Collaboration.COLLABORATION_TYPES.owner)
            return HttpResponseRedirect(reverse("bioresources:tool_view", args=[tool.id]))
    else:
        if "pk" in request.GET:
            resource = Tool.objects.get(id=request.GET["pk"])
            form = ToolForm(instance=resource)
        else:
            form = ToolForm()

    data = {'form': form}
    if "pk" in request.GET:
        data["pk"] = request.GET["pk"]

    return render(request, 'submission/tool_submission.html', data)
