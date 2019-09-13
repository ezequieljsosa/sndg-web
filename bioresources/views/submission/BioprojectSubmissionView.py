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
        if Tool.objects.filter(name=cleaned_data["name"]).exists():
            self._errors['name'] = self._errors.get('name', [])
            self._errors['name'].append(__("%s already exists") % cleaned_data["name"])


@login_required
def ToolSubmissionView(request):
    if request.method == 'POST':
        form = ToolForm(request.POST)

        if form.is_valid():
            tool = form.save()
            from bioresources.graph import Tool as gTool
            node = gTool(rid=tool.id, title=tool.name, tool_type=tool.tool_type)
            node.save()
            return HttpResponseRedirect('/tool/' + str(tool.id))
    else:
        form = ToolForm()

    return render(request, 'submission/tool_submission.html', {'form': form})
