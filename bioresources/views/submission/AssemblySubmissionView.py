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

from bioresources.models.Assembly import Assembly


class AssemblyForm(forms.ModelForm):
    release_date = forms.DateField(required=True, widget=forms.SelectDateWidget(years=range(1990, datetime.now().year)))

    class Meta:
        model = Assembly
        fields = ["name", "description", "intraspecific_name", "species_name","level"]

    def __init__(self, *args, **kwargs):
        super(AssemblyForm, self).__init__(*args, **kwargs)
        self.helper = FormHelper()
        self.helper.add_input(Submit('submit', 'Submit'))
        if self.instance.id:
            self.fields['name'].widget.attrs['readonly'] = True



    def clean(self):
        cleaned_data = super(self.__class__, self).clean()
        qs = self._meta.model.objects.filter(name=cleaned_data["name"])
        if "pk" in self.data:
            if qs.exclude(id=self.data["pk"]).exists():
                self._errors['name'] = self._errors.get('name', [])
                self._errors['name'].append(__("%s already exists") % cleaned_data["name"])
        else:
            if qs.exists():
                self._errors['name'] = self._errors.get('name', [])
                self._errors['name'].append(__("%s already exists") % cleaned_data["name"])


@login_required
def AssemblySubmissionView(request):
    form_class = AssemblyForm
    model_class = Assembly
    if request.method == 'POST':
        if "pk" in request.GET:
            obj = model_class.objects.get(id=request.GET["pk"])
            form = form_class(request.POST,instance=obj)
        else:
            form = form_class(request.POST)

        if form.is_valid():
            with transaction.atomic():
                obj = form.save()
                if not Collaboration.objects.filter(resource=obj,person=request.user.person).exists():
                    Collaboration.objects.create(resource=obj,person=request.user.person,type=Collaboration.COLLABORATION_TYPES.owner)
            return HttpResponseRedirect(reverse("bioresources:" + obj.type_name()  + "_view", args=[obj.id]))
    else:
        if "pk" in request.GET:
            resource = model_class.objects.get(id=request.GET["pk"])
            form = form_class(instance=resource)
        else:
            form = form_class()

    data = {'form': form}
    if "pk" in request.GET:
        data["pk"] = request.GET["pk"]



    return render(request, 'submission/tool_submission.html', {'form': form})

# from modeltranslation.translator import translator, TranslationOptions
# class NewsTranslationOptions(TranslationOptions):
#     fields = ('title', 'text',)
#
# translator.register(Assembly, NewsTranslationOptions)