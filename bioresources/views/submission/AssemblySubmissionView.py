# -*- coding: utf-8 -*-
from __future__ import unicode_literals
from django.utils.translation import gettext_lazy as __
from datetime import  datetime

from django.shortcuts import redirect, reverse
from django.shortcuts import render
from django import forms

from crispy_forms.helper import FormHelper
from crispy_forms.layout import Submit

from bioresources.models.Assembly import Assembly

class AssemblyForm(forms.ModelForm):
#     # name = forms.CharField(max_length=350,required=True)
#     # description = forms.CharField(widget=forms.Textarea,required=False)
#     # intraspecific_name = forms.CharField(max_length=250, required=False)
#     # species_name = forms.CharField(max_length=200, required=False)
#     # level = forms.CharField(max_length=50, required=True)
#     # # ncbi_org = forms.CharField(max_length=200, null=True)
    release_date = forms.DateField(required=True, widget=forms.SelectDateWidget(years=range(1990,datetime.now().year)))
#     # # update_date = forms.DateField(null=True)
#     # assembly_type = forms.ChoiceField( required=True, choices=(
#     #     ("haploid", "haploid"), ("diploid", "diploid"), ("other", "other")))
    class Meta:
        model = Assembly
        fields = ["name","description","intraspecific_name","species_name","level","release_date","assembly_type",]

    def __init__(self, *args, **kwargs):
        super(AssemblyForm, self).__init__(*args, **kwargs)
        self.helper = FormHelper()
        self.helper.add_input(Submit('submit', 'Submit'))


    def clean(self):
        cleaned_data = super(AssemblyForm, self).clean()
        if Assembly.objects.filter(name=cleaned_data["name"]).exists():
            self._errors['name'] = self._errors.get('name', [])
            self._errors['name'].append("%s already exists" % cleaned_data["name"])

# @login_required
def AssemblySubmissionView(request):
    if request.method == 'POST':
        # create a form instance and populate it with data from the request:
        form = AssemblyForm(request.POST)

        if form.is_valid():
            assembly = form.save()

            return HttpResponseRedirect('/assembly/' + str(assembly.id))


    else:
        form = AssemblyForm()

    return render(request, 'forms/assembly_new.html', {'form': form})