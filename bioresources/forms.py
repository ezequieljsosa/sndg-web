import datetime

from captcha.fields import CaptchaField
from crispy_forms.helper import FormHelper
from crispy_forms.layout import Submit
from django import forms
from django.contrib.auth.forms import UserCreationForm
from django.contrib.auth.models import User
from haystack.forms import SearchForm

from .models import Resource, Assembly


class AssemblyForm(forms.ModelForm):
    # name = forms.CharField(max_length=350,required=True)
    # description = forms.CharField(widget=forms.Textarea,required=False)
    # intraspecific_name = forms.CharField(max_length=250, required=False)
    # species_name = forms.CharField(max_length=200, required=False)
    # level = forms.CharField(max_length=50, required=True)
    # # ncbi_org = forms.CharField(max_length=200, null=True)
    release_date = forms.DateField(required=True, widget=forms.SelectDateWidget(years=range(1990,datetime.datetime.now().year)))
    # # update_date = forms.DateField(null=True)
    # assembly_type = forms.ChoiceField( required=True, choices=(
    #     ("haploid", "haploid"), ("diploid", "diploid"), ("other", "other")))
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




class PublicationForm(forms.Form):
    name = forms.CharField(max_length=100)
    type = forms.CharField(max_length=100)
    start_date = forms.DateField()
    end_date = forms.DateField()

    def clean(self):
        cleaned_data = super(PublicationForm, self).clean()
        start_date = cleaned_data.get("start_date")
        end_date = cleaned_data.get("end_date")

        if start_date and end_date and (start_date > end_date):
            self._errors['start_date'] = self._errors.get('start_date', [])
            self._errors['start_date'].append("Start date must be before end date.")

        return cleaned_data


class SignUpForm(UserCreationForm):
    first_name = forms.CharField(max_length=30, required=False, help_text='Optional.')
    last_name = forms.CharField(max_length=30, required=False, help_text='Optional.')
    email = forms.EmailField(max_length=254, help_text='Required. Inform a valid email address.')

    class Meta:
        model = User
        fields = ('username', 'first_name', 'last_name', 'email', 'password1', 'password2',)


class AllauthSignupForm(forms.Form):
    captcha = CaptchaField()

    def signup(self, request, user):
        """ Required, or else it throws deprecation warnings """
        pass


