from django import forms
from django.contrib.auth.forms import UserCreationForm
from django.contrib.auth.models import User
from captcha.fields import CaptchaField
from haystack.forms import SearchForm

from .models import Resource


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


class BioSearchForm(SearchForm):
    models = [Resource]

    type = forms.TextInput()
    authors = forms.TextInput()
    affiliations = forms.TextInput()
    taxon = forms.TextInput()

    # start_date = forms.DateField(required=False)
    # end_date = forms.DateField(required=False)

    def search(self):
        # First, store the SearchQuerySet received from other processing.
        sqs = super(BioSearchForm, self)

        if not self.is_valid():
            return self.no_query_found()

        if not self.cleaned_data.get('q'):
            return self.no_query_found()
        from haystack.inputs import AutoQuery, Exact

        params = {"content": AutoQuery(self.cleaned_data['q']),
                  "type": Exact(self.data["type"])}
        for x in ["authors", "affiliations", "taxon"] + Resource.facet_dict.get(self.data["type" ], []):
            if x in self.data and self.data[x]:
                params[x] = self.data [x]
        sqs = self.searchqueryset.filter(**params)

        if self.load_all:
            sqs = sqs.load_all()

        # # Check to see if a start_date was chosen.
        # if self.cleaned_data['start_date']:
        #     sqs = sqs.filter(pub_date__gte=self.cleaned_data['start_date'])
        #
        # # Check to see if an end_date was chosen.
        # if self.cleaned_data['end_date']:
        #     sqs = sqs.filter(pub_date__lte=self.cleaned_data['end_date'])

        return sqs
