from datetime import date

from django import forms
from django.http import HttpResponse
from haystack.forms import SearchForm
from haystack.generic_views import SearchView

from .models import Note


class DateRangeSearchForm(SearchForm):
    models = [Note]

    start_date = forms.DateField(required=False)
    end_date = forms.DateField(required=False)

    def search(self):
        # First, store the SearchQuerySet received from other processing.
        sqs = super(DateRangeSearchForm, self)
        #sqs. ("gobierno")
        sqs = sqs.search()
        if not self.is_valid():
            return self.no_query_found()

        # Check to see if a start_date was chosen.
        if self.cleaned_data['start_date']:
            sqs = sqs.filter(pub_date__gte=self.cleaned_data['start_date'])

        # Check to see if an end_date was chosen.
        if self.cleaned_data['end_date']:
            sqs = sqs.filter(pub_date__lte=self.cleaned_data['end_date'])

        return sqs


class MySearchView(SearchView):
    """My custom search view."""

    form_class = DateRangeSearchForm
    #context_object_name = 'page_obj'

    def get_queryset(self):

        queryset = super(MySearchView, self).get_queryset()
        # further filter queryset based on some set of criteria
        return queryset.filter(pub_date__gte=date(2015, 1, 1))

    def get_context_data(self, *args, **kwargs):
        context = super(MySearchView, self).get_context_data(*args, **kwargs)
        # do something

        context["suggestions"] = context["view"].queryset.query._raw.spellcheck["suggestions"]
        return context


def index(request):
    return HttpResponse("Hello, world. You're at the polls index.")
