from django.urls import path

from . import views


app_name = 'bioresources'
urlpatterns = [
    path('', views.index, name='index'),
    # path('publications', views.FilteredPersonListView.as_view(), name='publications'),
    path('search/', views.BioSearchView.as_view(), name='search_view'),

    path('publication/<int:publication_id>', views.publication, name='publication_view'),
]