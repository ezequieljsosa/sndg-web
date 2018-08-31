from django.urls import path

from . import views


app_name = 'bioresources'
urlpatterns = [
    path('', views.index, name='index'),
    # path('publications', views.FilteredPersonListView.as_view(), name='publications'),
    path('search/', views.BioSearchView.as_view(), name='search_view'),
    path('structure/<int:pk>', views.structure, name='structure_view'),
    path('publication/<int:pk>', views.publication, name='pmc_view'),
    path('tool/<int:pk>', views.tool, name='tool_view'),
    path('barcode/<int:pk>', views.barcode, name='barcode_view'),
    path('bioproject/<int:pk>', views.bioproject, name='bioproject_view'),
    path('organization/<int:pk>', views.organization, name='org_view'),
    path('person/<int:pk>', views.person, name='person_view'),
    path('expression/<int:pk>', views.expression, name='gds_view'),
    path('assembly/<int:pk>', views.assembly, name='assembly_view'),

    path('genome/<int:genome_id>', views.publication, name='genome_view'),
]