from django.urls import path

from . import views

app_name = 'biosql'
urlpatterns = [
    path('about', views.AboutView.as_view(), name='index'),
    path('tax/<int:pk>', views.TaxView.as_view(), name='tax_view'),
    path('seq/<int:pk>', views.sequence_view, name='seq_view'),
    path('assembly/<int:pk>', views.assembly_view, name='assembly_view'),
    path('sequence/<int:pk>', views.TaxView.as_view(), name='sequence_view'),
    path('structure/<int:pk>', views.StructureView.as_view(), name='structure_view'),
    path('structure_raw/<str:pdbid>', views.structure_raw, name='structure_raw_view'),


    path('variant/<int:pk>', views.TaxView.as_view(), name='variant_view'),

    path('blast', views.AboutView.as_view(), name='blast_view'),
    path('msa/', views.AboutView.as_view(), name='msa_view'),
    path('primer/', views.AboutView.as_view(), name='primer_view'),



    path('analysis/<int:pk>', views.TaxView.as_view(), name='analysis_view'), # tree / msa / blast



    # path('pathway/<int:pk>', views.TaxView.as_view(), name='tax_view'),

    # path('publications', views.FilteredPersonListView.as_view(), name='publications'),
    # path('search/', views.BioSearchView.as_view(), name='search_view'),


]