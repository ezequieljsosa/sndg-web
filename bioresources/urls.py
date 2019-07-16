from django.urls import path
from .views import ResumableUploadView

from bioresources.views import UploadView
from . import views
from .models import Resource

app_name = 'bioresources'
urlpatterns = [
    path('', views.index, name='index'),
    # path('publications', views.filteredpersonlistview.as_view(), name='publications'),
    path('search/', views.BioSearchView.as_view(), name='search_view'),
    path('structure/<int:pk>', views.structure, name='structure_view'),
    path('publication/<int:pk>', views.publication, name='publication_view'),

    path('tool/<int:pk>', views.tool, name='tool_view'),
    path('reads/<int:pk>', views.reads, name='reads_view'),

    path('sample/<int:pk>', views.sample_view, name='sample_view'),
    path('sample/<str:pk>', views.sample_view, name='sample_view'),

    path('barcode/<int:pk>', views.barcode, name='barcode_view'),

    path('bioproject/<int:pk>', views.bioproject, name='bioproject_view'),
    path('organization/<int:pk>', views.organization, name='organization_view'),
    path('person/<int:pk>', views.person, name='person_view'),
    path('expression/<int:pk>', views.expression, name='expression_view'),
    path('assembly/<int:pk>', views.assembly, name='assembly_view'),
    path('assembly/<str:pk>', views.assembly, name='assembly_view2'),

    path('genome/<int:genome_id>', views.publication, name='genome_view'),

    # forms
    path('ra/new', views.tool, name='ra_new'),
    path('sample/new', views.tool, name='sample_new'),


    # path('upload/', view=views.ExampleView.as_view(), name='home'),
    # path('upload_api/', view=views.NotConcurrentUploaderView.as_view(), name='upload'),

    path('upload/', view=ResumableUploadView.as_view(), name='upload'),
    path('submission/', view=views.submission, name='submission'),
    path('submission/assembly', views.assembly_new, name='assembly_new'),
    path('submission/import', views.import_resource, name='import_resource'),
    # path('upload_api/', view=ResumableUploadView.as_view(),  name='upload_api'),


    path('available_tools/', view=ResumableUploadView.as_view(), name='available_tools'),
    path('stats/', view=ResumableUploadView.as_view(), name='stats'),
    path('about/', view=ResumableUploadView.as_view(), name='about'),

]
