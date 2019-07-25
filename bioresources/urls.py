from django.urls import path

from .views import index
from .views.models.AssemblyView import assembly_view
from .views.models.BioProjectView import bioproject
from .views.models.OrganizationView import organization
from .views.models.PersonView import person
from .views.models.PublicationView import publication
from .views.models.SampleView import sample_view
from .views.models.StructureView import structure
from .views.models.ToolView import tool
from .views.models.BarcodeView import barcode
from .views.models.ReadsAchiveView import reads
from .views.models.ExpressionView import expression
from .views.models.ProteinView import ProteinView
from .views.models.NucleotideView import NucleotideView
from .views.models.TaxView import TaxView

from .views.search.BioSearch import BioSearchView
from .views.search import search_redirect
from .views.jobs.BlastView import blast,blast_result



# from .views import ResumableUploadView

app_name = 'bioresources'
urlpatterns = [

    # Main page and resources search
    path('', index.index, name='index'),
    path('blast/', blast, name='available_tools'),
    path('job/<int:jid>', blast_result, name='job'),
    path('', index.index, name='stats'),

    path('search/', BioSearchView.as_view(), name='search_view'),
    path('search_redirect/<str:acctype>/<str:acc>',  search_redirect , name='search_redirect'),

    # Models
    path('bioproject/<int:pk>',  bioproject, name='bioproject_view'),
    path('organization/<int:pk>', organization, name='organization_view'),
    path('person/<int:pk>', person, name='person_view'),
    path('expression/<int:pk>', expression, name='expression_view'),
    path('assembly/<int:pk>', assembly_view, name='assembly_view'),
    path('genome/<int:genome_id>', publication, name='genome_view'),
    path('sample/<int:pk>', sample_view, name='sample_view'),
    path('structure/<int:pk>', structure, name='structure_view'),
    path('publication/<int:pk>', publication, name='publication_view'),
    path('barcode/<int:pk>', barcode, name='barcode_view'),
    path('tool/<int:pk>', tool, name='tool_view'),
    path('reads/<int:pk>', reads, name='reads_view'),
    path('nucleotide/<int:pk>', NucleotideView, name='nucleotide_view'),
    path('protein/<int:pk>', ProteinView, name='protein_view'),
    path('tax/<int:pk>',  TaxView.as_view(), name='tax_view'),




    # path('sample/<str:pk>', sample_view, name='sample_view2'),
    # path('assembly/<str:pk>', assembly, name='assembly_view2'),

    # # Jobs
    # path('jobs/', view=ResumableUploadView.as_view(), name='jobs'),
    # path('blast/', view=ResumableUploadView.as_view(), name='blast'),
    #
    #
    # # Upload
    # path('upload/', view=ResumableUploadView.as_view(), name='upload'),
    # path('submission/', view=views.submission, name='submission'),
    # path('submission/assembly', views.assembly_new, name='assembly_new'),
    # path('submission/import', views.import_resource, name='import_resource'),
    #
    #
    # # Submission
    # path('ra/new', views.tool, name='ra_new'),
    # path('sample/new', views.tool, name='sample_new'),
    #
    #
    # # Other
    # path('stats/', view=ResumableUploadView.as_view(), name='stats'),


]
urlpatterns += [path('search/', BioSearchView.as_view(), name='search_view'),
                ]