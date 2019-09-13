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
from .views.models.ReadsArchiveView import reads
from .views.models.ExpressionView import expression
from .views.models.ProteinView import ProteinView
from .views.models.NucleotideView import NucleotideView
from .views.models.TaxView import TaxView

from .views.search.BioSearch import BioSearchView
from .views.search.BioSearchRelated import BioSearchRelatedView
from .views.search import search_redirect
from .views.jobs.BlastView import blast, job_view

from .views.submission.SubmissionNewView import SubmissionNewView
from .views.submission.ToolSubmissionView import ToolSubmissionView
from .views.submission.BioprojectSubmissionView import BioprojectSubmissionView

from .views.submission.SubmissionImportView import SubmissionImportView, SubmitImportView
from .views.submission.SubmissionRelatedView import SubmissionRelatedView, claim_resource, mark_to_relate, \
    claim_identity
from .views.submission import submission_start

from .views.user.UserResourcesView import UserResourcesView

### from .views import ResumableUploadView

app_name = 'bioresources'
urlpatterns = [

    # Main page and resources search
    path('', index.index, name='index'),
    path('blast/', blast, name='available_tools'),
    path('job/<int:jid>', job_view, name='job_view'),
    path('', index.index, name='stats'),

    path('search/', BioSearchView.as_view(), name='search_view'),
    path('search/<str:rtype_src>/<int:rid>/<str:rtype_dst>', BioSearchRelatedView.as_view(), name='search_view'),
    path('search_redirect/<str:acctype>/<str:acc>', search_redirect, name='search_redirect'),

    # Models
    path('bioproject/<int:pk>', bioproject, name='bioproject_view'),
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
    path('tax/<int:pk>', TaxView.as_view(), name='tax_view'),

    # # Jobs
    # path('jobs/', view=ResumableUploadView.as_view(), name='jobs'),
    # path('blast/', view=ResumableUploadView.as_view(), name='blast'),
    #
    #
    # # Upload
    # path('upload/', view=ResumableUploadView.as_view(), name='upload'),

    path('submission/', view=submission_start, name='submission'),
    path('submission/new', view=SubmissionNewView, name='submission_new'),
    path('submission/tool', view=ToolSubmissionView, name='tool_submission'),
    path('submission/bioproject', view=BioprojectSubmissionView, name='bioproject_submission'),
    path('submission/import', view=SubmissionImportView, name='submission_import'),
    path('submission/import/submit', view=SubmitImportView, name='submission_import_submit'),
    path('submission/related', view=SubmissionRelatedView, name='submission_related'),

    path('user/resources', view=UserResourcesView, name='user_resources'),
    path('user/claim_identity/<int:person_id>', view=claim_identity, name='claim_identity'),
    path('relate/<int:src_id>/<int:dst_id>', view=SubmissionRelatedView, name='relate_resources'),
    path('relate/<int:resource_id>', view=claim_resource, name='claim_resource'),
    path('relate/<int:resource_id>/mark', view=mark_to_relate, name='mark_to_relate'),

    #
    #     # # Submission
    #     # path('ra/new', views.tool, name='ra_new'),
    #     # path('sample/new', views.tool, name='sample_new'),
    #     #
    #     #
    #     # # Other
    #     # path('stats/', view=ResumableUploadView.as_view(), name='stats'),
    #
]
# urlpatterns += [path('search/', BioSearchView.as_view(), name='search_view'),   ]
