from django.shortcuts import redirect,reverse
from bioseq.models.Bioentry import Bioentry

def search_redirect(request,acctype,acc):
    if acctype == "locus_tag":
        be = Bioentry.objects.get(accession=acc)
        return redirect(reverse("bioresources:protein_view",args=[be.bioentry_id]))
    else:
        raise Exception("Incorrect redirection")