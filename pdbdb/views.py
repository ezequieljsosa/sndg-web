from django.shortcuts import render
from django.http import HttpResponse
from django.views.generic import TemplateView, DetailView

from django.db.models import Q
import json

from .models import PDB,ResidueSet,PDBResidueSet,Property,ResidueSetProperty


class StructureView(TemplateView):
    # http://nglviewer.org/ngl/api/class/src/stage/stage.js~Stage.html#instance-method-loadFile
    template_name = "pdbdb/structure_detail.html"

    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)
        context["pdb"] = self.kwargs["pdbid"]

        pdbobj = PDB.objects.prefetch_related("residues").get(code=self.kwargs["pdbid"])


        context["chains"] = [{"name":x} for x in set([r.chain for r in pdbobj.residues.all() if r.chain.strip()]) ]
        context["layers"] = ["water", "hetero"]

        ds = Property.objects.get(name="druggability_score")
        rs = ResidueSet.objects.get(name="FPocketPocket")

        # sq = ResidueSetProperty.objects.select_related(pdbresidue_set)\
        #     .filter(property=ds,value__gte=0.2,pdbresidue_set=OuterRef("id"))

        context["pockets"] = PDBResidueSet.objects.prefetch_related("properties__property","residues__atoms__atom").filter(
            Q(pdb=pdbobj) , Q(residue_set=rs) , Q(properties__property=ds) & Q(properties__value__gte=0.2)).all()
        for p in context["pockets"]:
            p.druggability = [x.value for x in p.properties.all() if x.property == ds][0]
            p.atoms = []
            for r in p.residues.all():
                for a in r.atoms.all():
                    p.atoms.append(a.atom.serial)
        # context["residuesets"] = [{"name": "csa", "residues": range(700, 750)}]
        return context


def structure_raw(request, pdbid):


    pdbobj = PDB.objects.prefetch_related("residues__atoms").get(code=pdbid)
    pdb2 = "\n".join(pdbobj.lines()) + "\n"

    # with open("/tmp/pepe/pepe.ent","w") as h:
    #      h.write(pdb)

    # pdb = open(
    # "/data/databases/pdb/divided/" + pdbid[1:3] + "/pdb" + pdbid + ".ent").read() + "\n"


# "\n" + "\n".join(
#     [x for x in open("/tmp/fpocket/4zu4_out/4zu4_out.pdb").readlines()
#      if x[17:20].strip() == "STP" ])

    return HttpResponse(pdb2  )
