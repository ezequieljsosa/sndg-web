import os
import sys
import traceback
import warnings
import json


from Bio import BiopythonWarning, BiopythonParserWarning, BiopythonDeprecationWarning, BiopythonExperimentalWarning
from django.core.management.base import BaseCommand, CommandError
from django.db import transaction
from tqdm import tqdm

from pdbdb.models import PDB, Property, PDBResidueSet, ResidueSet, ResidueSetProperty, ResidueSetResidue, \
    AtomResidueSet, Residue, Atom

warnings.simplefilter('ignore', RuntimeWarning)
warnings.simplefilter('ignore', BiopythonWarning)
warnings.simplefilter('ignore', BiopythonParserWarning)
warnings.simplefilter('ignore', BiopythonDeprecationWarning)
warnings.simplefilter('ignore', BiopythonExperimentalWarning)

sys.path.append("/home/eze/workspace/sndg-bio")

from SNDG.Structure.FPocket import FPocket, fpocket_properties_map

pocket_prop_map = {v: k for k, v in fpocket_properties_map.items()}


def iterpdbs(pdbs_dir, pdb_extention=".ent"):
    for index_dir in os.listdir(pdbs_dir):
        if len(index_dir) == 2:
            for x in os.listdir(pdbs_dir + "/" + index_dir):
                if x.endswith(pdb_extention):
                    yield x[3:7], pdbs_dir + "/" + index_dir + "/" + x


class Command(BaseCommand):
    help = 'Loads the pdb files to the database'

    def add_arguments(self, parser):
        parser.add_argument('--pdb', required=False)
        parser.add_argument('--tmp', default="/tmp/fpocket")

    def handle(self, *args, **options):

        tmp = os.path.abspath(options['tmp'])
        if not os.path.exists(tmp):
            os.makedirs(tmp)

        total = PDB.objects.count()

        properties = {name: Property.objects.get_or_create(name=name, description=desc)[0]
                      for name, desc in fpocket_properties_map.items()}

        rsfpocker = ResidueSet.objects.get_or_create(name="FPocketPocket", description="")[0]

        from django.db.models import Max
        class Struct:
            def __init__(self, **entries):
                self.__dict__.update(entries)

        with tqdm(PDB.objects.all(), total=total) as pbar:
            for pdb in pbar:
                pbar.set_description(pdb.code)
                try:
                    pdb = PDB.objects.prefetch_related("residues__atoms").get(id=pdb.id)
                    pdb_path = tmp + "/" + pdb.code + ".pdb"
                    with open(pdb_path, "w") as h:
                        h.write(pdb.text)

                    Residue.objects.filter(pdb=pdb, resname="STP").delete()
                    PDBResidueSet.objects.filter(pdb=pdb,residue_set__name="FPocketPocket").delete()
                    pockets_path = os.path.abspath(tmp) + "/" + pdb.code + "pockets.json"
                    if not os.path.exists(pockets_path):
                        res = FPocket(pdb_path, tmp).hunt_pockets()
                        res.save(pockets_path )
                    with open(pockets_path) as h:
                        res_pockets = json.load(h)
                    rss = []
                    nro_atm = Atom.objects.filter(residue__pdb=pdb).aggregate(Max("serial"))["serial__max"]
                    with transaction.atomic():
                        for pocket in self.tqdm(res_pockets):
                            pocket = Struct(**pocket)
                            rs = PDBResidueSet(name="%i" % pocket.number, pdb=pdb, residue_set=rsfpocker)
                            rss.append(rs)
                            rs.save()

                            res_alpha = {}
                            for stp_line in pocket.as_lines:
                                nro_atm += 1
                                resid = int(stp_line[22:26])
                                if resid in res_alpha:
                                    r = res_alpha[resid]
                                else:
                                    r = Residue(pdb=pdb, chain=stp_line[22:23], resid=resid,
                                            type="",
                                            resname="STP", disordered=1)
                                    r.save()
                                    res_alpha[resid] = r
                                Atom(residue=r, serial=nro_atm, name=stp_line[12:16],
                                     x=float(stp_line[30:38].strip()), y=float(stp_line[38:46].strip()),
                                     z=float(stp_line[46:54].strip()),
                                     occupancy=float(stp_line[54:60].strip()), bfactor=float(stp_line[60:66].strip()),
                                     element="").save()

                            atoms = Atom.objects.select_related("residue").filter(residue__pdb=pdb,
                                serial__in=[int(x) for x in pocket.atoms])
                            residues = set([x.residue for x in atoms])
                            rs_dict = {}
                            for r in residues:
                                rsr = ResidueSetResidue(residue=r, pdbresidue_set=rs)
                                rsr.save()
                                rs_dict[r.id] = rsr

                            for atom in atoms:
                                AtomResidueSet(atom=atom, pdb_set=rs_dict[atom.residue.id]).save()

                            for k, v in pocket.properties.items():
                                prop = pocket_prop_map[k]
                                prop_model = properties[prop]
                                ResidueSetProperty(pdbresidue_set=rs, property=prop_model, value=v).save()

                    # res.delete_dir()


                except IOError as ex:
                    traceback.print_exc()
                    self.stderr.write("error processing pockets from %s: %s" % (pdb.code,str(ex)) )

                except Exception as ex:
                    traceback.print_exc()
                    raise CommandError(ex)
