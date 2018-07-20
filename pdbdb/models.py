from django.db import models
from datetime import datetime


class CustomBinaryCharField(models.CharField):
    def db_type(self, connection):
        return super(CustomBinaryCharField, self).db_type(connection) + ' binary'


class PDB(models.Model):
    id = models.AutoField(primary_key=True)
    code = models.CharField(max_length=100)
    resolution = models.FloatField(default=20)
    experiment = models.CharField(max_length=100, null=True)
    tax = models.CharField(max_length=255, null=True)
    deprecated = models.BooleanField(default=False)
    date = models.DateTimeField(default=datetime.now)

    class Meta:
        unique_together = (('code'),)

    def lines(self):
        # http://www.wwpdb.org/documentation/file-format-content/format33/sect9.html#ATOM
        # http://www.wwpdb.org/documentation/file-format-content/format33/sect9.html#HETATM
        lines = []
        for r in self.residues.all():
            for a in r.atoms.all():
                if r.type == "R":
                    line = "ATOM  "  # 1 -  6        Record name   "ATOM  "
                else:
                    line = "HETATM"  # 1 -  6        Record name   "ATOM  "

                line += str(a.serial).rjust(5)  # 7 - 11        Integer       serial       Atom  serial number.
                if r.resname == "STP":
                    line += " "
                    line += a.name.ljust(4)  # 13 - 16        Atom          name         Atom name.
                else:

                    line += "  "
                    line += a.name.ljust(3)  # 13 - 16        Atom          name         Atom name.
                line += a.altLoc  # 17             Character     altLoc       Alternate location indicator.
                line += r.resname.rjust(3)  # 18 - 20        Residue name  resName      Residue name.
                line += " "
                line += r.chain  # 22             Character     chainID      Chain identifier.
                line += str(r.resid).rjust(4)  # 23 - 26        Integer       resSeq       Residue sequence number.
                line += r.icode.rjust(
                    4)  # 27             AChar         iCode        Code for insertion of residues.
                line += ("%.3f" % a.x).rjust(
                    8)  # 31 - 38        Real(8.3)     x            Orthogonal coordinates for X in Angstroms.
                line += ("%.3f" % a.y).rjust(
                    8)  # 39 - 46        Real(8.3)     y            Orthogonal coordinates for Y in Angstroms.
                line += ("%.3f" % a.z).rjust(
                    8)  # 47 - 54        Real(8.3)     z            Orthogonal coordinates for Z in Angstroms.
                line += ("%.2f" % a.occupancy).rjust(6)  # 55 - 60        Real(6.2)     occupancy    Occupancy.
                line += ("%.2f" % a.bfactor).rjust(
                    6)  # 61 - 66        Real(6.2)     tempFactor   Temperature  factor.
                line += "".rjust(10)
                line += a.element.rjust(2)  # 77 - 78        LString(2)    element      Element symbol, right-justified.
                line += "  "  # 79 - 80        LString(2)    charge       Charge  on the atom.

                lines.append(line)
        return lines


class Residue(models.Model):
    HETATOM = 'H'
    ATOM = 'A'

    id = models.AutoField(primary_key=True)
    pdb = models.ForeignKey(PDB, related_name='residues',
                            db_column="pdb_id", on_delete=models.DO_NOTHING)
    chain = CustomBinaryCharField(max_length=20)
    resname = models.CharField(max_length=4)
    resid = models.IntegerField()
    icode = models.CharField(max_length=2, default="")

    type = models.CharField(max_length=10)
    disordered = models.BooleanField(default=False)

    modelable = models.BooleanField(default=True)
    seq_order = models.IntegerField(null=True)

    class Meta:
        unique_together = (('pdb', "chain", "resid", "icode", "type"),)


class Atom(models.Model):
    id = models.AutoField(primary_key=True)

    serial = models.IntegerField()
    name = models.CharField(max_length=50)
    altLoc = models.CharField(max_length=1,default=" ")

    residue = models.ForeignKey(Residue, related_name='atoms',
                                db_column="residue_id", on_delete=models.DO_NOTHING)
    x = models.FloatField()
    y = models.FloatField()
    z = models.FloatField()
    occupancy = models.FloatField()
    bfactor = models.FloatField()
    anisou = models.FloatField(null=True)
    element = models.CharField(max_length=10)

    # class Meta:
    #     unique_together = (('residue', "serial"),)


class ResidueSet(models.Model):
    id = models.AutoField(primary_key=True)
    name = models.CharField(max_length=255)  # Pocket / CSA
    description = models.TextField(blank=True, default="")

    class Meta:
        unique_together = (("name"),)


class PDBResidueSet(models.Model):
    id = models.AutoField(primary_key=True)
    pdb = models.ForeignKey(PDB, related_name='residue_sets',
                            db_column="pdb_id", on_delete=models.DO_NOTHING)
    residue_set = models.ForeignKey(ResidueSet, related_name='residuesets',
                                    db_column="residueset_id", on_delete=models.DO_NOTHING)
    name = models.CharField(max_length=100,default="")

    class Meta:
        unique_together = (('pdb', "residue_set","name"),)

class ResidueSetResidue(models.Model):
    id = models.AutoField(primary_key=True)
    residue = models.ForeignKey(Residue, related_name='residue_sets',
                                db_column="residue_id", on_delete=models.DO_NOTHING)
    pdbresidue_set = models.ForeignKey(PDBResidueSet, related_name='residues',
                                       db_column="pdbresidueset_id", on_delete=models.DO_NOTHING)

    class Meta:
        unique_together = (('residue', "pdbresidue_set"),)

class AtomResidueSet(models.Model):
    id = models.AutoField(primary_key=True)
    atom = models.ForeignKey(Atom, db_column="atom_id", on_delete=models.DO_NOTHING)
    pdb_set = models.ForeignKey(ResidueSetResidue, related_name='atoms',
                                db_column="residuesetresidue_id", on_delete=models.DO_NOTHING)


    class Meta:
        unique_together = (('pdb_set', "atom"),)


class Property(models.Model):
    id = models.AutoField(primary_key=True)
    name = models.CharField(max_length=255)
    description = models.TextField(blank=True, default="")

    class Meta:
        unique_together = (('name'),)


class PropertyTag(models.Model):
    id = models.AutoField(primary_key=True)
    property = models.ForeignKey(Property, related_name='tags',
                                 db_column="residue_id", on_delete=models.DO_NOTHING)
    tag = models.CharField(max_length=50)
    description = models.TextField(blank=True, default="")

    class Meta:
        unique_together = (('property', 'tag'),)


class PDBProperty(models.Model):
    id = models.AutoField(primary_key=True)
    pdb = models.ForeignKey(PDB, related_name='properties',
                            db_column="pdb_id", on_delete=models.DO_NOTHING)
    property = models.ForeignKey(PDB, db_column="property_id", on_delete=models.DO_NOTHING)

    value = models.FloatField(null=True)
    tag = models.ForeignKey(PropertyTag, related_name='pdbs',
                            db_column="propertytag_id", null=True, on_delete=models.DO_NOTHING)

    class Meta:
        unique_together = (('pdb', "property", "tag"),)


class ResidueProperty(models.Model):
    id = models.AutoField(primary_key=True)
    residue = models.ForeignKey(Residue, related_name='properties',
                                db_column="residue_id", on_delete=models.DO_NOTHING)
    property = models.ForeignKey(PDB, db_column="property_id", on_delete=models.DO_NOTHING)
    value = models.FloatField(null=True)
    tag = models.ForeignKey(PropertyTag, related_name='residues',
                            db_column="propertytag_id", null=True, on_delete=models.DO_NOTHING)

    class Meta:
        unique_together = (('residue', "property", "tag"),)


class ResidueSetProperty(models.Model):
    id = models.AutoField(primary_key=True)
    pdbresidue_set = models.ForeignKey(PDBResidueSet, related_name='properties',
                                       db_column="pdbresidueset_id", on_delete=models.DO_NOTHING)
    property = models.ForeignKey(Property, related_name='residuesets', db_column="property_id", on_delete=models.DO_NOTHING)
    value = models.FloatField(null=True)
    tag = models.ForeignKey(PropertyTag, related_name='residue_sets',
                            db_column="propertytag_id", null=True, on_delete=models.DO_NOTHING)

    class Meta:
        unique_together = (('pdbresidue_set', "property", "tag"),)


class ChainProperty(models.Model):
    id = models.AutoField(primary_key=True)
    pdb = models.ForeignKey(PDB, related_name='chain_props',
                            db_column="pdb_id", on_delete=models.DO_NOTHING)
    chain = CustomBinaryCharField(max_length=10)
    property = models.ForeignKey(Property, related_name="chains", db_column="property_id", on_delete=models.DO_NOTHING)
    value = models.FloatField(null=True)
    tag = models.ForeignKey(PropertyTag, related_name='chains',
                            db_column="propertytag_id", null=True, on_delete=models.DO_NOTHING)

    class Meta:
        unique_together = (('pdb', "chain", "property", "tag"),)


class AtomProperty(models.Model):
    id = models.AutoField(primary_key=True)
    atom = models.ForeignKey(Atom, related_name='properties',
                             db_column="atom_id", on_delete=models.DO_NOTHING)
    property = models.ForeignKey(Property, db_column="property_id", on_delete=models.DO_NOTHING, related_name="atoms")
    value = models.FloatField(null=True)
    tag = models.ForeignKey(PropertyTag, related_name='atoms',
                            db_column="propertytag_id", null=True, on_delete=models.DO_NOTHING)

    class Meta:
        unique_together = (('atom', "property", "tag"),)
