# Generated by Django 2.0 on 2018-07-19 15:35

import datetime
from django.db import migrations, models
import django.db.models.deletion


class Migration(migrations.Migration):

    initial = True

    dependencies = [
    ]

    operations = [
        migrations.CreateModel(
            name='Atom',
            fields=[
                ('id', models.AutoField(primary_key=True, serialize=False)),
                ('serial', models.IntegerField()),
                ('name', models.CharField(max_length=3)),
                ('x', models.FloatField()),
                ('y', models.FloatField()),
                ('z', models.FloatField()),
                ('occupancy', models.FloatField()),
                ('bfactor', models.FloatField()),
                ('anisou', models.FloatField(null=True)),
                ('element', models.CharField(max_length=10)),
            ],
        ),
        migrations.CreateModel(
            name='AtomProperty',
            fields=[
                ('id', models.AutoField(primary_key=True, serialize=False)),
                ('value', models.FloatField(null=True)),
                ('atom', models.ForeignKey(db_column='atom_id', on_delete=django.db.models.deletion.DO_NOTHING, related_name='properties', to='pdbdb.Atom')),
            ],
        ),
        migrations.CreateModel(
            name='ChainProperty',
            fields=[
                ('id', models.AutoField(primary_key=True, serialize=False)),
                ('chain', models.CharField(max_length=10)),
                ('value', models.FloatField(null=True)),
            ],
        ),
        migrations.CreateModel(
            name='PDB',
            fields=[
                ('id', models.AutoField(primary_key=True, serialize=False)),
                ('code', models.CharField(max_length=100)),
                ('resolution', models.FloatField(default=20)),
                ('experiment', models.CharField(max_length=100, null=True)),
                ('tax', models.CharField(max_length=255, null=True)),
                ('deprecated', models.BooleanField(default=False)),
                ('date', models.DateTimeField(default=datetime.datetime.now)),
            ],
        ),
        migrations.CreateModel(
            name='PDBProperty',
            fields=[
                ('id', models.AutoField(primary_key=True, serialize=False)),
                ('value', models.FloatField(null=True)),
                ('pdb', models.ForeignKey(db_column='pdb_id', on_delete=django.db.models.deletion.DO_NOTHING, related_name='properties', to='pdbdb.PDB')),
                ('property', models.ForeignKey(db_column='property_id', on_delete=django.db.models.deletion.DO_NOTHING, to='pdbdb.PDB')),
            ],
        ),
        migrations.CreateModel(
            name='PDBResidueSet',
            fields=[
                ('id', models.AutoField(primary_key=True, serialize=False)),
                ('pdb', models.ForeignKey(db_column='pdb_id', on_delete=django.db.models.deletion.DO_NOTHING, related_name='residue_sets', to='pdbdb.PDB')),
            ],
        ),
        migrations.CreateModel(
            name='Property',
            fields=[
                ('id', models.AutoField(primary_key=True, serialize=False)),
                ('name', models.CharField(max_length=255)),
                ('description', models.TextField(blank=True, default='')),
            ],
        ),
        migrations.CreateModel(
            name='PropertyTag',
            fields=[
                ('id', models.AutoField(primary_key=True, serialize=False)),
                ('tag', models.CharField(max_length=50)),
                ('description', models.TextField(blank=True, default='')),
                ('property', models.ForeignKey(db_column='residue_id', on_delete=django.db.models.deletion.DO_NOTHING, related_name='tags', to='pdbdb.Property')),
            ],
        ),
        migrations.CreateModel(
            name='Residue',
            fields=[
                ('id', models.AutoField(primary_key=True, serialize=False)),
                ('chain', models.CharField(max_length=20)),
                ('resname', models.CharField(max_length=4)),
                ('resid', models.IntegerField()),
                ('icode', models.CharField(default='', max_length=2)),
                ('type', models.CharField(choices=[('H', 'HETATOM'), ('A', 'ATOM')], max_length=1)),
                ('disordered', models.BooleanField(default=False)),
                ('modelable', models.BooleanField(default=True)),
                ('seq_order', models.IntegerField(null=True)),
                ('pdb', models.ForeignKey(db_column='pdb_id', on_delete=django.db.models.deletion.DO_NOTHING, related_name='residues', to='pdbdb.PDB')),
            ],
        ),
        migrations.CreateModel(
            name='ResidueProperty',
            fields=[
                ('id', models.AutoField(primary_key=True, serialize=False)),
                ('value', models.FloatField(null=True)),
                ('property', models.ForeignKey(db_column='property_id', on_delete=django.db.models.deletion.DO_NOTHING, to='pdbdb.PDB')),
                ('residue', models.ForeignKey(db_column='residue_id', on_delete=django.db.models.deletion.DO_NOTHING, related_name='properties', to='pdbdb.Residue')),
                ('tag', models.ForeignKey(db_column='propertytag_id', null=True, on_delete=django.db.models.deletion.DO_NOTHING, related_name='residues', to='pdbdb.PropertyTag')),
            ],
        ),
        migrations.CreateModel(
            name='ResidueSet',
            fields=[
                ('id', models.AutoField(primary_key=True, serialize=False)),
                ('name', models.CharField(max_length=255)),
                ('description', models.TextField(blank=True, default='')),
            ],
        ),
        migrations.CreateModel(
            name='ResidueSetProperty',
            fields=[
                ('id', models.AutoField(primary_key=True, serialize=False)),
                ('value', models.FloatField(null=True)),
                ('pdbresidue_set', models.ForeignKey(db_column='pdbresidueset_id', on_delete=django.db.models.deletion.DO_NOTHING, related_name='properties', to='pdbdb.PDBResidueSet')),
                ('property', models.ForeignKey(db_column='property_id', on_delete=django.db.models.deletion.DO_NOTHING, related_name='residuesets', to='pdbdb.PDB')),
                ('tag', models.ForeignKey(db_column='propertytag_id', null=True, on_delete=django.db.models.deletion.DO_NOTHING, related_name='residue_sets', to='pdbdb.PropertyTag')),
            ],
        ),
        migrations.CreateModel(
            name='ResidueSetResidue',
            fields=[
                ('id', models.AutoField(primary_key=True, serialize=False)),
                ('pdbresidue_set', models.ForeignKey(db_column='pdbresidueset_id', on_delete=django.db.models.deletion.DO_NOTHING, related_name='residues', to='pdbdb.PDBResidueSet')),
                ('residue', models.ForeignKey(db_column='residue_id', on_delete=django.db.models.deletion.DO_NOTHING, related_name='residue_sets', to='pdbdb.Residue')),
            ],
        ),
        migrations.AlterUniqueTogether(
            name='residueset',
            unique_together={('name',)},
        ),
        migrations.AlterUniqueTogether(
            name='property',
            unique_together={('name',)},
        ),
        migrations.AddField(
            model_name='pdbresidueset',
            name='residue_set',
            field=models.ForeignKey(db_column='residueset_id', on_delete=django.db.models.deletion.DO_NOTHING, related_name='pdbs', to='pdbdb.ResidueSet'),
        ),
        migrations.AddField(
            model_name='pdbproperty',
            name='tag',
            field=models.ForeignKey(db_column='propertytag_id', null=True, on_delete=django.db.models.deletion.DO_NOTHING, related_name='pdbs', to='pdbdb.PropertyTag'),
        ),
        migrations.AlterUniqueTogether(
            name='pdb',
            unique_together={('code',)},
        ),
        migrations.AddField(
            model_name='chainproperty',
            name='pdb',
            field=models.ForeignKey(db_column='pdb_id', on_delete=django.db.models.deletion.DO_NOTHING, related_name='chain_props', to='pdbdb.PDB'),
        ),
        migrations.AddField(
            model_name='chainproperty',
            name='property',
            field=models.ForeignKey(db_column='property_id', on_delete=django.db.models.deletion.DO_NOTHING, related_name='chains', to='pdbdb.PDB'),
        ),
        migrations.AddField(
            model_name='chainproperty',
            name='tag',
            field=models.ForeignKey(db_column='propertytag_id', null=True, on_delete=django.db.models.deletion.DO_NOTHING, related_name='chains', to='pdbdb.PropertyTag'),
        ),
        migrations.AddField(
            model_name='atomproperty',
            name='property',
            field=models.ForeignKey(db_column='property_id', on_delete=django.db.models.deletion.DO_NOTHING, related_name='atoms', to='pdbdb.PDB'),
        ),
        migrations.AddField(
            model_name='atomproperty',
            name='tag',
            field=models.ForeignKey(db_column='propertytag_id', null=True, on_delete=django.db.models.deletion.DO_NOTHING, related_name='atoms', to='pdbdb.PropertyTag'),
        ),
        migrations.AddField(
            model_name='atom',
            name='residue',
            field=models.ForeignKey(db_column='residue_id', on_delete=django.db.models.deletion.DO_NOTHING, related_name='atoms', to='pdbdb.Residue'),
        ),
        migrations.AlterUniqueTogether(
            name='residuesetresidue',
            unique_together={('residue', 'pdbresidue_set')},
        ),
        migrations.AlterUniqueTogether(
            name='residuesetproperty',
            unique_together={('pdbresidue_set', 'property', 'tag')},
        ),
        migrations.AlterUniqueTogether(
            name='residueproperty',
            unique_together={('residue', 'property', 'tag')},
        ),
        migrations.AlterUniqueTogether(
            name='residue',
            unique_together={('pdb', 'chain', 'resid', 'icode', 'type')},
        ),
        migrations.AlterUniqueTogether(
            name='propertytag',
            unique_together={('property', 'tag')},
        ),
        migrations.AlterUniqueTogether(
            name='pdbresidueset',
            unique_together={('pdb', 'residue_set')},
        ),
        migrations.AlterUniqueTogether(
            name='pdbproperty',
            unique_together={('pdb', 'property', 'tag')},
        ),
        migrations.AlterUniqueTogether(
            name='chainproperty',
            unique_together={('pdb', 'chain', 'property', 'tag')},
        ),
        migrations.AlterUniqueTogether(
            name='atomproperty',
            unique_together={('atom', 'property', 'tag')},
        ),
        migrations.AlterUniqueTogether(
            name='atom',
            unique_together={('residue', 'serial')},
        ),
    ]
