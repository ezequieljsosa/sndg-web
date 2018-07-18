# Generated by Django 2.0 on 2018-07-06 19:40

from django.db import migrations, models
import django.db.models.deletion


class Migration(migrations.Migration):

    dependencies = [
        ('bioresources', '0002_auto_20180628_1652'),
    ]

    operations = [
        migrations.CreateModel(
            name='Assembly',
            fields=[
                ('resource_ptr', models.OneToOneField(auto_created=True, on_delete=django.db.models.deletion.CASCADE, parent_link=True, primary_key=True, serialize=False, to='bioresources.Resource')),
                ('level', models.CharField(max_length=50, null=True)),
                ('release_date', models.DateField(null=True)),
                ('assembly_type', models.CharField(choices=[('haploid', 'haploid'), ('diploid', 'diploid'), ('other', 'other')], max_length=50, null=True)),
            ],
            bases=('bioresources.resource',),
        ),
        migrations.CreateModel(
            name='Expression',
            fields=[
                ('resource_ptr', models.OneToOneField(auto_created=True, on_delete=django.db.models.deletion.CASCADE, parent_link=True, primary_key=True, serialize=False, to='bioresources.Resource')),
                ('pdat', models.DateField(null=True)),
                ('gdstype', models.CharField(max_length=50, null=True)),
                ('submitters', models.TextField(null=True)),
            ],
            bases=('bioresources.resource',),
        ),
        migrations.CreateModel(
            name='Structure',
            fields=[
                ('resource_ptr', models.OneToOneField(auto_created=True, on_delete=django.db.models.deletion.CASCADE, parent_link=True, primary_key=True, serialize=False, to='bioresources.Resource')),
                ('pdbClass', models.CharField(max_length=50, null=True)),
                ('deposit_date', models.DateField(null=True)),
                ('method', models.CharField(max_length=50, null=True)),
                ('org_list', models.TextField(null=True)),
            ],
            bases=('bioresources.resource',),
        ),
        migrations.AddField(
            model_name='bioproject',
            name='method',
            field=models.CharField(max_length=200, null=True),
        ),
        migrations.AddField(
            model_name='bioproject',
            name='submitters',
            field=models.TextField(null=True),
        ),
        migrations.AddField(
            model_name='bioproject',
            name='target',
            field=models.CharField(max_length=200, null=True),
        ),
        migrations.AlterField(
            model_name='affiliation',
            name='author',
            field=models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, related_name='authors', to='bioresources.Person'),
        ),
        migrations.AlterField(
            model_name='affiliation',
            name='publication',
            field=models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, related_name='affiliations', to='bioresources.Publication'),
        ),
        migrations.AlterField(
            model_name='bioproject',
            name='capture',
            field=models.CharField(choices=[('whole', 'whole'), ('exome', 'exome'), ('barcode', 'barcode'), ('TargetedLocusLoci', 'TargetedLocusLoci')], max_length=200, null=True),
        ),
        migrations.AlterField(
            model_name='bioproject',
            name='material',
            field=models.CharField(choices=[('genome', 'genome'), ('metagenome', 'metagenome'), ('chromosome', 'chromosome'), ('transcriptome', 'transcriptome'), ('reagent', 'reagent'), ('proteome', 'proteome')], max_length=20, null=True),
        ),
        migrations.AlterField(
            model_name='bioproject',
            name='sample_scope',
            field=models.CharField(choices=[('monoisolate', 'monoisolate'), ('multi-species', 'multi-species'), ('environment', 'environment'), ('synthetic', 'synthetic'), ('other', 'other')], max_length=20, null=True),
        ),
    ]
