# Generated by Django 2.2.2 on 2019-08-28 12:57

from django.db import migrations, models
import django.db.models.deletion


class Migration(migrations.Migration):

    dependencies = [
        ('bioresources', '0011_auto_20190827_1443'),
    ]

    operations = [
        migrations.CreateModel(
            name='UnprocessedResource',
            fields=[
                ('resource_ptr', models.OneToOneField(auto_created=True, on_delete=django.db.models.deletion.CASCADE, parent_link=True, primary_key=True, serialize=False, to='bioresources.Resource')),
                ('future_type', models.PositiveIntegerField(choices=[(0, 'PUBLICATION'), (1, 'BIOPROJECT'), (2, 'SEQUENCE'), (3, 'ASSEMBLY'), (4, 'GENOME'), (5, 'READS'), (6, 'STRUCTURE'), (7, 'EXPRESSION'), (8, 'BARCODE'), (9, 'SAMPLE'), (10, 'TOOL'), (50, 'UNPROCESSED'), (40, 'PROTEIN'), (30, 'ORGANIZATION'), (20, 'PERSON')])),
            ],
            bases=('bioresources.resource',),
        ),
        migrations.AlterField(
            model_name='resource',
            name='type',
            field=models.PositiveIntegerField(choices=[(0, 'PUBLICATION'), (1, 'BIOPROJECT'), (2, 'SEQUENCE'), (3, 'ASSEMBLY'), (4, 'GENOME'), (5, 'READS'), (6, 'STRUCTURE'), (7, 'EXPRESSION'), (8, 'BARCODE'), (9, 'SAMPLE'), (10, 'TOOL'), (50, 'UNPROCESSED'), (40, 'PROTEIN'), (30, 'ORGANIZATION'), (20, 'PERSON')]),
        ),
    ]