# Generated by Django 2.2.2 on 2019-08-02 18:01

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('bioresources', '0007_auto_20190726_1954'),
    ]

    operations = [
        migrations.RenameField(
            model_name='affiliation',
            old_name='publication',
            new_name='resource',
        ),
        migrations.AlterField(
            model_name='resource',
            name='type',
            field=models.PositiveIntegerField(choices=[(0, 'PUBLICATION'), (1, 'BIOPROJECT'), (2, 'SEQUENCE'), (3, 'ASSEMBLY'), (4, 'GENOME'), (5, 'READS'), (6, 'STRUCTURE'), (7, 'EXPRESSION'), (8, 'BARCODE'), (9, 'SAMPLE'), (10, 'TOOL'), (40, 'PROTEIN'), (30, 'ORGANIZATION'), (20, 'PERSON')]),
        ),
    ]
