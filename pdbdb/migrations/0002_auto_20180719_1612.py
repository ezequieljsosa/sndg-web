# Generated by Django 2.0 on 2018-07-19 16:12

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('pdbdb', '0001_initial'),
    ]

    operations = [
        migrations.AlterField(
            model_name='atom',
            name='name',
            field=models.CharField(max_length=50),
        ),
        migrations.AlterField(
            model_name='residue',
            name='type',
            field=models.CharField(max_length=10),
        ),
    ]
