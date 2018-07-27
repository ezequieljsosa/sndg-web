# Generated by Django 2.0 on 2018-07-23 18:53

from django.db import migrations, models
import django.db.models.deletion


class Migration(migrations.Migration):

    dependencies = [
        ('bioresources', '0006_auto_20180713_1845'),
    ]

    operations = [
        migrations.AlterField(
            model_name='resourcerelation',
            name='source',
            field=models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, related_name='targets', to='bioresources.Resource'),
        ),
        migrations.AlterField(
            model_name='resourcerelation',
            name='target',
            field=models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, related_name='sources', to='bioresources.Resource'),
        ),
        # migrations.AlterUniqueTogether(
        #     name='resourcerelation',
        #     unique_together={('source', 'target', 'deprecated')},
        # ),
    ]
