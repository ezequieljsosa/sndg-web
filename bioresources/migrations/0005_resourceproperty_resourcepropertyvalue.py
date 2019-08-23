# Generated by Django 2.2.2 on 2019-07-22 16:05

from django.db import migrations, models
import django.db.models.deletion


class Migration(migrations.Migration):

    dependencies = [
        ('bioseq', '0001_initial'),
        ('bioresources', '0004_readsarchive_sample'),
    ]

    operations = [
        migrations.CreateModel(
            name='ResourceProperty',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('organization', models.ForeignKey(on_delete=django.db.models.deletion.DO_NOTHING, to='bioresources.Organization')),
                ('resource', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, related_name='properties', to='bioresources.Resource')),
                ('term', models.ForeignKey(on_delete=django.db.models.deletion.DO_NOTHING, to='bioseq.Term')),
            ],
        ),
        migrations.CreateModel(
            name='ResourcePropertyValue',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('value', models.CharField(max_length=200)),
                ('property', models.OneToOneField(on_delete=django.db.models.deletion.CASCADE, related_name='value', to='bioresources.ResourceProperty')),
            ],
        ),
    ]
