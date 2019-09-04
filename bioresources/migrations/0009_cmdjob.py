# Generated by Django 2.2.2 on 2019-08-27 12:54

from django.db import migrations, models
import django.db.models.deletion


class Migration(migrations.Migration):

    dependencies = [
        ('bioresources', '0008_auto_20190827_1253'),
    ]

    operations = [
        migrations.CreateModel(
            name='CmdJob',
            fields=[
                ('job_ptr', models.OneToOneField(auto_created=True, on_delete=django.db.models.deletion.CASCADE, parent_link=True, primary_key=True, serialize=False, to='bioresources.Job')),
                ('command', models.TextField()),
            ],
            bases=('bioresources.job',),
        ),
    ]