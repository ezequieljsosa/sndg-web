# Generated by Django 2.2.2 on 2019-08-02 18:48

from django.db import migrations


class Migration(migrations.Migration):

    dependencies = [
        ('bioresources', '0011_auto_20190802_1821'),
    ]

    operations = [
        migrations.RemoveField(
            model_name='affiliation',
            name='organizations',
        ),
    ]
