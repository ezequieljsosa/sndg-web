# Generated by Django 2.2.2 on 2019-08-02 18:21

from django.db import migrations


class Migration(migrations.Migration):

    dependencies = [
        ('bioresources', '0010_auto_20190802_1820'),
    ]

    operations = [
        migrations.RenameField(
            model_name='affiliation',
            old_name='resource',
            new_name='publication',
        ),
    ]
