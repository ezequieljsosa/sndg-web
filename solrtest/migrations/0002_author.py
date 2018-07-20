# Generated by Django 2.0 on 2018-07-04 17:18

from django.db import migrations, models
import django.db.models.deletion


class Migration(migrations.Migration):

    dependencies = [
        ('solrtest', '0001_initial'),
    ]

    operations = [
        migrations.CreateModel(
            name='Author',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('name', models.CharField(max_length=200)),
                ('num', models.IntegerField()),
                ('note', models.ForeignKey(null=True, on_delete=django.db.models.deletion.SET_NULL, related_name='notes', to='solrtest.Note')),
            ],
        ),
    ]