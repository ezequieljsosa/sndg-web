from django.contrib.auth.models import User
from django.db import models


class Note(models.Model):
    user = models.ForeignKey(User,null=True,on_delete=models.SET_NULL)
    pub_date = models.DateTimeField()
    title = models.CharField(max_length=200)
    body = models.TextField()

    def author_names(self):
        auths = []
        for author in self.notes.all():
             auths.append( author.name)
        return auths

    def __unicode__(self):
        return self.title

class Author(models.Model):
    note = models.ForeignKey(Note,null=True,on_delete=models.SET_NULL, related_name="notes")
    name =  models.CharField(max_length=200)
    num = models.IntegerField()
