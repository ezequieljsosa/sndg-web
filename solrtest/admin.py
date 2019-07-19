from django.contrib import admin

from .models import Note, Author

admin.site.register(Note)
admin.site.register(Author)

