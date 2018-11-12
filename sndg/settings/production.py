# noinspection PyUnresolvedReferences
from .base import *


DEBUG = False
ALLOWED_HOSTS = ['localhost', ]
SECRET_KEY = env('DJANGO_SECRET_KEY')
