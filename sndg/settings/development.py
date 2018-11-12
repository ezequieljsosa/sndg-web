# noinspection PyUnresolvedReferences
from .base import *
DEBUG = True

INSTALLED_APPS += [
    'debug_toolbar',
]

MIDDLEWARE += ['debug_toolbar.middleware.DebugToolbarMiddleware', ]

EMAIL_BACKEND = 'django.core.mail.backends.console.EmailBackend'

DEBUG_TOOLBAR_CONFIG = {
    'JQUERY_URL': '',
}

SECRET_KEY = env('DJANGO_SECRET_KEY', default='w9q=dji%@_zrg$z-t&+hd2y^dkw2j56x!!ayw310nqd^v22ttb')