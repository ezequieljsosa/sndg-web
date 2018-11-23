"""
Django settings for sndg project.

Generated by 'django-admin startproject' using Django 1.11.13.

For more information on this file, see
https://docs.djangoproject.com/en/1.11/topics/settings/

For the full list of settings and their values, see
https://docs.djangoproject.com/en/1.11/ref/settings/
"""
# noinspection PyUnresolvedReferences
import environ
import os

env = environ.Env()
environ.Env.read_env()

# CREATE DATABASE mydatabase CHARACTER SET utf8mb4 COLLATE utf8mb4_unicode_ci;
# DATABASE_URL=mysql://user:%23password@127.0.0.1:3306/dbname
DATABASES = {
    'default': env.db(),
}

# Build paths inside the project like this: os.path.join(BASE_DIR, ...)
BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# Quick-start development settings - unsuitable for production
# See https://docs.djangoproject.com/en/1.11/howto/deployment/checklist/

ALLOWED_HOSTS = []

# Application definition

INSTALLED_APPS = [
    'django.contrib.admin',
    'django.contrib.auth',
    'django.contrib.contenttypes',
    'django.contrib.sessions',
    'django.contrib.messages',
    'django.contrib.staticfiles',
    'django_extensions',

    'biosql',
    'pdbdb',
    # 'vardb',
    'bioresources',  # importante que este antes de allauth

    'django_select2',
    'easy_select2',

    'haystack',
    # 'solrtest',

    'django_filters',
    'bootstrap3',
    'django_tables2',
    'crudbuilder',

    'django.contrib.sites',

    'allauth',
    'allauth.account',
    'allauth.socialaccount',
    'allauth.socialaccount.providers.google',
    'captcha',
    'crispy_forms',

    "resumable",

]

MIDDLEWARE = [
    'django.middleware.security.SecurityMiddleware',
    'django.contrib.sessions.middleware.SessionMiddleware',
    'django.middleware.common.CommonMiddleware',
    'django.middleware.csrf.CsrfViewMiddleware',
    'django.contrib.auth.middleware.AuthenticationMiddleware',
    'django.contrib.messages.middleware.MessageMiddleware',
    'django.middleware.clickjacking.XFrameOptionsMiddleware',

]

ROOT_URLCONF = 'sndg.urls'

TEMPLATES = [
    {
        'BACKEND': 'django.template.backends.django.DjangoTemplates',
        'DIRS': [os.path.join(BASE_DIR, 'templates')],
        'APP_DIRS': True,
        'OPTIONS': {
            'context_processors': [
                'django.template.context_processors.debug',
                'django.template.context_processors.request',
                'django.contrib.auth.context_processors.auth',
                'django.contrib.messages.context_processors.messages',

            ],
        },
    },
]

WSGI_APPLICATION = 'sndg.wsgi.application'

# Password validation
# https://docs.djangoproject.com/en/1.11/ref/settings/#auth-password-validators

AUTH_PASSWORD_VALIDATORS = [
    {
        'NAME': 'django.contrib.auth.password_validation.UserAttributeSimilarityValidator',
    },
    {
        'NAME': 'django.contrib.auth.password_validation.MinimumLengthValidator',
    },
    {
        'NAME': 'django.contrib.auth.password_validation.CommonPasswordValidator',
    },
    {
        'NAME': 'django.contrib.auth.password_validation.NumericPasswordValidator',
    },
]

# Internationalization
# https://docs.djangoproject.com/en/1.11/topics/i18n/

LANGUAGE_CODE = 'en-us'

TIME_ZONE = 'UTC'

USE_I18N = True

USE_L10N = True

USE_TZ = True

# Static files (CSS, JavaScript, Images)
# https://docs.djangoproject.com/en/1.11/howto/static-files/
STATIC_ROOT = '/tmp/'
STATIC_URL = '/static/'
SITE_ROOT = os.path.dirname(os.path.realpath(__file__))

STATICFILES_DIRS = (
    os.path.abspath(os.path.join(SITE_ROOT, "../static/")),
    ("jbrowse", "/data/xomeq/JBrowse-1.14.2/"),
)

LOGIN_REDIRECT_URL = '/'

# AUTH_USER_MODEL = "bioresources.User"

AUTHENTICATION_BACKENDS = [
    # Needed to login by username in Django admin, regardless of `allauth`
    'django.contrib.auth.backends.ModelBackend',

    # `allauth` specific authentication methods, such as login by e-mail
    # 'allauth.account.auth_backends.AuthenticationBackend',
]
ACCOUNT_EMAIL_REQUIRED = True
ACCOUNT_USERNAME_REQUIRED = False
ACCOUNT_SIGNUP_FORM_CLASS = 'bioresources.forms.AllauthSignupForm'
SITE_ID = 1

EMAIL_BACKEND = 'django.core.mail.backends.smtp.EmailBackend'
EMAIL_USE_TLS = True
EMAIL_HOST = 'smtp.gmail.com'
EMAIL_PORT = 587
EMAIL_HOST_USER = DEFAULT_FROM_EMAIL = 'xxx@gmail.com'
EMAIL_HOST_PASSWORD = 'xxx'

HAYSTACK_CONNECTIONS = {
    # SNDG_SOLR=solr://127.0.0.1:8984/solr/Notes
    'default': {
        'ENGINE': 'haystack.backends.solr_backend.SolrEngine',
        'URL': env('SNDG_SOLR'),  # 'http://127.0.0.1:8984/solr/Notes'
        'INCLUDE_SPELLING': True,
        # ...or for multicore...
        # 'URL': 'http://127.0.0.1:8983/solr/mysite',
        'EXCLUDED_INDEXES': ['bioresources.search_indexes.PublicationIndexOAI'],
    },
    'oai': {
        'ENGINE': 'haystack.backends.solr_backend.SolrEngine',
        'URL': 'http://127.0.0.1:8984/solr/oai',
        'INCLUDE_SPELLING': False,
        'EXCLUDED_INDEXES': ['bioresources.search_indexes.PublicationIndex',
                             'bioresources.search_indexes.StructureIndex',
                             'bioresources.search_indexes.AssemblyIndex',
                             'bioresources.search_indexes.ExpressionIndex',
                             'bioresources.search_indexes.ToolIndex',
                             'bioresources.search_indexes.BioProjectIndex',
                             'bioresources.search_indexes.PersonIndex',
                             'bioresources.search_indexes.OrganizationIndex',
                             'bioresources.search_indexes.BarcodeIndex',
                             ],

    },
}


HAYSTACK_ID_FIELD = env("HAYSTACK_ID_FIELD",default='id')
HAYSTACK_DJANGO_CT_FIELD = env("HAYSTACK_DJANGO_CT_FIELD",default='django_ct')
HAYSTACK_DJANGO_ID_FIELD = env("HAYSTACK_DJANGO_ID_FIELD",default='django_id')
HAYSTACK_DOCUMENT_FIELD = env("HAYSTACK_DOCUMENT_FIELD",default='text')

CRISPY_TEMPLATE_PACK = 'bootstrap4'
FILE_UPLOAD_TEMP_DIR = "/tmp/pepe"

OAIPMH_DOMAIN = "sndg.qb.fcen.uba.ar"

# LOCALE_PATHS = os.path.abspath(os.path.join(SITE_ROOT, "../locale")),
