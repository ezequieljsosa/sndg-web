"""sndg URL Configuration

The `urlpatterns` list routes URLs to views. For more information please see:
    https://docs.djangoproject.com/en/1.11/topics/http/urls/
Examples:
Function views
    1. Add an import:  from my_app import views
    2. Add a URL to urlpatterns:  url(r'^$', views.home, name='home')
Class-based views
    1. Add an import:  from other_app.views import Home
    2. Add a URL to urlpatterns:  url(r'^$', Home.as_view(), name='home')
Including another URLconf
    1. Import the include() function: from django.conf.urls import url, include
    2. Add a URL to urlpatterns:  url(r'^blog/', include('blog.urls'))
"""
from django.conf.urls import url
from django.contrib import admin
from django.urls import path,include
from django.contrib.auth import views as auth_views


# from bioresources import views as bioviews

urlpatterns = [
    url(r'^admin/', admin.site.urls),
    path('', include('bioresources.urls')),
    path('bioseq/', include('biosql.urls')),
    path('pdbdb/', include('pdbdb.urls')),
    path('strainq/', include('vardb.urls')),
    url(r'^crud/',  include('crudbuilder.urls')),


    url(r'^accounts/', include('allauth.urls')),
    url(r'^captcha/', include('captcha.urls')),

    url(r'^select2/', include('django_select2.urls')),


    # url(r'^login/$', auth_views.login, name='login'),
    # url(r'^logout/$', auth_views.logout,{'next_page': '/'}, name='logout'),
    # url(r'^signup/$', bioviews.signup, name='signup'),
    # url(r'^search/', include('haystack.urls')),
    # url(r'^solr/', include('solrtest.urls')),



]
