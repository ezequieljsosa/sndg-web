from django.urls import path

from . import views

app_name = 'vardb'
urlpatterns = [
    path('var/<int:pk>', views.variant, name='var_view'),



]