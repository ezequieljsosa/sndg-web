# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.test import TestCase
from .models.Organization import Organization

# from bioresources.io.adapters import NCBIGDSAdapter


class AdaptersTestCase(TestCase):
    def setUp(self):
        pass

    def test_gds_adapter(self):
        Organization.objects.get(name="NCBI")
        # adapter = NCBIGDSAdapter()
        # summaryData = adapter.fetch(GSE107376)
        # adapter.save(summaryData)
        #
        # self.assertEqual(lion.speak(), 'The lion says "roar"')
        # self.assertEqual(cat.speak(), 'The cat says "meow"')
