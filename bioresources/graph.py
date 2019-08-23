from neomodel import StructuredNode, StringProperty, DateProperty, One, RelationshipTo,IntegerProperty

from django_neomodel import DjangoNode


class Tax(StructuredNode):
    rid = IntegerProperty(unique_index=True)
    name = StringProperty()
    ncbi_id = StringProperty(unique_index=True)


class Country(StructuredNode):
    name = StringProperty(unique_index=True)


class Organization(StructuredNode):
    rid = IntegerProperty(unique_index=True)
    name = StringProperty(unique_index=True)
    level = DateProperty()
    location = RelationshipTo('Country', 'LOCATION')


class Journal(StructuredNode):
    name = StringProperty(unique_index=True)


class Person(StructuredNode):
    rid = IntegerProperty(unique_index=True)
    name = StringProperty()


class Species(StructuredNode):
    name = StringProperty(unique_index=True)


class Resource(StructuredNode):
    rid = IntegerProperty(unique_index=True)
    title = StringProperty(unique_index=True)
    # description = StringProperty()
    species = RelationshipTo('Species', 'SPECIES')
    tax = RelationshipTo('Tax', 'TAX')


class Publication(Resource):
    title = StringProperty(unique_index=True)
    # published = DateProperty()
    authors = RelationshipTo('Person', 'AUTHOR')
    organizations = RelationshipTo('Organization', 'AFF')
    journal = RelationshipTo('Journal', 'PUBLISHER', cardinality=One)
    resources = RelationshipTo('Resource', 'USED')


class Expression(Resource):
    pdat = StringProperty()
    gdstype = StringProperty()


class Assembly(Resource):
    intraspecific_name = StringProperty()
    level = StringProperty()
    reads = RelationshipTo('Reads', 'READS')
    samples = RelationshipTo('Sample', 'SAMPLE')


class Barcodes(Resource):
    country = RelationshipTo('Country', 'COUNTRY')
    subdivision = StringProperty()
    marker = StringProperty()


class Structure(Resource):
    deposit_date = DateProperty()
    method = StringProperty()


class Sample(Resource):
    country = RelationshipTo('Country', 'COUNTRY')
    subdivision = StringProperty()
    collection_date = DateProperty()

class Reads(Resource):
    sample = RelationshipTo('Sample', 'SOURCE', cardinality=One)



class Tool(Resource):
    tool_type = StringProperty()

