from neomodel import StructuredNode, StringProperty, DateProperty, One, RelationshipTo, IntegerProperty

from django_neomodel import DjangoNode
from neomodel import db

"""MATCH (a:Artist),(b:Album)
WHERE a.Name = "Strapping Young Lad" AND b.Name = "Heavy as a Really Heavy Thing"
CREATE (a)-[r:RELEASED]->(b)
RETURN r"""


def connect_nodes(r1, r2, reltype: str = "USES"):
    from bioresources.models.Resource import Resource
    resource_class = {
        Resource.RESOURCE_TYPES.STRUCTURE: "Structure",
        Resource.RESOURCE_TYPES.ASSEMBLY: "Assembly",
        Resource.RESOURCE_TYPES.SAMPLE: "Sample",
        Resource.RESOURCE_TYPES.EXPRESSION: "Expression",
        Resource.RESOURCE_TYPES.PUBLICATION: "Publication",
        Resource.RESOURCE_TYPES.TOOL: "Tool",
        Resource.RESOURCE_TYPES.READS: "Reads",
        Resource.RESOURCE_TYPES.PERSON: "Person",

    }
    query = """ MATCH (a:%(src)s),(b:%(dst)s)
                WHERE a.rid = %(src_id)s AND b.rid = %(dst_id)s 
                CREATE (a)-[r:%(reltype)s]->(b)
                RETURN r""" % {
        "src": resource_class[r1.type],
        "dst": resource_class[r2.type],
        "src_id": r1.id, "dst_id": r2.id, "reltype": reltype
    }
    return db.cypher_query(query, {})


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
