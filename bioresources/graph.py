from neomodel import StructuredNode, StringProperty, DateProperty, One, RelationshipTo, IntegerProperty

from django_neomodel import DjangoNode
from neomodel import db

"""MATCH (a:Artist),(b:Album)
WHERE a.Name = "Strapping Young Lad" AND b.Name = "Heavy as a Really Heavy Thing"
CREATE (a)-[r:RELEASED]->(b)
RETURN r"""


def get_resource_class():
    from bioresources.models.Resource import Resource
    return {
        Resource.RESOURCE_TYPES.STRUCTURE: "Structure",
        Resource.RESOURCE_TYPES.ASSEMBLY: "Assembly",
        Resource.RESOURCE_TYPES.SAMPLE: "Sample",
        Resource.RESOURCE_TYPES.EXPRESSION: "Expression",
        Resource.RESOURCE_TYPES.PUBLICATION: "Publication",
        Resource.RESOURCE_TYPES.TOOL: "Tool",
        Resource.RESOURCE_TYPES.READS: "Reads",
        Resource.RESOURCE_TYPES.PERSON: "Person",

    }


def delete_edge(r1, r2, reltype: str = "USES"):
    resource_class = get_resource_class()
    query = """ MATCH (a:%(src)s)-[e:%(reltype)]-(b:%(dst)s)
                WHERE a.rid = %(src_id)s AND b.rid = %(dst_id)s                 
                DELETE e""" % {
        "src": resource_class[r1.type],
        "dst": resource_class[r2.type],
        "src_id": r1.id, "dst_id": r2.id, "reltype": reltype
    }
    return db.cypher_query(query, {})


def connect_nodes(r1, r2, reltype: str = "USES"):
    resource_class = get_resource_class()

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

    @classmethod
    def from_resource(cls, organization):
        return cls(rid=organization.id, name=organization.name)


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

    @classmethod
    def from_resource(cls, publication):

        gPublication = Publication(rid=publication.id, title=publication.name,
                                   publication=publication.date_of_publication)
        gPublication.save()
        for aff in publication.affiliations.all():
            per = Person.nodes.get(rid=aff.author_id)
            gPublication.authors.connect(per)
            for org in aff.organizations.all():
                gOrg = Organization.nodes.get(rid=org.id)
                gPublication.organizations.connect(gOrg)


class Expression(Resource):
    pdat = StringProperty()
    gdstype = StringProperty()

    @classmethod
    def from_resource(cls, expresion):
        r = Expression(rid=expresion.id, title=expresion.name, pdat=expresion.pdat, gdstype=expresion.gdstype)
        r.save()


class Assembly(Resource):
    intraspecific_name = StringProperty()
    level = StringProperty()
    reads = RelationshipTo('Reads', 'READS')
    samples = RelationshipTo('Sample', 'SAMPLE')

    @classmethod
    def from_resource(cls, assembly):
        r = Assembly(rid=assembly.id, title=assembly.name, intraspecific_name=assembly.intraspecific_name)
        r.save()
        s = Species.nodes.get(name=assembly.species_name)
        r.species.connect(s)


class Barcodes(Resource):
    country = RelationshipTo('Country', 'COUNTRY')
    subdivision = StringProperty()
    marker = StringProperty()

    # @classmethod
    # def from_resource(cls, sample):


class Structure(Resource):
    deposit_date = DateProperty()
    method = StringProperty()

    @classmethod
    def from_resource(cls, sample):
        g = cls(rid=sample.id, title=sample.name, subdivision=sample.subdivision,
                collection_date=sample.collection_date)
        g.save()
        if sample.country:
            c = Country.nodes.get(name=sample.country)
            r.country.connect(c)


class Sample(Resource):
    country = RelationshipTo('Country', 'COUNTRY')
    subdivision = StringProperty()
    collection_date = DateProperty()

    @classmethod
    def from_resource(cls, sample):
        g = cls(rid=sample.id, title=sample.name, subdivision=sample.subdivision,
                collection_date=sample.collection_date)
        g.save()
        if sample.country:
            c = Country.nodes.get(name=sample.country)
            r.country.connect(c)


class Reads(Resource):
    sample = RelationshipTo('Sample', 'SOURCE', cardinality=One)

    @classmethod
    def from_resource(cls, reads):
        g = cls(rid=reads.id, title=reads.name)
        g.save()


class Tool(Resource):
    tool_type = StringProperty()

    @classmethod
    def from_resource(cls, tool):
        g = cls(rid=tool.id, title=tool.name, tool_type=tool.tool_type)
        g.save()


from django.db.models.signals import post_save, post_delete, post_init
from django.dispatch import receiver
from bioresources.models.Resource import Collaboration

from bioresources.models.Person import Person as rPerson


@receiver(post_save, sender=Collaboration)
def my_handler(sender, **kwargs):
    c = Collaboration.objects.prefetch_related("person", "resource").get(id=kwargs["instance"].id)
    c.person.type = rPerson.TYPE
    connect_nodes(c.person, c.resource, reltype=Collaboration.rev_types[c.type])


@receiver(post_delete, sender=Collaboration)
def my_handler2(sender, **kwargs):
    c = Collaboration.objects.prefetch_related("person", "resource").get(id=kwargs["instance"].id)
    c.person.type = Person.TYPE
    delete_edge(c.person, c.resource, reltype=Collaboration.rev_types[c.type])


from bioresources.models.Resource import Resource
from bioresources.models.ReadsArchive import ReadsArchive as rReadsArchive
from bioresources.models.Sample import Sample as rSample
from bioresources.models.Assembly import Assembly as rAssembly
from bioresources.models.Expression import Expression as rExpression
from bioresources.models.Barcode import Barcode as rBarcode
from bioresources.models.Structure import Structure as rStructure
from bioresources.models.Tool import Tool as rTool

gclass_dict = {
    Resource.RESOURCE_TYPES.STRUCTURE: Structure,
    Resource.RESOURCE_TYPES.ASSEMBLY: Assembly,
    Resource.RESOURCE_TYPES.SAMPLE: Sample,
    Resource.RESOURCE_TYPES.EXPRESSION: Expression,
    Resource.RESOURCE_TYPES.PUBLICATION: Publication,
    Resource.RESOURCE_TYPES.TOOL: Tool,
    Resource.RESOURCE_TYPES.READS: Reads
}


# TODO habilitar barcodes
@receiver(post_save, sender=[rTool, rReadsArchive, rSample, rAssembly, rExpression, rStructure])  # rBarcode
def my_handler3(sender, **kwargs):
    r = kwargs["instance"]

    if len(gclass.nodes.filter(rid=r.id)) == 0:
        gclass = gclass_dict[sender.TYPE]
        gclass.from_resource(r)
