import datetime
import xmltodict
from Bio import Entrez
from bioresources.models import Sample, Publication, ExternalId, Resource, ResourceRelation, Organization, Person, \
    ReadsArchive, Assembly, Structure, Expression, Affiliation, ResourceProperty, ResourcePropertyValue
from biosql.models import Taxon, TaxonName, Term, Ontology
from django.db import transaction


def ncbi_link(dbfrom, id, linkname):
    pass


def scopus_publication(article):
    if Publication.objects.filter(scopus_id=article['dc:identifier']).exists():
        publication = Publication.objects.get(scopus_id=article['dc:identifier'])
    elif ("prism:doi" in article) and \
            Publication.objects.filter(doi=article['prism:doi']).exists():
        publication = Publication.objects.get(doi=article['prism:doi'])
    elif Publication.objects.filter(name=article["dc:title"][:350]).exists():
        publication = Publication.objects.get(name=article["dc:title"][:350])
    else:
        publication = Publication(
            type=Resource.PUBLICATION,
            name=article["dc:title"][:350],
            date_of_publication=datetime.datetime.strptime(article['prism:coverDate'], "%Y-%m-%d"),
        )
        publication.scopus_id = article['dc:identifier']
        publication.electronic_id = article["eid"]
        if "dc:description" in article:
            publication.description = article["dc:description"]
        if 'pubmed-id' in article:
            publication.pubmed_id = article['pubmed-id']
        if "prism:doi" in article:
            publication.doi = article["prism:doi"]
        if 'prism:issn' in article:
            publication.issn = article['prism:issn']
        publication.save()
    return publication


def scopus_affiliation(affiliation):
    afcountry = affiliation["affiliation-country"]

    org = Organization(name=affiliation["affilname"],
                       country=afcountry,
                       city=affiliation["affiliation-city"]
                       )

    # TODO country detection based on the city
    if not affiliation["affiliation-country"]:
        if any([x in str(org.city) for x in [
            "CABA", "Buenos Aires", "Rosario"
        ]]):
            org.country = "Argentina"

    if "afid" in affiliation:
        org.scopus_id = affiliation['afid']
        if Organization.objects.filter(scopus_id=org.scopus_id).exists():
            org = Organization.objects.get(scopus_id=org.scopus_id)
        else:
            org.save()

    return org


def scopus_author(author, publication, arg):
    person = Person(surname=author["surname"],
                    name=author["given-name"] if author["given-name"] else "",
                    scopus_id=author["authid"])
    if Person.objects.filter(scopus_id=person.scopus_id).exists():
        person = Person.objects.get(scopus_id=person.scopus_id)
    else:
        person.save()

    if ("afid" in author):

        aff = Affiliation(publication=publication, author=person)
        aff.save()
        for affdict in author["afid"]:
            aff.organizations.add(Organization.objects.get(scopus_id=affdict["$"]))
        aff.save()

        if [x for x in author["afid"] if x["$"] in arg]:
            person.save()
    return person


def scopus_extended_publication(article):
    publication = scopus_publication(article)
    arg = []  # TODO criterio pais?
    for affiliation in article["affiliation"]:
        org = scopus_affiliation(affiliation)
        if org.country == "Argentina" and org.scopus_id:
            arg.append(org.scopus_id)

    if not arg:
        # for a in article["affiliation"]:
        #     if not a["affiliation-country"]:
        #         if a["affiliation-city"]:
        #             pepe.append(a["affiliation-city"])
        return

    if "author" in article:
        for author in article["author"]:
            scopus_author(author, publication, arg)
    return publication


class NCBIBioSampleAdapter:

    def fetch(self, ncbi_id):
        biosampledata = Entrez.read(Entrez.esummary(db="biosample", id=ncbi_id))["DocumentSummarySet"][
            "DocumentSummary"][0]
        return xmltodict.parse(biosampledata["SampleData"])["BioSample"]

    def save(self, summaryData, ncbi_id):

        acc = summaryData["@accession"]
        desc = summaryData["Description"]["Title"]

        qs = Sample.objects.filter(external_ids__identifier=acc, external_ids__type="accession")
        if qs.exists():
            return qs.first()
        ncbi_org = Organization.objects.get(name="NCBI")

        ontology = Ontology.objects.get_or_create(name="NCBI sample", definition="Attributes of an NCBI Sample")[0]

        s = Sample(type=Resource.SAMPLE, name=acc, description=desc)

        if ("Organism" in summaryData["Description"]) and ("@taxonomy_id" in summaryData["Description"]["Organism"]):
            tax = summaryData["Description"]["Organism"]['@taxonomy_id']
            s.ncbi_tax = Taxon.objects.filter(ncbi_taxon_id=int(tax)).first()

        try:
            s.publication_date = datetime.datetime.strptime(summaryData["@publication_date"].split("T")[0], "%Y-%m-%d")
        except ValueError:
            pass
        try:
            s.update_date = datetime.datetime.strptime(summaryData["@last_update"].split("T")[0], "%Y-%m-%d")
        except ValueError:
            pass

        s.save()
        s.publishers.add(ncbi_org)

        for x in summaryData["Attributes"]["Attribute"]:
            if x["#text"].strip() and (x["#text"] != "."):
                term = Term.objects.get_or_create(
                    ontology=ontology, name=x["@attribute_name"], identifier=x["@attribute_name"][:255])[0]
                prop = ResourceProperty.objects.create(term=term, resource=s, organization=ncbi_org)
                ResourcePropertyValue.objects.create(property=prop, value=x["#text"][:200])

        ExternalId(resource=s, organization=ncbi_org,
                   identifier=acc, type="accession").save()
        ExternalId(resource=s, organization=ncbi_org,
                   identifier=ncbi_id, type="identifier").save()
        return s


class NCBISRAAdapter:

    def fetch(self, ncbi_id):
        return Entrez.read(Entrez.esummary(db="sra", id=ncbi_id))[0]

    def save(self, summaryData, ncbi_id):
        ncbi_org = Organization.objects.get(name="NCBI")
        expData = xmltodict.parse("<xml>" + Entrez.read(
            Entrez.esummary(db="sra", id=ncbi_id))[0]["ExpXml"] + "</xml>")["xml"]
        acc = expData["Experiment"]["@acc"]

        qs = ReadsArchive.objects.filter(external_ids__identifier=acc, external_ids__type="accession")
        if qs.exists():
            return qs.first()

        s = ReadsArchive(type=Resource.READS, name=acc, description=expData["Summary"]["Title"])
        try:
            s.release_date = datetime.datetime.strptime(summaryData["CreateDate"].split(" ")[0], "%Y/%m/%d")
        except ValueError:
            pass
        try:
            s.update_date = datetime.datetime.strptime(summaryData["UpdateDate"].split(" ")[0], "%Y/%m/%d")
        except ValueError:
            pass

        if ("Organism" in expData) and ("@taxid" in expData["Organism"]):
            s.ncbi_tax = Taxon.objects.filter(ncbi_taxon_id=int(expData["Organism"]["@taxid"])).first()

        s.save()
        s.publishers.add(ncbi_org)

        ExternalId(resource=s, organization=ncbi_org,
                   identifier=acc, type="accession").save()
        ExternalId(resource=s, organization=ncbi_org,
                   identifier=ncbi_id, type="identifier").save()
        return s


class NCBIAssemblyAdapter:

    def fetch(self, ncbi_id):
        return Entrez.read(Entrez.esummary(db="assembly", id=ncbi_id))["DocumentSummarySet"][
            "DocumentSummary"][0]

    def save(self, summaryData, ncbi_id):

        name = summaryData["AssemblyName"]
        acc = summaryData["AssemblyAccession"]
        qs = Assembly.objects.filter(external_ids__identifier=acc, external_ids__type="accession")
        if qs.exists():
            return qs.first()
        ncbi_org = Organization.objects.get(name="NCBI")

        tax = Taxon.objects.filter(ncbi_taxon_id=int(summaryData["Taxid"])).first()

        s = Assembly(type=Resource.ASSEMBLY, name=acc + "_" + name, description=summaryData["AssemblyDescription"],
                     ncbi_tax=tax,
                     ncbi_org=summaryData["SubmitterOrganization"],
                     level=summaryData["AssemblyStatus"],
                     assembly_type=summaryData["AssemblyType"],
                     species_name=summaryData["SpeciesName"])

        try:
            s.intraspecific_name = str(
                summaryData["Biosource"]["InfraspeciesList"][0]["Sub_type"]) + " " + \
                                   summaryData["Biosource"]["InfraspeciesList"][0]["Sub_value"]
        except IndexError:
            pass

        try:
            s.release_date = datetime.datetime.strptime(summaryData["SeqReleaseDate"].split(" ")[0],
                                                        "%Y/%m/%d")
        except ValueError:
            pass
        try:
            s.update_date = datetime.datetime.strptime(summaryData["LastUpdateDate"].split(" ")[0],
                                                       "%Y/%m/%d")
        except ValueError:
            pass

        s.save()
        s.publishers.add(ncbi_org)

        ExternalId(resource=s, organization=ncbi_org,
                   identifier=acc, type="accession").save()
        ExternalId(resource=s, organization=ncbi_org,
                   identifier=ncbi_id, type="identifier").save()
        return s


class NCBIStructureAdapter:

    def fetch(self, ncbi_id):
        return Entrez.read(Entrez.esummary(db="structure", id=ncbi_id))[0]

    def save(self, summaryData, ncbi_id):

        acc = summaryData["PdbAcc"]
        qs = Structure.objects.filter(external_ids__identifier=acc, external_ids__type="accession")
        if qs.exists():
            return qs.first()
        ncbi_org = Organization.objects.get(name="NCBI")

        s = Structure(type=Resource.STRUCTURE, name=acc, description=summaryData["PdbDescr"],
                      method=summaryData["ExpMethod"])
        s.publishers.add(ncbi_org)
        try:
            s.deposit_date = datetime.datetime.strptime(summaryData["PdbDepositDate"],
                                                        "%Y/%m/%d %H:%M")  # 2017/08/10 00:00
        except ValueError:
            pass
        if "OrganismList" in summaryData and summaryData["OrganismList"]:
            tax = TaxonName.objects.get(name=summaryData["OrganismList"][0])
            if tax:
                tax = tax.taxon
                s.ncbi_tax = tax

        s.save()
        s.publishers.add(ncbi_org)

        ExternalId(resource=s, organization=ncbi_org,
                   identifier=acc, type="accession").save()
        ExternalId(resource=s, organization=ncbi_org,
                   identifier=ncbi_id, type="identifier").save()
        return s


class NCBIGDSAdapter:

    def fetch(self, ncbi_id):
        return Entrez.read(Entrez.esummary(db="gds", id=ncbi_id))[0]

    def save(self, summaryData, ncbi_id):

        acc = summaryData["Accession"]
        qs = Expression.objects.filter(external_ids__identifier=acc, external_ids__type="accession")
        if qs.exists():
            return qs.first()
        ncbi_org = Organization.objects.get(name="NCBI")

        s = Expression(type=Resource.EXPRESSION, name=acc,
                       description=summaryData["title"] + "." + summaryData["summary"],
                       gdstype=summaryData["gdsType"])

        if "OrganismList" in summaryData:
            s.ncbi_org = "||".join(summaryData["OrganismList"])
        if "ExpMethod" in summaryData:
            s.method = str(summaryData["ExpMethod"])

        try:
            s.pdat = datetime.datetime.strptime(summaryData["PDAT"], "%Y/%m/%d")  # 2017/08/10 00:00
        except ValueError:
            pass

        tax = TaxonName.objects.get(name=summaryData["taxon"].split(";")[0])
        if tax:
            tax = tax.taxon
            s.ncbi_tax = tax

        s.save()
        s.publishers.add(ncbi_org)

        ExternalId(resource=s, organization=ncbi_org,
                   identifier=acc, type="accession").save()
        ExternalId(resource=s, organization=ncbi_org,
                   identifier=ncbi_id, type="identifier").save()
        return s
