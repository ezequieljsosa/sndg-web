from Bio import Entrez

from .adapters import NCBIAssemblyAdapter, NCBIPubmedAdapter, NCBIPmcAdapter, NCBIBioSampleAdapter, NCBIBioProject, \
    NCBIStructureAdapter, NCBIGDSAdapter, NCBISRAAdapter, NCBINuccoreAdapter, NCBIGeneAdapter, NCBIProteinAdapter

from bioresources.models.Resource import Resource


class NCBISearch():
    allowed_dbs = ["gds", "pmc", "pubmed", "bioproject", "biosample", "sra", "taxonomy", "structure", "genome",
                   "protein", "gene", "nuccore", "assembly"]
    database_map = {"assembly": NCBIAssemblyAdapter(), "pmc": NCBIPmcAdapter(), "pubmed": NCBIPubmedAdapter(),
                    "gene": NCBIGeneAdapter(), "protein": NCBIProteinAdapter(), "structure": NCBIStructureAdapter(),
                    "gds": NCBIGDSAdapter(), "bioproject": NCBIBioProject(),
                    "nuccore": NCBINuccoreAdapter(), "sra": NCBISRAAdapter(), "biosample": NCBIBioSampleAdapter()

                    }
    db_type = {
        "assembly": Resource.RESOURCE_TYPES.ASSEMBLY, "pmc": Resource.RESOURCE_TYPES.PUBLICATION,
        "pubmed": Resource.RESOURCE_TYPES.PUBLICATION,
        "gene": 40, "protein": 40, "structure": Resource.RESOURCE_TYPES.STRUCTURE,
        "gds": Resource.RESOURCE_TYPES.EXPRESSION, "bioproject": Resource.RESOURCE_TYPES.BIOPROJECT,
        "nuccore": 40, "sra": Resource.RESOURCE_TYPES.READS, "biosample": Resource.RESOURCE_TYPES.SAMPLE
    }

    def search_all(self, search):
        result = {}
        data = Entrez.read(Entrez.egquery(term=search))["eGQueryResult"]
        for dbResult in data:
            if dbResult["DbName"] in NCBISearch.allowed_dbs:
                try:
                    result[dbResult["DbName"]] = int(dbResult["Count"])
                except:
                    result["DbName"] = 0
        result["assembly"] = len(
            Entrez.read(Entrez.esearch(db="assembly", term=search + "[All Names] OR " + search + "[All Uids]"))[
                'IdList'])
        return result

    def search_database(self, database, search, retmax=20, retstart=0):
        records = []
        adapter = NCBISearch.database_map[database]
        esearch = Entrez.read(Entrez.esearch(db=database, term=search, retmax=retmax, retstart=retstart))
        id_list = esearch["IdList"]
        for ncbi_id, data in zip(id_list, adapter.fetch_list(",".join(id_list))):
            record = adapter.adapt(data, ncbi_id)
            records.append(record)
        return records, int(esearch["Count"])
