from Bio import Entrez

from .adapters import NCBIAssemblyAdapter


class NCBISearch():

    def search(self,search,retmax=20,retstart=0):
        records = []
        esearch = Entrez.read(Entrez.esearch(db="assembly", term=search, retmax=retmax,retstart=retstart))
        adapter = NCBIAssemblyAdapter()

        for data in adapter.fetch_list(",".join(esearch["IdList"])):
            record = adapter.adapt(data)
            records.append(record)
        return records,int(esearch["Count"])