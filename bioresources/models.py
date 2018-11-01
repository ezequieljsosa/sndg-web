# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.contrib.auth.models import User
from django.db import models
from biosql.models import Taxon,Term



def get_class(kls):
    parts = kls.split('.')
    module = ".".join(parts[:-1])
    m = __import__(module)
    for comp in parts[1:]:
        m = getattr(m, comp)
    return m


class ProcessStatus(models.Model):
    name = models.CharField(max_length=200)
    created_at = models.DateTimeField(auto_now_add=True)

    def step(self, key):
        for step in self.steps.all():
            if step.name == key:
                return step
        raise IndexError("'%s' step not found" % key)

    def __str__(self):
        return "ProcessStatus(%s)" % self.name

    def __repr__(self):
        return self.__str__()

class ProcessStatusStep(models.Model):
    name = models.CharField(max_length=200)
    class_name = models.CharField(max_length=200)
    class_identifier = models.CharField(max_length=200)
    cast_int_identifier = models.BooleanField(default=True)
    created_at = models.DateTimeField(auto_now_add=True)
    completed = models.BooleanField(default=False)
    process_status = models.ForeignKey(ProcessStatus, related_name="steps", on_delete=models.CASCADE)

    def __contains__(self, key):
        return len(self.units.filter(process_identifier=key))

    def append(self, db_identifier, process_identifier=None):
        if not process_identifier:
            process_identifier = db_identifier
        ProcessStatusStepProcessUnit(process_status_step=self,
                                     db_identifier=db_identifier, process_identifier=process_identifier).save()

    def results(self):
        clazz = get_class(self.class_name)
        for unit in self.units.all():
            if unit.db_identifier:
                identifier = int(unit.db_identifier) if self.cast_int_identifier else unit.db_identifier
                yield clazz.objects.get(**{self.class_identifier: identifier})

    def __str__(self):
        return "ProcessStatusStep(%s)" % self.name

    def __repr__(self):
        return self.__str__()

class ProcessStatusStepProcessUnit(models.Model):
    db_identifier = models.CharField(max_length=30,null=True)
    process_identifier = models.CharField(max_length=50)
    created_at = models.DateTimeField(auto_now_add=True)
    process_status_step = models.ForeignKey(ProcessStatusStep, related_name="units", on_delete=models.CASCADE)

    def __str__(self):
        return "PSStepPUnit('%s','%s')" % (self.db_identifier,self.process_identifier)

    def __repr__(self):
        return self.__str__()


class Person(models.Model):
    surname = models.CharField(max_length=200, blank=False)
    name = models.CharField(max_length=200, default="")
    scopus_id = models.CharField(max_length=200, null=True)
    scopus_names = models.TextField(null=True)
    email = models.EmailField()
    arg_affiliation = models.BooleanField(default=False)

    def rtype(self):
        return "person"

    def organizations(self):
        organizations = []
        for aff in self.affiliations.all():
            for org in aff.organizations.all():
                if org not in organizations:
                    organizations.append(org)
        return organizations

    def related_org_names(self):
        return [x.name for x in self.organizations()]

    def complete_name(self):
        return self.surname + ", " + self.name

    def type(self):
        return "person"

    def __str__(self):
        return self.name + " " + self.surname


class Identity(models.Model):
    person = models.ForeignKey(Person, on_delete=models.CASCADE)
    identifier = models.CharField(max_length=200)
    email = models.EmailField(null=True)
    url = models.URLField(null=True)
    authority = models.CharField(max_length=200)
    created_at = models.DateTimeField(auto_now_add=True)
    updated_at = models.DateTimeField()
    ends = models.DateTimeField()


class Organization(models.Model):
    name = models.CharField(max_length=200)
    url = models.URLField(null=True)
    country = models.CharField(max_length=200, null=True)
    city = models.CharField(max_length=200, null=True)
    scopus_id = models.CharField(max_length=200, null=True)
    scopus_names = models.TextField(null=True)

    def rtype(self):
        return "org"

    def __str__(self):
        return " ".join([x for x in [self.name, "|", self.country, self.city] if x])


class Resource(models.Model):
    PUBLICATION = 'publication'
    BIOPROJECT = 'bioproject'
    SEQUENCE = 'nuccore'
    ASSEMBLY = 'assembly'
    GENOME = 'genome'
    READS = 'ra'
    STRUCTURE = 'structure'
    EXPRESSION = 'expression'
    BARCODE = 'barcode'
    SAMPLE = 'sample'

    facet_dict = {
        "assembly": ["species_name", "level", "assembly_type"],
        "gds": ["pdat", "gdstype"],
        "bioproject": ["sample_scope", "material"],  # , "capture_target", "method"
        "barcode": ["subdivision", "marker"],

    }

    RESOURCE_TYPES = (
        (PUBLICATION, PUBLICATION),
        (BIOPROJECT, BIOPROJECT),
        (SEQUENCE, SEQUENCE),
        (ASSEMBLY, ASSEMBLY),
        (GENOME, GENOME),
        (READS, READS),
        (STRUCTURE, STRUCTURE),
        (EXPRESSION, EXPRESSION),
        (BARCODE, BARCODE),
    )

    id = models.AutoField(primary_key=True)

    type = models.CharField(max_length=10, choices=RESOURCE_TYPES)
    name = models.CharField(max_length=350, blank=False)
    description = models.TextField(blank=True)

    ncbi_tax = models.ForeignKey(Taxon, db_column="ncbi_tax", to_field="ncbi_taxon_id",
                                 on_delete=models.SET_NULL, null=True, related_name="bioresources")

    deprecated = models.BooleanField(default=False)
    index_updated = models.BooleanField(default=False)

    created_at = models.DateTimeField(auto_now_add=True)
    updated_at = models.DateTimeField(auto_now=True)

    def compile(self):
        return """<metadata xmlns=\"http://www.lyncode.com/xoai\">\n    <element name=\"dc\">\n        \n        <element name=\"contributor\">\n            <element name=\"advisor1\">\n                \n            <field name=\"value\">Godoy, Moacir Fernandes de</field>\n         </element>\n            <element name=\"advisor1ID\">\n                \n            <field name=\"value\">CPF:57095256853</field>\n         </element>\n            <element name=\"advisor1Lattes\">\n                \n            <field name=\"value\">http://buscatextual.cnpq.br/buscatextual/visualizacv.do?id=K4787885D5&amp;dataRevisao=null</field>\n         </element>\n            <element name=\"advisor-co1\">\n                \n            <field name=\"value\">Azoubel, Reinaldo</field>\n         </element>\n            <element name=\"advisor-co1ID\">\n                \n            <field name=\"value\">CPF:01542680891</field>\n         </element>\n            <element name=\"advisor-co1Lattes\">\n                \n            <field name=\"value\">http://buscatextual.cnpq.br/buscatextual/visualizacv.do?id=K4780809Y6&amp;dataRevisao=null</field>\n         </element>\n            <element name=\"referee1\">\n                \n            <field name=\"value\">Prates, José Carlos</field>\n         </element>\n            <element name=\"referee1ID\">\n                \n            <field name=\"value\">CPF:00000000004</field>\n         </element>\n            <element name=\"referee1Lattes\">\n                \n            <field name=\"value\">http://buscatextual.cnpq.br/buscatextual/visualizacv.do?id=K4780499D8&amp;dataRevisao=null</field>\n         </element>\n            <element name=\"referee2\">\n                \n            <field name=\"value\">Taboga, Sebastião Roberto</field>\n         </element>\n            <element name=\"referee2ID\">\n                \n            <field name=\"value\">CPF:05643691876</field>\n         </element>\n            <element name=\"referee2Lattes\">\n                \n            <field name=\"value\">http://buscatextual.cnpq.br/buscatextual/visualizacv.do?id=K4785629H9&amp;dataRevisao=null</field>\n         </element>\n            <element name=\"referee3\">\n                \n            <field name=\"value\">Braile, Domingo Marcolino</field>\n         </element>\n            <element name=\"referee3ID\">\n                \n            <field name=\"value\">CPF:01172786887</field>\n         </element>\n            <element name=\"referee3Lattes\">\n                \n            <field name=\"value\">http://buscatextual.cnpq.br/buscatextual/visualizacv.do?id=K4767966J3&amp;dataRevisao=null</field>\n         </element>\n            <element name=\"referee4\">\n                \n            <field name=\"value\">Fattini, Carlo Américo</field>\n         </element>\n            <element name=\"referee4ID\">\n                \n            <field name=\"value\">CPF:00000000003</field>\n         </element>\n        <element name=\"authorID\">\n            <field name=\"value\">CPF:20267167814</field>\n         </element>\n         <element name=\"authorLattes\">\n            <field name=\"value\">http://buscatextual.cnpq.br/buscatextual/visualizacv.do?id=K4735493J9&amp;dataRevisao=null</field>\n         </element>\n         <element name=\"author\">\n            <field name=\"value\">Batigália, Fernando</field>\n         </element>\n      </element>\n        <element name=\"date\">\n            <element name=\"accessioned\">\n                \n            <field name=\"value\">2016-01-26T12:51:07Z</field>\n         </element>\n            <element name=\"available\">\n                \n            <field name=\"value\">2006-09-04</field>\n         </element>\n            <element name=\"issued\">\n                \n            <field name=\"value\">2003-03-14</field>\n         </element>\n        </element>\n        <element name=\"identifier\">\n            <element name=\"citation\">\n                \n            <field name=\"value\">BATIGÁLIA, Fernando. Stereological study of vasa vasorum in coronary arteries with different degrees of atherosclerosis. 2003. 485 f. Tese (Doutorado em Medicina Interna; Medicina e Ciências Correlatas) - Faculdade de Medicina de São José do Rio Preto, São José do Rio Preto, 2003.</field>\n         </element>\n            <element name=\"uri\">\n                \n            <field name=\"value\">http://bdtd.famerp.br/handle/tede/5</field>\n         </element>\n        </element>\n        <element name=\"description\">\n            \n            <element name=\"abstract\">\n                <element name=\"eng\">\n                    <field name=\"value\">Introduction: About half of the cases of atherosclerotic coronary disease (a mean of 15% of women deaths and 25% of men) cannot be explained by most of the known risk factors. Coronary vasa vasorum are associated with coronary artery disease; however, their anatomy and physiopathology are not well clear.\n\nObjective: The aim of this study was to carry out a post mortem stereological study of adventitial vasa vasorum in different histopathological degrees of coronary atherosclerosis intending to correlate vasa vasorum, myocardial infarction physiopathology and histopathological degrees of atherosclerosis.\n\nMethod: Ten consecutive autopsies of adults (5 men, 5 women, from 35 to 83 years-old, frozen at 4o C) were performed. Six proximal, medium and distal biopsies of the anterior and posterior interventricular coronary branches (at intervals of 1.5 cm) were performed per autopsy (a total of 60 coronary biopsy fragments). Fragments were processed by histological routine technique and cut in 4 fragments of 4 mm thickness. The first two consecutive histological fragments were stained by hematoxylin-eosin, and the two remaining by Masson´s trichrome. The fragments were histopathologically analysed according to Stary´s coronary atherosclerosis classification and examined by Zeiss Jenaâ, a light microscope with a bright chamber attached a Zeissâ micrometer scale, to outline adventitial vasa vasorum as well as to measure the coronary intraluminal diameter and the medial thickness. Intersection points of vasa vasorum with Merzâ´s grille were manually counted. For all types of vasa vasorum, points on Merzâ´s grille were counted to obtain the following stereological parameters of vasa vasorum: diameter, wall thickness, volumetric and superficial density, and adventitial connective tissue density. Parametric data were analysed by Pearson s linear correlation and principal component analysis. Agreement in determining coronary atherosclerosis degree in laminas stained by hematoxylin-eosin or Masson´s trichrome was assessed by kappa statistics. Differences among variables at each atherosclerosis degree was assessed by analysis of variance or Kruskal-Wallis test. \n\nResults: Coronary intraluminal diameter correlated negatively with coronary medial thickness and number of adventitial vasa vasorum (r&gt;0.50; P-value&lt;0.05). These correlations may be explained by sex, age and coronary atherosclerosis degree (r&gt;0.50; P-value&lt;0.05). All stereological parameters of vasa vasorum correlated negatively with coronary intraluminal diameter and positively with medial thickness, both explained by sex and atherosclerosis degree (r&gt;0.50; P-value&lt;0.05). The size of all types of vasa vasorum augmented proportionally to atherosclerosis histopathological degree aggravation. Kappa statistics for hematoxylin-eosin and Masson´s trichrome presented agreements varying from  substantial or good  to  almost perfect or fine . All variables presented significant differences since the degree  II  of atherosclerosis. \n\nConclusions: Coronary medial thickness and number of vasa vasorum correlated negatively with coronary intraluminal diameter. These correlations may be explained by sex and coronary atherosclerosis degree. Stereological parameters of vasa vasorum (except coronary adventitial connective tissue density) correlated positively with number of adventitial vasa vasorum as well as medial thickness, and negatively with coronary intraluminal diameter. Both correlations were determined by degree of atherosclerosis. Venular rupture in vulnerable atherosclerotic plaques may be associated with myocardial infarction arising since the degree  II  of atherosclerosis.</field>\n                </element>\n            <element name=\"por\">\n               <field name=\"value\">Introdução: Cerca de 50% dos casos de doença arterial coronária aterosclerótica (15% das mortes masculinas e 25% femininas) não são explicados pelos clássicos fatores de risco. Vasa vasorum coronários podem associar-se à aterosclerose coronariana; contudo, sua anatomia e fisiopatologia não estão completamente elucidadas.\n\nObjetivos: Realizar estudo estereológico dos vasa vasorum da túnica externa de artérias coronárias autopsiadas buscando correlação entre vasa vasorum, fisiopatologia do infarto do miocárdio e graus histopatológicos de aterosclerose. \n\nMaterial e Método: Em dez autópsias consecutivas de adultos (5 homens, 5 mulheres, 35 a 83 anos, congelados a 4º C) efetuaram-se biópsias proximal, média e distal dos ramos interventriculares anterior e posterior em cada autópsia, a cada 1,5 cm, totalizando 6 biópsias por autópsia. Cada fragmento foi processado histologicamente com cortes sucessivos de 4 cm de espessura, totalizando 4 fragmentos (24 fragmentos por autópsia). Os dois primeiros fragmentos foram corados em hematoxilina-eosina e os dois últimos em tricrômico de Masson. Lâminas histológicas foram diagnosticadas histopatologicamente quanto ao grau de aterosclerose coronariana pela classificação de Stary e examinadas em microscópio de luz Zeiss Jena com câmara clara e escala micrométrica Zeissâ para delineamento dos vasa vasorum da túnica externa coronária e mensuração do diâmetro do lúmen coronário e da espessura da túnica média coronária. Pontos dos delineamentos dos vasa vasorum sobre a grade estereológica foram manualmente contados obtendo-se diâmetro do lúmen, espessura da parede e densidades volumétrica e superficial dos vasa vasorum, e densidade do tecido conectivo da túnica externa. Dados paramétricos foram analisados por correlação linear de Pearson e análise de componentes principais, concordância no diagnóstico do grau aterosclerótico pela estatística kappa, e diferenças entre valores das variáveis a cada grau aterosclerótico por análise de variância ou teste de Kruskal-Wallis \n\nResultados: Diâmetro do lúmen coronário correlacionou-se negativamente com espessura da túnica média e número de vasa vasorum da túnica externa (r&gt;0,50; valor-p&lt;0,05), com correlações explicadas pelo sexo e grau de aterosclerose (r&gt;0,50; valor-p&lt;0,05). Parâmetros estereológicos dos vasa vasorum correlacionaram-se negativamente com diâmetro do lúmen coronário e positivamente com espessura da túnica média, com correlações explicadas pelo grau de aterosclerose (r&gt;0,50; valor-p&lt;0,05). Tamanho de todos os tipos de vasa vasorum aumentou proporcionalmente ao agravamento das lesões ateroscleróticas. Estatística kappa para lâminas histológicas coradas em hematoxilina-eosina ou em tricrômico de Masson apresentou concordâncias variando de  substancial ou boa  a  quase perfeita ou ótima . Todas as variáveis envolvidas apresentaram diferenças significativas a partir do grau  II  de aterosclerose coronariana.\n\n\nConclusões: Espessura da túnica média e número de vasa vasorum da túnica externa correlacionaram-se negativamente com diâmetro do lúmen coronário, com correlações explicadas pelo sexo e grau de aterosclerose. Parâmetros estereológicos (exceto densidade do tecido conectivo da túnica externa) variaram proporcionalmente com espessura da túnica média e com número de vasa vasorum, e inversamente com diâmetro do lúmen coronário, com correlações explicadas pelo grau histopatológico de aterosclerose. Tamanho de cada tipo de vasa vasorum aumentou proporcionalmente à progressão dos graus ateroscleróticos.  Ruptura venular precoce em placas ateroscleróticas vulneráveis poderia propiciar infarto do miocárdio desde o grau  II  ou  III  de aterosclerose.</field>\n            </element>\n         </element>\n            <element name=\"provenance\">\n                \n            <field name=\"value\">Made available in DSpace on 2016-01-26T12:51:07Z (GMT). No. of bitstreams: 1\nbatigalia_tese_parte1.pdf: 948858 bytes, checksum: b48f5591911812d8f298686183d84a2a (MD5)\n  Previous issue date: 2003-03-14</field>\n         </element>\n            <element name=\"sponsorship\">\n                \n            <field name=\"value\">Sociedade de Cardiologia do Estado de São Paulo</field>\n         </element>\n        </element>\n        <element name=\"format\">\n            \n        <element name=\"none\">\n            <field name=\"value\">application/pdf</field>\n         </element>\n      </element>\n        <element name=\"language\">\n            \n        <element name=\"iso\">\n            <field name=\"value\">por</field>\n         </element>\n      </element>\n        <element name=\"publisher\">\n            \n            \n            \n            \n            \n        <element name=\"none\">\n            <field name=\"value\">Faculdade de Medicina de São José do Rio Preto</field>\n         </element>\n         <element name=\"program\">\n            <field name=\"value\">Programa de Pós-Graduação em Ciências da Saúde</field>\n         </element>\n         <element name=\"initials\">\n            <field name=\"value\">FAMERP</field>\n         </element>\n         <element name=\"country\">\n            <field name=\"value\">BR</field>\n         </element>\n         <element name=\"department\">\n            <field name=\"value\">Medicina Interna; Medicina e Ciências Correlatas</field>\n         </element>\n      </element>\n        \n        <element name=\"subject\">\n            <element name=\"por\">\n                <field name=\"value\">Coronariopatia</field>\n                <field name=\"value\">Vasa Vasorum</field>\n                <field name=\"value\">Estereologia</field>\n                <field name=\"value\">Túnica Externa</field>\n                <field name=\"value\">Coronária</field>\n                <field name=\"value\">Aterosclerose</field>\n            </element>\n            <element name=\"eng\">\n                <field name=\"value\">Stereology</field>\n                <field name=\"value\">External Layer</field>\n                <field name=\"value\">Adventitia</field>\n                <field name=\"value\">Coronary</field>\n                <field name=\"value\">Atherosclerosis</field>\n            </element>\n            <element name=\"cnpq\">\n                \n            <field name=\"value\">CNPQ::CIENCIAS DA SAUDE::MEDICINA::CLINICA MEDICA::CARDIOLOGIA</field>\n         </element>\n        </element>\n        <element name=\"title\">\n            <element name=\"por\">\n                <field name=\"value\">Estudo estereológico dos vasa vasorum em artérias coronárias com diferentes graus de aterosclerose</field>\n            </element>\n            <element name=\"alternative\">\n                <element name=\"eng\">\n                    <field name=\"value\">Stereological study of vasa vasorum in coronary arteries with different degrees of atherosclerosis</field>\n                </element>\n            </element>\n        </element>\n        <element name=\"type\">\n            \n        <element name=\"status\">\n            <field name=\"value\">info:eu-repo/semantics/publishedVersion</field>\n         </element>\n         <element name=\"driver\">\n            <field name=\"value\">info:eu-repo/semantics/doctoralThesis</field>\n         </element>\n      </element>\n      <element name=\"rights\">\n         <element name=\"driver\">\n            <field name=\"value\">info:eu-repo/semantics/openAccess</field>\n         </element>\n      </element>\n      <element name=\"source\">\n         <element name=\"none\">\n            <field name=\"value\">reponame:Biblioteca Digital de Teses e Dissertações da FAMERP</field>\n            <field name=\"value\">instname:FAMERP</field>\n            <field name=\"value\">instacron:FAMERP</field>\n         </element>\n      </element>\n   </element>\n    <element name=\"bundles\">\n        <element name=\"bundle\">\n            <field name=\"name\">ORIGINAL</field>\n            <element name=\"bitstreams\">\n                <element name=\"bitstream\">\n                    <field name=\"name\">batigalia_tese_parte1.pdf</field>\n                    <field name=\"format\">application/pdf</field>\n                    <field name=\"size\">948858</field>\n                    \n                    <field name=\"checksum\">b48f5591911812d8f298686183d84a2a</field>\n                    <field name=\"checksumAlgorithm\">MD5</field>\n                    <field name=\"sid\">1</field>\n                <field name=\"url\">http://bdtd.famerp.br/bitstream/tede/5/1/batigalia_tese_parte1.pdf</field>\n            </element>\n            </element>\n        </element>\n    </element>\n    <element name=\"others\">\n        <field name=\"handle\">tede/5</field>\n        <field name=\"identifier\">oai:localhost:tede/5</field>\n        <field name=\"lastModifyDate\">2016-01-26 10:51:07.75</field>\n    </element>\n    <element name=\"repository\">\n        <field name=\"name\">TEDE</field>\n        <field name=\"mail\"/>\n    </element>\n</metadata>"""

    class Meta:
        unique_together = ('type', 'name',)

    def ncbi_tax_keywords(self):
        if self.ncbi_tax:
            return self.ncbi_tax.keywords.first().text
        return None

    def taxon_name(self):
        if self.ncbi_tax:
            return [x for x in self.ncbi_tax.names.all() if x.name_class == "scientific name"][0].name
        return None

    def oai_public(self):
        return True

    def oai_collections(self):
        return ["collection1"]

    def oai_communities(self):
        return ["community1"]

    def oai_submitter(self):
        return "sndg"

    def permalink(self):
        return "oai:" + str(self.id)

    def metadata_dc_language(self):
        return ["en"]

    def metadata_dc_rights(self):
        return ["info:eu-repo/semantics/openAccess"]

    def metadata_dc_format(self):
        return ["application/pdf"]

    def metadata_dc_publisher(self):
        return ["Faculdade de Medicina de São José do Rio Preto",
                "Programa de Pós-Graduação em Ciências da Saúde",
                "FAMERP",
                "BR",
                "Medicina Interna; Medicina e Ciências Correlatas"]

    def __str__(self):
        return self.name


class ResourceRelation(models.Model):
    source = models.ForeignKey(Resource, on_delete=models.CASCADE, related_name="targets")
    target = models.ForeignKey(Resource, on_delete=models.CASCADE, related_name="sources")
    role = models.CharField(max_length=200, blank=False)
    deprecated = models.BooleanField(default=False)

    class Meta:
        unique_together = (('source', 'target', 'deprecated'),)

    def __str__(self):
        return ("(" + self.source.type + ":" + str(self.source.id) +
                ") -> " + "(" + self.target.type + ":" + str(self.target.id) + ")")


class BioProject(Resource):
    """
    https://www.ncbi.nlm.nih.gov/books/NBK169438/
    """
    ss_monoisolate = 'monoisolate'
    ss_multispecies = 'multi-species'
    ss_environment = 'environment'
    ss_synthetic = 'synthetic'
    ss_other = 'other'

    SAMPLE_SCOPE_TYPES = (
        (ss_monoisolate, ss_monoisolate),
        (ss_multispecies, ss_multispecies),
        (ss_environment, ss_environment),
        (ss_synthetic, ss_synthetic),
        (ss_other, ss_other),

    )

    m_genome = 'genome'
    m_metagenome = 'metagenome'
    m_chromosome = 'chromosome'
    m_transcriptome = 'transcriptome'
    m_reagent = 'reagent'
    m_proteome = 'proteome'

    MATERIAL_TYPES = (
        (m_genome, m_genome),
        (m_metagenome, m_metagenome),
        (m_chromosome, m_chromosome),
        (m_transcriptome, m_transcriptome),
        (m_reagent, m_reagent),
        (m_proteome, m_proteome),

    )

    c_whole = "whole"
    c_exome = "exome"
    c_barcode = "barcode"
    c_TargetedLocusLoci = "TargetedLocusLoci"

    CAPTURE_TYPES = (
        (c_whole, c_whole),
        (c_exome, c_exome),
        (c_barcode, c_barcode),
        (c_TargetedLocusLoci, c_TargetedLocusLoci),
    )

    sample_scope = models.CharField(max_length=20, choices=SAMPLE_SCOPE_TYPES, null=True)
    material = models.CharField(max_length=20, choices=MATERIAL_TYPES, null=True)
    capture = models.CharField(max_length=200, choices=CAPTURE_TYPES, null=True)

    target = models.CharField(max_length=200, null=True)
    submitters = models.TextField(null=True)

    method = models.CharField(max_length=200, null=True)

    # objetive = models.CharField(max_length=200, blank=False)

    def related_author_names(self):
        ran = []
        for affiliation in self.targets.filter(source__type="pmc").all():
            ran += affiliation.source.author_names()
        return list(set(ran))

    def related_org_names(self):
        ran = []
        for affiliation in self.targets.filter(source__type="pmc").all():
            ran += affiliation.source.affiliation_names()
        return list(set(ran))


class Publication(Resource):
    doi = models.CharField(max_length=100)
    date_of_publication = models.DateField(max_length=200)
    pubmed_id = models.CharField(max_length=50, null=True)
    electronic_id = models.CharField(max_length=50, null=True)
    scopus_id = models.CharField(max_length=50, null=True)
    issn = models.CharField(max_length=50, null=True)

    def affiliation_names(self, all=False):
        affs = []
        for affiliation in self.affiliations.all():
            qs = affiliation.organizations
            if not all:
                qs = qs.filter(country="Argentina")
            for org in qs.all():
                if org.name not in affs:
                    affs.append(org.name)
        return affs

    def author_names(self):
        authors = []
        for affiliation in self.affiliations.filter(organizations__country="Argentina").all():
            if affiliation.author.complete_name() not in authors:
                authors.append(affiliation.author.complete_name())
        return authors

    def author_names_idx(self):
        return " ".join(self.author_names())

    def affiliation_names_idx(self):
        return " ".join(self.affiliation_names())


class Affiliation(models.Model):
    publication = models.ForeignKey(Publication, on_delete=models.CASCADE, related_name="affiliations")
    author = models.ForeignKey(Person, on_delete=models.CASCADE, related_name="affiliations")
    organizations = models.ManyToManyField(Organization)


class ExternalId(models.Model):
    organization = models.ForeignKey(Organization, on_delete=models.CASCADE)
    identifier = models.CharField(max_length=20)
    type = models.CharField(max_length=20)
    resource = models.ForeignKey(Resource, on_delete=models.CASCADE, related_name="external_ids")

    def __str__(self):
        return self.organization.name + ":" + self.identifier


class Structure(Resource):
    pdbClass = models.CharField(max_length=50, null=True)
    deposit_date = models.DateField(null=True)
    method = models.CharField(max_length=50, null=True)
    org_list = models.TextField(null=True)

    def __str__(self):
        return self.name

    def related_author_names(self):
        ran = []
        for affiliation in self.targets.all():
            ran += affiliation.source.author_names()
        return list(set(ran))

    def related_org_names(self):
        ran = []
        for affiliation in self.targets.all():
            ran += affiliation.source.affiliation_names()
        return list(set(ran))


class Assembly(Resource):
    intraspecific_name = models.CharField(max_length=250, null=True)
    species_name = models.CharField(max_length=200, null=True)
    level = models.CharField(max_length=50, null=True,choices=(("complete","complete"),
            ("chromosome","chromosome"),("scaffold","scaffold"),("contig","contig"),))
    ncbi_org = models.CharField(max_length=200, null=True)
    release_date = models.DateField(null=True)
    update_date = models.DateField(null=True)
    assembly_type = models.CharField(max_length=50, null=True, choices=(
        ("haploid", "haploid"), ("diploid", "diploid"), ("other", "other")))

    def related_author_names(self):
        ran = []
        for affiliation in self.targets.all():
            ran += affiliation.source.author_names()
        return list(set(ran))

    def related_org_names(self):
        ran = []
        for affiliation in self.targets.all():
            ran += affiliation.source.affiliation_names()
        return list(set(ran))


class Expression(Resource):
    pdat = models.DateField(null=True)
    gdstype = models.CharField(max_length=250, null=True)
    submitters = models.TextField(null=True)

    def related_author_names(self):
        ran = []
        for affiliation in self.targets.all():
            ran += affiliation.source.author_names()
        return list(set(ran))

    def related_org_names(self):
        ran = []
        for affiliation in self.targets.all():
            ran += affiliation.source.affiliation_names()
        return list(set(ran))


class Barcode(Resource):
    country = models.CharField(max_length=100)
    subdivision = models.CharField(max_length=150)
    marker = models.CharField(max_length=50, null=True)
    image_url = models.URLField(null=True)
    bold_org = models.CharField(max_length=255, null=True)
    collectors = models.CharField(max_length=255, null=True)



class ResourceProperty(models.Model):


    organization = models.ForeignKey(Organization, on_delete=models.DO_NOTHING)
    resource = models.ForeignKey(Resource, on_delete=models.CASCADE)
    term = models.ForeignKey(Term, on_delete=models.DO_NOTHING)

class ResourcePropertyValue(models.Model):
    property = models.ForeignKey(ResourceProperty, on_delete=models.DO_NOTHING)
    value = models.CharField(max_length=200)




class Sample(Resource):
    country = models.CharField(max_length=100)
    subdivision = models.CharField(max_length=150)
    collection_date = models.DateField(null=True)
    publication_date = models.DateField(null=True)
    update_date = models.DateField(null=True)


class ReadsArchive(Resource):
    release_date = models.DateField(null=True)
    update_date = models.DateField(null=True)
