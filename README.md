# sndg-web
SNDG Web Apps

/opt/solr-7.3.1/bin/solr start
./manage.py shell_plus --ipython --print-sql


pip install  git+https://github.com/ezequieljsosa/elsapy.git@complete_view_fields


./manage.py  build_solr_schema -u default > /data/sndg_solr/schema.xml

docker run -v $PWD/data/sndg_solr:/opt/solr/server/solr/Notes --name sndg_solr7 -d -p 127.0.0.1:8984:8983 -t solr:7 
docker exec -u root     sndg_solr7 bash -c 'chown 8983:8983 /opt/solr/server/solr/Notes'
curl "http://localhost:8984/solr/admin/cores?action=CREATE&name=Notes&schema=schema.xml"


HAYSTACK_ID_FIELD="item.id" HAYSTACK_DJANGO_CT_FIELD='metadata.django_ct'  HAYSTACK_DJANGO_ID_FIELD='item.id' HAYSTACK_DOCUMENT_FIELD='metadata.search'  ./manage.py rebuild_index -u oai


$ bioresources: ../manage.py makemessages -l es

