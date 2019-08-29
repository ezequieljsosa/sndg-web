# -*- coding: utf-8 -*-
from __future__ import unicode_literals
import os
import sys
from datetime import datetime
import gzip
import subprocess as sp
from tqdm import tqdm
from glob import glob

from django.utils.translation import gettext_lazy as _, ngettext as __
from django.db import models
from django.conf import settings

import Bio.SeqIO as bpio
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from bioseq.io.BioIO import BioIO
from bioseq.io.DB2JBrowse import DB2JBrowse
from bioresources.io.NCBISearch import NCBISearch
from bioresources.models.Job import Job
from bioresources.models.Assembly import Assembly
from bioseq.io.DB2JBrowse import DB2JBrowse


class LoadGenomeJob(Job):
    assembly = models.ForeignKey(Assembly, on_delete=models.CASCADE, related_name="jobs")

    def execute(self):
        with open(self.result, "w") as stdout, open(self.dev_error, "w") as stderr:
            sys.stdout = stdout
            sys.stderr = stderr

            NCBISearch.download_assembly(self.assembly.name, workdir=self.job_dir())
            input_file = glob(self.job_dir() + "*_genomic.gbff.gz")[0]

            io = BioIO(self.assembly.name, self.assembly.ncbi_tax.ncbi_taxon_id)

            grep_cmd = 'zgrep -c "FEATURES *Location/Qualifiers" "%s"' % input_file
            io.create_db()

            total = int(sp.check_output(grep_cmd, shell=True))
            with gzip.open(input_file, "rt") as h:
                io.process_record_list(bpio.parse(h, "gb"), total)

            io = DB2JBrowse(jbrowse_path=settings.ROOT_DIR, jbrowse_data_path=os.path.abspath(settings.SNDG_JBROWSE))
            io.ovewrite = True
            io.excluded.append("gene")
            io.db2fs(self.assembly.name)

        self.end = datetime.now()
        self.status = Job.STATUS.FINISHED
