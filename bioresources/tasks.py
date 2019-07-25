# -*- coding: utf-8 -*-
from __future__ import absolute_import, unicode_literals
import subprocess as sp
import os
from celery import shared_task


import traceback
from datetime import datetime


@shared_task
def execute_job(jobId):
    from .models.Job import Job
    if not os.path.exists("/tmp/sndg"): #TODO hacerlo configurable
        os.makedirs("/tmp/sndg")
    job = Job.objects.get(id=jobId)
    job.status = Job.STATUS.QUEUED
    job.start = datetime.now()
    job.save()
    stderr = "/tmp/sndg/{jid}.err".format(jid=jobId)
    stdout = "/tmp/sndg/{jid}.out".format(jid=jobId)
    try:
        with open(stdout,"w") as hstdout,open(stderr,"w") as hstderr:
           sp.call(job.command, shell=True, stdout=hstdout,
                stderr=hstderr)
    except:
        job = Job.objects.get(id=jobId)
        job.result = open(stderr).read()
        job.status = Job.STATUS.ERROR
        job.dev_error = traceback.format_exc()
    else:
        job = Job.objects.get(id=jobId)
        job.result = stdout
        job.status = Job.STATUS.FINISHED
    job.end = datetime.now()
    job.save()
