from django.shortcuts import render, redirect, reverse
from bioresources.models.Job import Job


def blast(request):
    ndatabases = [1, 2, 3]
    pdatabases = [4, 56]
    if "biodatabase" in request.GET:
        db = Biodatabase.objects.get(biodatabase_id=request.GET["biodatabase"])
        dbp = Biodatabase.objects.get(name=(db.name + "_prots"))
        ndatabases.append({"value": db.biodatabase_id, "label": "Genome - " + db.name})
        pdatabases.append({"value": dbp.biodatabase_id, "label": "Proteome - " + db.name})

    if request.method == 'POST':
        job = Job(command="ls /data", user=request.user)
        job.save()
        job.execute()
        return redirect(reverse("bioresources:job", kwargs={"jid": job.id}))

    return render(request, 'tools/blast.html', {"ndatabases": ndatabases, "pdatabases": pdatabases})


def blast_result(request, jid):
    job = Job.objects.get(id=jid)
    if job.status == Job.STATUS.FINISHED:
        return render(request, 'tools/blast_result.html', {})
    else:
        return render(request, 'tools/generic_result.html', {})
