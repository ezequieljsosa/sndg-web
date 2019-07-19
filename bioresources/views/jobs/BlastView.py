def blast(request):
    ndatabases = []
    pdatabases = []
    if "biodatabase" in request.GET:
        db = Biodatabase.objects.get(biodatabase_id=request.GET["biodatabase"])
        dbp = Biodatabase.objects.get(name=(db.name + "_prots"))
        ndatabases.append({"value": db.biodatabase_id, "label": "Genome - " + db.name})
        pdatabases.append({"value": dbp.biodatabase_id, "label": "Proteome - " + db.name})

    return render(request, 'biosql/tools/blast.html', {"ndatabases": ndatabases, "pdatabases": pdatabases})


def blast_result(request):
    return render(request, 'biosql/tools/blast_result.html', {})