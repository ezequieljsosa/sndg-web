# @login_required
def assembly_new(request):
    if request.method == 'POST':
        # create a form instance and populate it with data from the request:
        form = AssemblyForm(request.POST)

        if form.is_valid():
            assembly = form.save()

            return HttpResponseRedirect('/assembly/' + str(assembly.id))


    else:
        form = AssemblyForm()

    return render(request, 'forms/assembly_new.html', {'form': form})


def import_resource(request):
    current_user = request.user
    if request.method == 'POST':
        from .io.NCBISearch import NCBISearch
        search = request.POST["search"]
        ncbi_search = NCBISearch()
        page = Page.from_request(request)
        records, total = ncbi_search.search(search, retmax=page.size, retstart=page.offset)
        existing = [x for x in Assembly.objects.filter(name__in=[x.name for x in records])]

        if current_user:
            collaborated = [x.name for x in existing if
                            hasattr(current_user, "person") and current_user.person in x.collaborators.all()]
        else:
            collaborated = []
        collaborate = [x.name for x in existing if x.name not in collaborated]

        page.set_count(total)
        return render(request, 'forms/import_resource.html',
                      {'resource': "assembly", "search": search, "collaborate": collaborate,
                       "page_obj": page, "results": records, "collaborated": collaborated})
    else:
        search = request.GET.get("search", "")
        return render(request, 'forms/import_resource.html', {
            'resource': "assembly", "search": search})