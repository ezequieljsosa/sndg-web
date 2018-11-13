from django import template
from urllib.parse import urlparse, parse_qsl
from django.conf import settings
import re

register = template.Library()


@register.filter(name='split')
def split(value, sep_idx, ):
    sep, idx1, idx2 = sep_idx.split("|")
    return value.split(sep)[int(idx1):int(idx2)]


@register.filter(name='replace')
def replace(value, org_repl):
    org, repl = org_repl.split("|")
    return value.replace(org, repl)


@register.filter(name='url_without_parameter')
def url_without_parameter(arg1, arg2):
    """
    Removes an argument from the get URL.
    Use:
        {{ request|url_without_parameter:'page' }}

    Args:
        arg1: request
        arg2: parameter to remove
    """
    if arg1.GET.getlist(arg2):
        # get all parameters :
        full_url = urlparse(arg1.get_full_path())
        parameters = {}
        # reconstruct query
        for key, value in parse_qsl(full_url.query):
            if parameters.get(key, None) is None:
                parameters[key] = value
        # remove parameters
        if parameters.get(arg2) is not None:
            del parameters[arg2]
        arg1.META['QUERY_STRING'] = u'&'.join([
            u'{}={}'.format(key, value)
            for key, value in parameters.items()])

    return arg1.get_full_path()

@register.filter(name='qs_without_parameter')
def qs_without_parameter(arg1, arg2):
    """
    Removes an argument from the get URL.
    Use:
        {{ request|url_without_parameter:'page' }}

    Args:
        arg1: request
        arg2: parameter to remove
    """


    parameters = {}
    for key, value in arg1.items():
        if parameters.get(key, None) is None and arg2 != key:
            try:
                parameters[key] = value[0]
            except IndexError:
                parameters[key] = value

    return "&".join(
        [k + "=" + v
         for k,v in  parameters.items() ])




@register.filter(name='index')
def index(List, i):
    return List[int(i)]

numeric_test = re.compile("^\d+$")



@register.filter(name='getattribute')
def getattribute(value, arg):
    """Gets an attribute of an object dynamically from a string name"""

    if hasattr(value, str(arg)):
        return getattr(value, arg)
    elif hasattr(value, '__getitem__') and arg in value:
        return value[arg]
    elif numeric_test.match(str(arg)) and len(value) > int(arg):
        return value[int(arg)]
    else:
        return ""
