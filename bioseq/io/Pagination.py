# -*- coding: utf-8 -*-
from __future__ import unicode_literals
import math

class Page:

    @staticmethod
    def from_request(request):
        size = int(request.GET.get("page_size", 10))
        page = int(request.GET.get("page", 1))

        if size > 500:
            size = 500



        return Page(size=size, offset=size * (page - 1))

    def __init__(self, size=10, offset=0, count=None):
        self.count = count
        self.size = size
        self.offset = offset
        self.page = int(1 + offset / size)
        self.number = self.page
        self.set_count(count)

    def set_count(self, count):

        prange = range(1, math.ceil(count / self.size)) if count else []

        self.num_pages = count * 1.0 /  self.size if count else None
        self.previous_page_number = self.page - 1 if self.page > 1 else None
        self.next_page_number = self.page + 1 if self.num_pages and  (self.page < self.num_pages) else None
        self.has_previous = bool(self.previous_page_number)
        self.has_next = bool(self.next_page_number)

        self.paginator = {"show_pages": (prange[max(
            self.page - 3, 0):self.page + 2]),
                          "num_pages": math.ceil(count / self.size) if count else None, "count": count}
