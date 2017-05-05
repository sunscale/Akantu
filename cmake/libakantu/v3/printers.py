#!/usr/bin/env python
# encoding: utf-8
#
# Inspired from boost's pretty printers from
# Rüdiger Sonderfeld <ruediger@c-plusplus.de>
# and from Pretty-printers for libstc++ from Free Software Foundation, Inc.
#

import gdb
import re
# import libstdcxx.v6.printers as std

__use_gdb_pp__ = True
try:
    import gdb.printing
except ImportError:
    __use_gdb_pp__ = False


class AkantuPrinter(object):
    regex = None

    @classmethod
    def supports(cls, typename):
        return cls.regex.search(typename)

    @classmethod
    def get_basic_type(cls, value):
        """ Determines the type associated to a value"""

        _type = value.type
        # If it points to a reference, get the reference.
        if _type.code == gdb.TYPE_CODE_REF:
            _type = _type.target()

        # Get the unqualified type, stripped of typedefs.
        _type = _type.unqualified().strip_typedefs()

        return _type.tag

if __use_gdb_pp__:
    __akantu_pretty_printers__ = \
        gdb.printing.RegexpCollectionPrettyPrinter("libakantu-v3")
else:
    class AkantuPrettyPrinters(object):
        def __init__(self, name):
            super(AkantuPrettyPrinters, self).__init__()
            self.name = name
            self.printers = {}

        def add_printer(self, name, regex, printer):
            self.printers[name] = printer

        def __call__(self, val):
            typename = AkantuPrinter.get_basic_type(val)
            if not typename:
                return None

            for name, printer in self.printers.iteritems():
                if(printer.supports(typename)):
                    return printer

            return None

    __akantu_pretty_printers__ = AkantuPrettyPrinters("libakantu-v3")


def register_pretty_printer(pretty_printer):
    "Registers a Pretty Printer"

    __akantu_pretty_printers__.add_printer(pretty_printer.name,
                                           pretty_printer.regex,
                                           pretty_printer)

    return pretty_printer

@register_pretty_printer
class AkaArrayPrinter(AkantuPrinter):
    """Pretty printer for akantu::Array<T>"""
    regex = re.compile('^akantu::Array<(.*), (true|false)>$')
    name = 'akantu::Array'

    def __init__(self, value):
        self.typename = self.get_basic_type(value)
        self.value = value
        self.ptr = self.value['values']
        self.size = int(self.value['size'])
        self.nb_component = int(self.value['nb_component'])

    def display_hint(self):
        return 'array'

    def to_string(self):
        m = self.regex.search(self.typename)
        return 'Array<{0}>({1}, {2}) stored at {3}'.format(
            m.group(1), self.size, self.nb_component, self.ptr)

    def children(self):
        _ptr = self.ptr
        for i in range(self.size):
            _values = ["{0}".format((_ptr + j).dereference())
                       for j in range(self.nb_component)]
            _ptr = _ptr + self.nb_component
            yield ('[{0}]'.format(i),
                   ('{0}' if self.nb_component == 1 else '[{0}]').format(
                       ', '.join(_values)))


# @register_pretty_printer
class AkaTensorPrinter(AkantuPrinter):
    """Pretty printer for akantu::Tensor<T>"""
    regex = re.compile('^akantu::Tensor<(.*), +(.*), +(.*)>$')
    name = 'akantu::Tensor'

    value = None
    typename = ""
    ptr = None
    dims = []
    ndims = 0

    def pretty_print(self):
        def ij2str(i, j, m):
            return "{0}".format((self.ptr+m*j + i).dereference())

        def line(i, m, n):
            return "[{0}]".format(", ".join((ij2str(i, j, m) for j in
                                             range(n))))

        m = int(self.dims[0])
        if (self.ndims == 1):
            n = 1
        else:
            n = int(self.dims[1])

        return "[{0}]".format(", ".join(line(i, m, n) for i in range(m)))

    def __init__(self, value):
        self.typename = self.get_basic_type(value)
        self.value = value
        self.ptr = self.value['values']
        self.dims = self.value['n']

    def children(self):
        yield ('values', self.pretty_print())
        yield ('wrapped', self.value['wrapped'])

@register_pretty_printer
class AkaVectorPrinter(AkaTensorPrinter):
    """Pretty printer for akantu::Vector<T>"""
    regex = re.compile('^akantu::Vector<(.*)>$')
    name = 'akantu::Vector'
    n = 0
    ptr = 0x0

    def __init__(self, value):
        super(AkaVectorPrinter, self).__init__(value)
        self.ndims = 1

    def to_string(self):
        m = self.regex.search(self.typename)
        return 'Vector<{0}>({1}) [{2}]'.format(m.group(1), int(self.dims[0]),
                                               str(self.ptr))


@register_pretty_printer
class AkaMatrixPrinter(AkaTensorPrinter):
    """Pretty printer for akantu::Matrix<T>"""
    regex = re.compile('^akantu::Matrix<(.*)>$')
    name = 'akantu::Matrix'

    def __init__(self, value):
        super(AkaMatrixPrinter, self).__init__(value)
        self.ndims = 2

    def to_string(self):
        m = self.regex.search(self.typename)
        return 'Matrix<%s>(%d, %d) [%s]' % (m.group(1), int(self.dims[0]),
                                            int(self.dims[1]),
                                            str(self.ptr))


@register_pretty_printer
class AkaElementPrinter(AkantuPrinter):
    """Pretty printer for akantu::Element"""
    regex = re.compile('^akantu::Element$')
    name = 'akantu::Element'

    def __init__(self, value):
        self.typename = self.get_basic_type(value)
        self.value = value

        self.element = self.value['element']
        self.eltype = self.value['type']
        self.ghost_type = self.value['ghost_type']

    def to_string(self):
        return 'Element({0}, {1}, {2})'.format(self.element, self.eltype,
                                               self.ghost_type)


def register_akantu_printers(obj):
    "Register Akantu Pretty Printers."

    if __use_gdb_pp__:
        gdb.printing.register_pretty_printer(obj, __akantu_pretty_printers__)
    else:
        if obj is None:
            obj = gdb
        obj.pretty_printers.append(__akantu_pretty_printers__)
