#!/usr/bin/env python3
# encoding: utf-8
#
# Inspired from boost's pretty printers from
# Rüdiger Sonderfeld <ruediger@c-plusplus.de>
# and from Pretty-printers for libstc++ from Free Software Foundation, Inc.
#
import gdb
import re
import sys

__use_gdb_pp__ = True
try:
    import gdb.printing
except ImportError:
    __use_gdb_pp__ = False


class AkantuPrinter(object):
    regex = None

    @classmethod
    def supports(cls, typename):
        #print('{0} ~= {1}'.format(typename, cls.regex), file=sys.stderr)
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

        @property
        def enabled(self):
            return True

        def add_printer(self, name, regex, printer):
            self.printers[name] = printer

        def __call__(self, val):
            typename = AkantuPrinter.get_basic_type(val)
            if not typename:
                return None

            for name, printer in self.printers.items():
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
    regex = re.compile('^akantu::Array<(.*?), (true|false)>$')
    name = 'akantu::Array'

    def __init__(self, value):
        self.typename = self.get_basic_type(value)
        self.value = value
        self.ptr = self.value['values']
        self.size = int(self.value['size_'])
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

if sys.version_info[0] > 2:
    ### Python 3 stuff
    Iterator = object
else:
    ### Python 2 stuff
    class Iterator:
        """Compatibility mixin for iterators

        Instead of writing next() methods for iterators, write
        __next__() methods and use this mixin to make them work in
        Python 2 as well as Python 3.

        Idea stolen from the "six" documentation:
        <http://pythonhosted.org/six/#six.Iterator>
        """
        def next(self):
            return self.__next__()


class RbtreeIterator(Iterator):
    """
    Turn an RB-tree-based container (std::map, std::set etc.) into
    a Python iterable object.
    """

    def __init__(self, rbtree):
        self.size = rbtree['_M_t']['_M_impl']['_M_node_count']
        self.node = rbtree['_M_t']['_M_impl']['_M_header']['_M_left']
        self.count = 0

    def __iter__(self):
        return self

    def __len__(self):
        return int (self.size)

    def __next__(self):
        if self.count == self.size:
            raise StopIteration
        result = self.node
        self.count = self.count + 1
        if self.count < self.size:
            # Compute the next node.
            node = self.node
            if node.dereference()['_M_right']:
                node = node.dereference()['_M_right']
                while node.dereference()['_M_left']:
                    node = node.dereference()['_M_left']
            else:
                parent = node.dereference()['_M_parent']
                while node == parent.dereference()['_M_right']:
                    node = parent
                    parent = parent.dereference()['_M_parent']
                if node.dereference()['_M_right'] != parent:
                    node = parent
                    self.node = node
        return result

def get_value_from_Rb_tree_node(node):
    """Returns the value held in an _Rb_tree_node<_Val>"""
    try:
        member = node.type.fields()[1].name
        if member == '_M_value_field':
            # C++03 implementation, node contains the value as a member
            return node['_M_value_field']
        elif member == '_M_storage':
            # C++11 implementation, node stores value in __aligned_membuf
            valtype = node.type.template_argument(0)
            return get_value_from_aligned_membuf(node['_M_storage'], valtype)
    except:
        pass
    raise ValueError("Unsupported implementation for %s" % str(node.type))

@register_pretty_printer
class AkaElementTypeMapArrayPrinter(AkantuPrinter):
    """Pretty printer for akantu::ElementTypeMap<Array<T>>"""
    regex = re.compile('^akantu::ElementTypeMap<akantu::Array<(.*?), (true|false)>\*, akantu::(.*?)>$')
    name = 'akantu::ElementTypeMapArray'

    # Turn an RbtreeIterator into a pretty-print iterator.
    class _iter(Iterator):
        def __init__(self, rbiter, type):
            self.rbiter = rbiter
            self.count = 0
            self.type = type

        def __iter__(self):
            return self

        def __next__(self):
            if self.count % 2 == 0:
                n = next(self.rbiter)
                n = n.cast(self.type).dereference()
                n = get_value_from_Rb_tree_node(n)
                self.pair = n
                item = n['first']
            else:
                item = self.pair['second']
            result = ('[{0}]'.format(self.count), item.dereference())
            self.count = self.count + 1
            return result

    def __init__(self, value):
        self.typename = self.get_basic_type(value)
        self.value = value
        self.data = self.value['data']
        self.ghost_data = self.value['ghost_data']

    def to_string(self):
        m = self.regex.search(self.typename)
        return 'ElementTypMapArray<{0}> with {1} _not_ghost and {2} _ghost'.format(
            m.group(1), len(RbtreeIterator(self.data)), len(RbtreeIterator(self.ghost_data)))

    def children(self):
        m = self.regex.search(self.typename)
        try:
            _type = gdb.lookup_type("akantu::Array<{0}, {1}>".format(
                m.group(1), m.group(2))).strip_typedefs().pointer()
            yield ('[_not_ghost]', self._iter(RbtreeIterator(self.data), _type))
            yield ('[_ghost]', self._iter(RbtreeIterator(self.ghost_data), _type))
        except RuntimeError:
            pass

    def display_hint(self):
        return 'map'


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
