import gdb

# Load the pretty-printers.
from .printers import register_akantu_printers
register_akantu_printers(gdb.current_objfile())
