import sys as __aka_sys


def export(definition):
    """
    Decorator to export definitions from sub-modules to the top-level package

    :param definition: definition to be exported
    :return: definition
    """
    __module = __aka_sys.modules[definition.__module__]
    __pkg = __aka_sys.modules[__module.__package__]
    __pkg.__dict__[definition.__name__] = definition

    if '__all__' not in __pkg.__dict__:
        __pkg.__dict__['__all__'] = []

    __pkg.__all__.append(definition.__name__)

    return definition


try:
    from termcolor import colored
except ImportError:
    # noinspection PyUnusedLocal
    def colored(string, *args, **kwargs):
        return string
__all__ = ['colored']


from . import truss_fe        # NOQA: F401
from . import static_solver   # NOQA: F401
from . import dynamic_solver  # NOQA: F401
