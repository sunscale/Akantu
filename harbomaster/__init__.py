# ------------------------------------------------------------------------------
__author__ = "Nicolas Richart"
__copyright__ = "Copyright (C) 2019, EPFL (Ecole Polytechnique Fédérale" \
                " de Lausanne) Laboratory (LSMS - Laboratoire de Simulation" \
                " en Mécanique des Solides)"
__credits__ = ["Nicolas Richart", "Lucas Frerot"]
__license__ = "L-GPLv3"
__maintainer__ = "Nicolas Richart"
__email__ = "nicolas.richart@epfl.ch"
# ------------------------------------------------------------------------------

# for the module
import sys as __hbm_sys


def export(definition):
    """
    Decorator to export definitions from sub-modules to the top-level package

    :param definition: definition to be exported
    :return: definition
    """
    __module = __hbm_sys.modules[definition.__module__]
    __pkg = __hbm_sys.modules[__module.__package__]
    __pkg.__dict__[definition.__name__] = definition

    if '__all__' not in __pkg.__dict__:
        __pkg.__dict__['__all__'] = []

    __pkg.__all__.append(definition.__name__)

    return definition


from . import testresults   # noqa
from . import arclint       # noqa
from . import hbm           # noqa
