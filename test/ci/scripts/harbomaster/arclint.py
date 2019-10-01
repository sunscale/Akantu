# ------------------------------------------------------------------------------
__author__ = "Nicolas Richart"
__copyright__ = "Copyright (C) 2019, EPFL (Ecole Polytechnique Fédérale" \
                " de Lausanne) Laboratory (LSMS - Laboratoire de Simulation" \
                " en Mécanique des Solides)"
__credits__ = ["Nicolas Richart"]
__license__ = "L-GPLv3"
__maintainer__ = "Nicolas Richart"
__email__ = "nicolas.richart@epfl.ch"
# ------------------------------------------------------------------------------

import json
from .results import Results
from . import export


@export
class ARCLintJson:
    STATUS = {'passed': Results.PASS,
              'failed': Results.FAIL}

    def __init__(self, filename):
        self._file = open(filename, "r")
        self._json = json.load(self._file)

    def __iter__(self):
        self._last_path_lints = []
        self._lints = iter(self._json)
        return self

    def __next__(self):
        class Lint:
            def __init__(self, path, json):
                self.json = json
                self.json['path'] = path
                if 'name' in json and 'code' in json and \
                   json['name'] == json['code']:
                    self.json['name'] = json['description']
                    del self.json['description']

            def __getattr__(self, name):
                if name == 'json':
                    return self.json
                elif name in self.json:
                    return self.json[name]
                else:
                    return None

            def __str__(self):
                return f'{self.path} => {self.name}'

        lint = None
        while not lint:
            if len(self._last_path_lints) > 0:
                lint = Lint(self._last_path,
                            self._last_path_lints.pop(0))
                break

            json = next(self._lints)
            if type(json) != dict:
                raise RuntimeError("Wrong input type for the linter processor")

            self._last_path = list(json.keys())[0]
            self._last_path_lints = json[self._last_path]

        return lint

    def __exit__(self, exc_type, exc_value, traceback):
        self._file.close()

    def __enter__(self):
        return self
