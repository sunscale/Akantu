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

import xml.etree.ElementTree as xml_etree
from collections import namedtuple

from .results import Results
from . import export


class TestResults:
    STATUS = {'passed': Results.PASS,
              'failed': Results.FAIL}
    Test = namedtuple('Test', 'name path status duration reason'.split())

    def __init__(self, filename):
        self._file = open(filename, "r")
        self._etree = xml_etree.parse(self._file)
        self._root = self._etree.getroot()

    def __exit__(self, exc_type, exc_value, traceback):
        self._file.close()

    def __enter__(self):
        return self


@export
class CTestResults(TestResults):
    def __init__(self, filename):
        super().__init__(self, filename)
        self.test_format = 'CTest'

    def __iter__(self):
        self._tests = iter(self._root.findall('./Testing/Test'))
        return self

    def __next__(self):
        element = next(self._tests)

        name = element.find('Name').text
        path = element.find('FullName').text
        status = CTestResults.STATUS[element.attrib['Status']]
        duration = float(element.find(
            "./Results/NamedMeasurement[@name='Execution Time']/Value").text)
        reason = None

        if status == Results.FAIL:
            reason = element.find(
                "./Results/NamedMeasurement[@name='Exit Code']/Value").text
            if reason == "Timeout":
                status = Results.BROKEN
            else:
                reason = "{0} with exit code [{1}]\nSTDOUT:\n{2}".format(
                    reason,
                    element.find(
                        "./Results/NamedMeasurement[@name='Exit Value']/Value").text,  # noqa: E501
                    '\n'.join((el.text for el in element.findall("./Results/Measurement/Value"))),  # noqa: E501
                )

        return self.Test(name, path, status, duration, reason)


@export
class JUnitTestResults(TestResults):
    def __init__(self, filename):
        super().__init__(self, filename)
        self.test_format = 'JUnit'

    def __iter__(self):
        self._tests = iter(self._root.findall('testcase'))
        return self

    def __next__(self):
        element = next(self._tests)

        name = element.attrib['name']
        path = element.attrib['file'] \
            + ':{}'.format(element.attrib['line'])
        duration = element.attrib['time']

        failure = element.find('failure')
        error = element.find('error')

        if failure is not None and error is not None:
            status = Results.PASS

        elif error:
            status = Results.BROKEN
            reason = error.attrib['message'] + '\n' + error.attrib['type']
        elif failure:
            status = Results.FAIL
            reason = failure.attrib['message'] + '\n' \
                + failure.attrib['type']
        return self.Test(name, path, status, duration, reason)
