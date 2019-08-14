import xml.etree.ElementTree as xml_etree
from .results import Results
from . import export


@export
class CTestResults:
    STATUS = {'passed': Results.PASS,
              'failed': Results.FAIL}

    def __init__(self, filename):
        self._file = open(filename, "r")
        self._etree = xml_etree.parse(self._file)
        self._root = self._etree.getroot()
        self.test_format = 'CTest'

    def __iter__(self):
        self._tests = iter(self._root.findall('./Testing/Test'))
        return self

    def __next__(self):
        class Test:
            def __init__(self, element):
                self.name = element.find('Name').text
                self.path = element.find('FullName').text
                self.status = CTestResults.STATUS[element.attrib['Status']]
                self.duration = float(element.find(
                    "./Results/NamedMeasurement[@name='Execution Time']/Value").text)
                self.reason = None
                if self.status == Results.FAIL:
                    self.reason = element.find(
                        "./Results/NamedMeasurement[@name='Exit Code']/Value").text
                    if self.reason == "Timeout":
                        self.status = Results.BROKEN
                    else:
                        self.reason = "{0} with exit code [{1}]\nSTDOUT:\n{2}".format(
                            self.reason,
                            element.find(
                                "./Results/NamedMeasurement[@name='Exit Value']/Value").text,
                            '\n'.join((el.text for el in element.findall("./Results/Measurement/Value"))),  # noqa: E501
                        )

        test = next(self._tests)

        return Test(test)

    def __exit__(self, exc_type, exc_value, traceback):
        self._file.close()

    def __enter__(self):
        return self
