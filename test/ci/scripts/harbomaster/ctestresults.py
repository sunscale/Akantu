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
                self.duration = float(element.find("./Results/NamedMeasurement[@name='Execution Time']/Value").text)

            def __str__(self):
                return f'{self._name}: {self._status} in {self._duration}'

        test = next(self._tests)
        #print(test.find('Name'))
        #while not test.find('Name'):
        #    test = next(self._tests)
            
        return Test(test)

    def __exit__(self, exc_type, exc_value, traceback):
        self._file.close()

    def __enter__(self):
        return self
