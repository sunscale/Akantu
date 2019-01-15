from phabricator import Phabricator
from . import export
from .results import Results

def get_phabricator_instance(ctx=None):
    _phab = None
    try:
        _host = None
        _username = None
        _token = None
        if ctx:
            _host = ctx.pop('HOST', None)
            _username = ctx.pop('USERNAME', None)
            _token = ctx.pop('API_TOKEN', None)
            
        _phab = Phabricator(host=_host,
                            username=_username,
                            token=_token)
        _phab.update_interfaces()
        # this request is just to make an actual connection
        _phab.user.whoami()
    except Exception as e:
        _logger.error(
            'Could not connect to phabricator, either give the' +
            ' connection with the default configuration of arc' +
            ' or in the backend configuration of the configuration' +
            ' file:\n' +
            '  in/out:\n' +
            '    username: mylogin\n' +
            '    host: https://c4science.ch/\n' +
            '    token: cli-g3amff25kdpnnv2tqvigmr4omnn7\n')
        raise e
    return _phab

@export
class Harbormaster:
    STATUS = {Results.PASS: 'pass',
              Results.FAIL: 'fail'}
    
    def __init__(self, **kwargs):
        ctx = kwargs['ctx']
        self.__phid = ctx['BUILD_TARGET_PHID']
        self.__phab = get_phabricator_instance(**kwargs)

    def _send_message(self, results):
        self.__phab.harbormaster.sendmessage(buildTargetPHID=self.__phid,
                                             type=self.STATUS[results])

    def send_unit_tests(self, tests):
        _unit_tests = []
        _format = tests.test_format
        for _test in tests:
            _test_dict = {
                'name': _test.name,
                'result': self.STATUS[_test.status],
                'format': _format
            }
            if _test.duration:
                _test_dict['duration'] = _test.duration
            if _test.path:
                _test_dict['path'] = _test.path
                
            _unit_tests.append(_test_dict)

        _msg = {'buildTargetPHID': self.__phid,
                'type': 'work',
                'unit':_unit_tests}
        print(_msg)
        self.__phab.harbormaster.sendmessage(**_msg)

    def attach_uri(self, key, uri, name):
        self.__phab.harbormaster.createartifact(buildTargetPHID=self.__phid,
                                                artifactType='uri',
                                                artifactKey=name,
                                                artifactData={
                                                    'uri': uri,
                                                    'name': name,
                                                    'ui.external': True
                                                })
                                                
        
    def passed(self):
        self._send_message(Results.PASS)

    def failed(self):
        self._send_message(Results.PASS)

    
