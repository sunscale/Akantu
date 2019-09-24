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

from phabricator import Phabricator
import yaml
import base64
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
        print('Could not connect to phabricator, either give the' +
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
              Results.FAIL: 'fail',
              Results.BROKEN: 'broken',
              Results.SKIP: 'skip',
              Results.UNSTABLE: 'unsound'}

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
        _list_of_failed = {}
        try:
            _yaml = open(".tests_previous_state", 'r')
            _previously_failed = yaml.load(_yaml)
            if not _previously_failed:
                _previously_failed = {}
            _yaml.close()
        except OSError:
            _previously_failed = {}

        for _test in tests:
            status = self.STATUS[_test.status]
            if (_test.status != Results.PASS and
                (_test.name in _previously_failed and
                 (_previously_failed[_test.name] == self.STATUS[_test.status] or  # noqa: E501
                  _previously_failed[_test.name] == 'unsound'))):
                status = 'unsound'

            _test_dict = {
                'name': _test.name,
                'result': status,
                'format': _format
            }

            if _test.duration:
                _test_dict['duration'] = _test.duration
            if _test.path:
                _test_dict['path'] = _test.path
            if _test.reason:
                _test_dict['details'] = _test.reason

            if status != 'pass':
                _list_of_failed[_test.name] = status
            _unit_tests.append(_test_dict)

        with open(".tests_previous_state", 'w+') as _cache_file:
            yaml.dump(_list_of_failed,
                      _cache_file,
                      default_flow_style=False)

        _msg = {'buildTargetPHID': self.__phid,
                'type': 'work',
                'unit': _unit_tests}
        self.__phab.harbormaster.sendmessage(**_msg)

    def send_lint(self, linter_processor):
        _lints = []
        for lint in linter_processor:
            _lint = {}
            for key in ['code', 'name', 'severity', 'path', 'line',
                        'char', 'description']:
                val = getattr(lint, key)
                if val:
                    _lint[key] = val

            _lints.append(_lint)

        _msg = {'buildTargetPHID': self.__phid,
                'type': 'work',
                'lint': _lints}

        self.__phab.harbormaster.sendmessage(**_msg)

    def send_uri(self, key, uri, name):
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
        self._send_message(Results.FAIL)

    def upload_file(self, filename, name, view_phid=None):
        with open(filename, 'rb') as f:
            data = f.read()
            base64_data = base64.b64encode(data)
            _msg = {
                'data_base64': base64_data.decode('ascii'),
                'name': filename,
                }
            if view_phid:
                _msg['viewPolicy'] = view_phid

            _res = self.__phab.file.upload(**_msg)

            print(f"{name} -> {_res}")
            self.__phab.harbormaster.createartifact(
                buildTargetPHID=self.__phid,
                artifactType='file',
                artifactKey=name,
                artifactData={
                    'filePHID': _res.response
                })
