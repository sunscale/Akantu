import json
import subprocess
import sys
import re
import os
try:
    from termcolor import colored
except ImportError:
    def colored(text, color):
        return text

from command import Command
from issue_formatter import IssueFormatter
from workspace import Workspace


class Runner:
    CONFIG_FILE_PATH = '/config.json'

    """Runs clang-tidy, collects and reports results."""
    def __init__(self):
        self._config_file_path = self.CONFIG_FILE_PATH
        self._config = {}
        self._decode_config()
        self._ansi_escape = re.compile(r'\x1B(?:[@-Z\\-_]|\[[0-?]*[ -/]*[@-~])')
        self._issue_parse = re.compile(r'(?P<file>.*\.(cc|hh)):(?P<line>[0-9]+):(?P<column>[0-9]+): (warning|error): (?P<detail>.*) \[(?P<type>.*)\]')  # noqa

        self._issues_fpr = []

        self._workspace = Workspace(self._config.get('include_paths', []))
        self._files = self._workspace.files
        self._include_paths = self._workspace.include_paths

    def run(self):
        if not len(self._files) > 0:
            return

        self._print_debug(f'[clang-tidy] analyzing {len(self._files)} files')
        command = Command(self._config, self._workspace).build()

        self._print_debug(f'[clang-tidy] command: {command}')

        self._generate_issues(command)

    def _decode_config(self):
        self._print_debug(f"Decoding config file {self._config_file_path}")

        contents = ""
        with open(self._config_file_path, "r") as config:
            contents = config.read()

        self._config = json.loads(contents)
        self._print_debug(f'[clang-tidy] config: {self._config}')

    def _print_issue(self, issue):
        issue_ = IssueFormatter(issue).format()

        path = os.path.dirname(os.path.abspath(issue_["location"]["path"]))
        if path not in self._include_paths:
            return

        if issue_['fingerprint'] in self._issues_fpr:
            return

        self._issues_fpr.append(issue_['fingerprint'])
        print('{}\0'.format(json.dumps(issue_)))

    def _generate_issues(self, command):
        issue = None
        for line in self._run_command(command):
            clean_line = self._ansi_escape.sub('', line)
            match = self._issue_parse.match(clean_line)
            if match:
                if issue is not None:
                    self._print_issue(issue)
                issue = match.groupdict()
            elif issue:
                if 'content' in issue:
                    issue['content'].append(line)
                else:
                    issue['content'] = [line]
        self._print_issue(issue)

    def _run_command(self, command):
        popen = subprocess.Popen(command,
                                 stdout=subprocess.PIPE,
                                 universal_newlines=True)

        for stdout_line in iter(popen.stdout.readline, ""):
            self._print_debug(stdout_line)
            yield stdout_line

        popen.stdout.close()

        return_code = popen.wait()
        if return_code:
            self._print_debug(
                f"[clang-tidy] {command} ReturnCode {return_code}")
            # raise subprocess.CalledProcessError(return_code, command)

    def _print_debug(self, message):
        print(message, file=sys.stderr, flush=True)
