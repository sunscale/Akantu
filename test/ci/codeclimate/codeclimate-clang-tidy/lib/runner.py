import json
import subprocess
import sys
import re
import tempfile

from command import Command
from issue_formatter import IssueFormatter
from workspace import Workspace


class Runner:
    CONFIG_FILE_PATH = '/config.json'

    """Runs clang-tidy, collects and reports results."""
    def __init__(self):
        self._config_file_path = self.CONFIG_FILE_PATH
        pass

    def run(self):
        config = self._decode_config()
        self._print_debug(f'[clang-tidy] config: {config}')

        workspace = Workspace(config.get('include_paths'))
        workspace_files = workspace.calculate()

        if not len(workspace_files) > 0:
            return

        self._print_debug(f'[clang-tidy] analyzing {len(workspace_files)} files')

        plugin_config = config.get('config', {})
        command = Command(plugin_config, workspace_files).build()

        self._print_debug(f'[clang-tidy] command: {command}')

        issues = {}
        for file_ in workspace_files:
            results = self._run_command([*command, file_])
            issues_ = self._parse_results(results)
            issues.update(issues_)

        for issue in issues.values():
            # cppcheck will emit issues for files outside of the workspace,
            # like header files. This ensures that we only emit issues for
            # files in the workspace.
            if issue and workspace.should_include(issue["location"]["path"]):
                print('{}\0'.format(json.dumps(issue)))

    def _decode_config(self):
        self._print_debug(f"Decoding config file {self._config_file_path}")

        contents = ""
        with open(self._config_file_path, "r") as config:
            contents = config.read()

        return json.loads(contents)

    def _run_command(self, command):
        process = subprocess.Popen(command,
                                   stdout=subprocess.PIPE,
                                   stderr=subprocess.PIPE)

        stdout, stderr = process.communicate()

        return_vals = stdout.decode('utf-8').split('\n')
        return return_vals

    def _parse_results(self, results):
        ansi_escape = re.compile(r'''
        \x1B  # ESC
        (?:   # 7-bit C1 Fe (except CSI)
            [@-Z\\-_]
        |     # or [ for CSI, followed by a control sequence
            \[
            [0-?]*  # Parameter bytes
            [ -/]*  # Intermediate bytes
            [@-~]   # Final byte
        )''', re.VERBOSE)

        re_log_parse = re.compile(
            r'(?P<file>.*\.(cc|hh)):(?P<line>[0-9]+):(?P<column>[0-9]+): warning: (?P<detail>.*) \[(?P<type>.*)\]'  # noqa
        )
        issues = {}
        issue = None

        for line in results:
            clean_line = ansi_escape.sub('', line)
            match = re_log_parse.match(clean_line)
            if match:
                if issue:
                    issue_ = IssueFormatter(issue).format()
                    issues[issue_['fingerprint']] = issue_
                issue = match.groupdict()
            elif issue:
                if 'content' in issue:
                    issue['content'].append(line)
                else:
                    issue['content'] = [line]
            else:
                issue = None

        if issue:
            issue_ = IssueFormatter(issue).format()
            issues[issue_['fingerprint']] = issue_

        return issues

    def _print_debug(self, message):
        print(message, file=sys.stderr)
