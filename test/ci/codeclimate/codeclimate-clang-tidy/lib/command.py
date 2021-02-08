import json
import os
import re


class Command:
    """Returns command line arguments by parsing codeclimate config file."""
    def __init__(self, config, workspace):
        self.config = config
        self._workspace = workspace

    def build(self):
        command = ['/usr/src/app/bin/run-clang-tidy',
                   '-clang-tidy-binary',
                   '/usr/bin/clang-tidy']

        if 'checks' in self.config:
            command.extend(
                ['-checks', self.config["checks"]])

        if 'config' in self.config:
            command.extend(
                ['-config', self.config["config"]])

        if 'header-filter' in self.config:
            command.extend(
                ['-header-filter', self.config["header-file"]])

        extra_args = []
        if 'extra-arg' in self.config:
            tmp_extra_args = self.config['extra-arg']
            if not isinstance(extra_args, list):
                tmp_extra_args = [extra_args]

            extra_args = []
            includes_re = re.compile(r'-I(.*)')
            for arg in tmp_extra_args:
                match = includes_re.match(arg)
                if match:
                    path = os.path.abspath(match.group(1))
                    extra_args.append(f'-I{path}')
                else:
                    extra_args.append(arg)

            if 'compilation-database-path' in self.config:
                for arg in extra_args:
                    command.extend(['-extra-arg', arg])

        if 'compilation-database-path' in self.config:
            command.extend(
                ['-p', self.config['compilation-database-path']])
        else:
            include_flags = ' -I'.join(self._workspace.include_paths)

            compile_commands = []
            for file_ in self._workspace.files:
                cmd = {
                    'directory': os.path.dirname(file_),
                    'file': file_,
                    'command': f'/usr/bin/clang++ {include_flags} {" ".join(extra_args)} -c {file_} -o dummy.o', # noqa
                }
                compile_commands.append(cmd)

            location = '/tmp'
            compile_database = os.path.join(location, 'compile_commands.json')
            with open(compile_database, 'w') as db:
                json.dump(compile_commands, db)

            command.extend(['-p', location])

        command.extend([f'{path}*' for path in self._workspace.paths])

        return command
