import os
import json
import sys


class Command:
    """Returns command line arguments by parsing codeclimate config file."""
    def __init__(self, config, file_list):
        self.config = config
        self.file_list = file_list

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
            extra_args = self.config['extra-arg']
            if not isinstance(extra_args, list):
                extra_args = [extra_args]

            if 'compilation-database-path' in self.config:
                for arg in extra_args:
                    command.extend(['-extra-arg', arg])

        # if 'line-filter' in self.config:
        #     command.append(
        #         f'--line-filter={json.dumps(self.config["line-filter"])}')

        # if 'system-headers' in self.config:
        #     command.append('--system-headers')

        if 'compilation-database-path' in self.config:
            command.extend(
                ['-p', self.config['compilation-database-path']])
        else:
            include_paths = []
            for file_ in self.file_list:
                include_paths.append(os.path.dirname(file_))
            include_paths = list(set(include_paths))

            include_flags = ' -I'.join(include_paths)

            compile_commands = []
            for file_ in self.file_list:
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

        command.extend(self.file_list)

        return command
