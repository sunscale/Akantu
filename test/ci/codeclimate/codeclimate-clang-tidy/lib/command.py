import os


class Command:
    """Returns command line arguments by parsing codeclimate config file."""
    def __init__(self, config, file_list):
        self.config = config
        self.file_list = file_list

    def build(self):
        command = ['/usr/src/app/bin/run-clang-tidy',
                   '-clang-tidy-binary',
                   '/usr/bin/clang-tidy']

        if self.config.get('check'):
            command.append(
                '-checks {}'.format(self.config.get('check')))

        if self.config.get('config'):
            command.append(
                '-config {}'.format(self.config.get('project')))

        if self.config.get('header-filter'):
            command.append(
                '-header-filter {}'.format(self.config.get('language')))

        if self.config.get('compilation_database_path'):
            command.append(
                '-p {}'.format(self.config.get('compilation_database_path')))

        include_paths = []
        for file_ in self.file_list:
            include_paths.append(os.path.dirname(file_))
        include_paths = [f'--extra-arg -I{path}' for path in set(include_paths)]

        command.extend(include_paths)
        # command.extend(self.file_list)

        return command
