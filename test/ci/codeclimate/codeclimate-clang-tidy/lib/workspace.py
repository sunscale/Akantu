import os


class Workspace:
    def __init__(self, include_paths, suffixes=['.c', '.cpp', '.cc', '.cxx']):
        self._include_paths = include_paths
        self._suffixes = suffixes

    @property
    def files(self):
        paths = self._walk()
        return [path for path in paths
                if self._should_include(path)]

    @property
    def include_paths(self):
        paths = self._walk()
        return [path for path in paths if os.path.isdir(path)]

    @property
    def paths(self):
        return [path for path in self._include_paths
                if self._should_include(path) or os.path.isdir(path)]

    def _should_include(self, name):
        _, ext = os.path.splitext(name)
        return ext in self._suffixes

    def _walk(self):
        paths = []

        for path in self._include_paths:
            if os.path.isdir(path):
                for root, dirs, files in os.walk(path):
                    paths.extend([os.path.join(root, dir_) for dir_ in dirs])
                    paths.extend([os.path.join(root, file_) for file_ in files])
            else:
                paths.append(path)

        return [os.path.abspath(path) for path in paths]
