#!/usr/bin/env python3
import re
import os
import hashlib
import json
import sys


# 7-bit C1 ANSI sequences
ansi_escape = re.compile(r'''
    \x1B  # ESC
    (?:   # 7-bit C1 Fe (except CSI)
        [@-Z\\-_]
    |     # or [ for CSI, followed by a control sequence
        \[
        [0-?]*  # Parameter bytes
        [ -/]*  # Intermediate bytes
        [@-~]   # Final byte
    )
''', re.VERBOSE)

re_log_parse = re.compile(
    r'(?P<file>.*\.(cc|hh)):(?P<line>[0-9]+):(?P<column>[0-9]+): warning: (?P<detail>.*) \[(?P<type>.*)\]'  # noqa
)


categories = {
    "bugprone": "Bug Risk",
    "modernize": "Clarity",
    "mpi": "Bug Risk",
    "openmp": "Bug Risk",
    "performance": "Performance",
    "readability": "Clarity",
}

issues = {}

with open('clang-tidy-all-out.log', 'r') as log:
    for line in log:
        clean_line = ansi_escape.sub('', line)
        match = re_log_parse.match(clean_line)
        if match:
            line_dict = match.groupdict()
            line_dict['file'] = os.path.relpath(line_dict['file'])

            line_dict['fingerprint'] = hashlib.md5(
                '{file}:{line}:{column}:{type}'.format(**line_dict).encode()
            ).hexdigest()

            issue = {
                'type': 'issue',
                'description': line_dict['detail'],
                'fingerprint': line_dict['fingerprint'],
                'location': {
                    "path": line_dict['file'],
                    "lines": {
                        "begin": int(line_dict['line']),
                        "end": int(line_dict['line']),
                    },
                    "positions": {
                        "begin": {
                            "line": int(line_dict['line']),
                            "column": int(line_dict['column']),
                        },
                    },

                },
                'severity': 'minor',
            }

            category = line_dict['type'].split('-')[0]
            if category in categories:
                issue['category'] = categories[category]

            # use a dictionnary to avoid duplicates
            issues[issue['fingerprint']] = issue

issues = list(issues.values())

json.dump(issues, sys.stdout)
