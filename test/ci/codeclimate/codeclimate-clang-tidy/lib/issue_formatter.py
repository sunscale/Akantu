import hashlib
import os


class IssueFormatter:
    CLASSIFICATIONS = {
        'bugprone': {
            'categories': ['Bug Risk'],
            'severity': 'major',
        },
        'modernize': {
            'categories': ['Clarity', 'Compatibility', 'Style'],
            'severity': 'info'
        },
        'mpi': {
            'categories': ['Bug Risk', 'Performance'],
            'severity': 'critical',
        },
        'openmp': {
            'categories': ['Bug Risk', 'Performance'],
            'severity': 'critical',
        },
        'performance': {
            'categories': ['Performance'],
            'severity': 'minor',
        },
        'readability': {
            'categories': ['Clarity', 'Style'],
            'severity': 'info'
        },
    }

    def __init__(self, issue):
        self.issue_dict = issue

    def format(self):
        self.issue_dict['file'] = os.path.relpath(self.issue_dict['file'])
        issue = {
            'type': 'issue',
            'check_name': self.issue_dict['type'],
            'description': self.issue_dict['detail'],
            'location': {
                "path": self.issue_dict['file'],
                "positions": {
                    "begin": {
                        "line": int(self.issue_dict['line']),
                        "column": int(self.issue_dict['column']),
                    },
                    "end": {
                        "line": int(self.issue_dict['line']),
                        "column": int(self.issue_dict['column']),
                    },
                },
            },
        }

        if 'content' in self.issue_dict:
            issue['content'] = {
                'body': '```\n' +
                '\n'.join(self.issue_dict['content']) +
                '\n```'
            }
        
        issue['fingerprint'] = hashlib.md5(
            '{file}:{line}:{column}:{type}'.format(**self.issue_dict).encode()
        ).hexdigest()

        type_ = self.issue_dict['type'].split('-')[0]
        if type_ in self.CLASSIFICATIONS:
            issue['categories'] = self.CLASSIFICATIONS[type_]['categories']
            issue['severity'] = self.CLASSIFICATIONS[type_]['severity']

        return issue
