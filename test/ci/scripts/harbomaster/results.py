from enum import Enum


class Results(Enum):
    PASS = 0
    FAIL = 1
    BROKEN = 2
    SKIP = 3
    UNSTABLE = 4
