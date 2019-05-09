import py11_akantu
private_keys = set(dir(py11_akantu)) - set(dir())

for k in private_keys:
    globals()[k] = getattr(py11_akantu, k)


def initialize(*args, **kwargs):
    raise RuntimeError("No need to call initialize,"
                       " use parseInput to read an input file")


def finalize(*args, **kwargs):
    raise RuntimeError("No need to call finalize")


FromStress = py11_akantu.FromHigherDim
FromTraction = py11_akantu.FromSameDim
py11_akantu.__initialize()