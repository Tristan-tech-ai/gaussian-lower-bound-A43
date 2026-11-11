"""Where the scripts of this folder find their input files: the folder computations/data of the repository.
find(name) returns one input file; state(name) returns one stored optimal state, named as in count_law_data.json, or None."""
import os

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.abspath(os.path.join(HERE, "..", "data"))
FIGSTYLE_DIRS = [HERE]


def find(name):
    p = os.path.join(DATA, name)
    if not os.path.exists(p):
        raise FileNotFoundError("%s is not in %s" % (name, DATA))
    return p


def state(name):
    p = os.path.join(DATA, "stored-optima", os.path.basename(name))
    return p if os.path.exists(p) else None
