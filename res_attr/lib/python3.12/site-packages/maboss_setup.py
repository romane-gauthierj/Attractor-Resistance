
from ctypes.util import find_library
from colomoto.setup_helper import setup
setup({"pkg": "colomoto/maboss",
    "check_progs": ["MaBoSS"]},
    {"pkg": "colomoto/libsbml-plus-packages",
        "check_install": lambda: find_library("sbml")})
