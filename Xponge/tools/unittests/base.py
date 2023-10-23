"""
    This **package** includes unittests of the basic functions
"""
__all__ = ["test_run",
           "test_walk"]

def test_run():
    from Xponge import Xprint
    Xprint("Success!")

def test_walk():
    Xprint("bad!")
