
import testing
from testing import value_eq,object_eq,text_eq


def test_import():
    from rmg_analyzer import RmgAnalyzer
#end def test_import



def test_empty_init():
    from rmg_analyzer import RmgAnalyzer

    RmgAnalyzer()
#end def test_empty_init
