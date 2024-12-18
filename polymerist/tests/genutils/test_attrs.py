'''Unit tests for `attrs` package'''

__author__ = 'Timotej Bernat'
__email__ = 'timotej.bernat@colorado.edu'

from polymerist.genutils.attrs import compile_argfree_getable_attrs
import pytest


class ArgfreeGettableAttrTest():
    '''Dummy class for testing that dynamic attribute inspection works properly'''
    FOO = 'bar'

    def __init__(self, answer : int=42, spam : str='eggs') -> None:
        self.answer = answer
        self.spam = spam

    def get_answer(self) -> int:
        '''Prototypical getter which returns an attribute of self and is prefixed by "get"'''
        return self.answer
    
    def stringy_answer(self) -> str:
        '''A getter whose name does NOT contain any form of "get", but which nevertheless qualifies as argument-free'''
        return str(self.answer)
    
    def GetSpam(self) -> str:
        '''Camel-case getter to check that this is correctly handled by regex'''
        return self.spam

    @property
    def get_answer_prop(self) -> str:
        '''Property getter, to test that inspection excludes this'''
        return self.answer
    
    @classmethod
    def get_foo(cls) -> str:
        '''Class attr-based getter to test how shared-namespace attrs are handled'''
        return cls.FOO

    def _get_spam(self) -> int:
        '''Getter which is marked "private" and should not be returned '''
        return self.spam
    
    def get_echo(self, val : str) -> str:
        '''Getter which requires an argument passed and is therefore not argument-free'''
        return val

@pytest.fixture
def testobj() -> ArgfreeGettableAttrTest:
    return ArgfreeGettableAttrTest(answer=57, spam='ham')


# NOTE: not checking with getter_re='get' as this also returns __getstate__, which is somewhat unexpected to handle

def test_getable_attrs_regex_any(testobj):
    '''Test if getable arg regex match hits anywhere in method name'''
    attrs = compile_argfree_getable_attrs(testobj, getter_re='answer') # regex match should also hit in MIDDLE of string
    assert attrs == {'get_answer': testobj.answer, 'stringy_answer' : str(testobj.answer)}

def test_getable_attrs_regex_start(testobj):
    '''Test if getable arg regex match hits only at the start of in method name'''
    attrs = compile_argfree_getable_attrs(testobj, getter_re='^get') # regex match should also hit in MIDDLE of string
    assert attrs == {'get_answer': testobj.answer, 'get_foo' : testobj.FOO}

def test_getable_attrs_regex_case_insensitive(testobj):
    '''Test if getable arg regex match hits find method with regex only at the start of in method name irrespective of case'''
    attrs = compile_argfree_getable_attrs(testobj, getter_re='^[gG]et')#, getter_re='get')
    assert attrs == {'GetSpam': testobj.spam, 'get_answer': testobj.answer, 'get_foo': testobj.FOO}

def test_getable_attrs_regex_camel(testobj):
    '''Test if getable arg regex match hits camel-cased '''
    attrs = compile_argfree_getable_attrs(testobj, getter_re='^Get[A-Z].*') # regex match should also hit in MIDDLE of string
    assert attrs == {'GetSpam' : testobj.spam}

def test_getable_attrs_regex_repl(testobj):
    '''Test if regex replacement on returned method names is performed as expected'''
    repl = 'TEST'
    attrs = compile_argfree_getable_attrs(testobj, getter_re='^get', repl_str=repl) # regex match should also hit in MIDDLE of string
    assert attrs == {f'{repl}_answer': testobj.answer, f'{repl}_foo': testobj.FOO}