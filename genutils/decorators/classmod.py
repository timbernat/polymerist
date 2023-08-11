'''Decorators for modifying classes'''

from typing import Callable, Iterable, Optional, TypeVar


C = TypeVar('C') # generic type for a to-be-decorated class

def generate_repr(cls : Optional[C]=None, disp_attrs : Optional[Iterable[str]]=None, lookup_attr : Optional[str]=None):
    '''
    Class decorator for auto-generating __repr__ methods

    By default (i.e. with no arguments), generated repr simply returns the name of the class
    If collection of "disp_attrs" is provided, will display the values of <disp_attrs> for the object instance being represented in series
    If "disp_attrs" is NOT provided but "lookup_attr" is, disp_attrs will be looked up from the modified parent class
    '''
    if disp_attrs is None: 
        disp_attrs = [] # set to empty list to avoid mutable default
    
    def class_decorator(cls : C) -> C:
        '''The actual (argument-free) class decorator'''
        nonlocal disp_attrs # avoids multiple scope issues in disparate use cases (refers to the variable in the outermost scope)
        if not disp_attrs and lookup_attr: # only use lookup if one is explicitly provided and no display attributes are provided
            assert(hasattr(cls, lookup_attr))
            disp_attrs = getattr(cls, lookup_attr) # if a lookup attribute is provided, lookup the attribute names within the class being modified

        def _repr_generic(self) -> str:
            attr_str = ', '.join(f'{attr}={getattr(self, attr)}' for attr in disp_attrs)
            return f'{cls.__name__}({attr_str})'
        setattr(cls, '__repr__', _repr_generic)

        return cls
    
    if cls is None: # null case (i.e. call without parens), return factory call
        return class_decorator
    return class_decorator(cls) # return literal class decorator call
        
def register_subclasses(cls : Optional[C]=None, key_attr : str='__name__', reg_attr : str='subclass_registry') -> Callable[[C], C]:
    '''
    Parametric class decorator for automatically generating a registry of subclasses of a target class
    Binds registry to the "registry" class property in the target class
    
    Subclasses are keyed by lookup of a target attribute <key_attr> in the child classes (by default just the name of the subclass),
    while the resulting registry class property is bound to the <reg_attr> attribute of the parent class
    '''
    def class_decorator(cls : C) -> C:
        '''The actual (argument-free) class decorator'''
        @classmethod # method should be accessibel class-wide
        @property    # make property to allow for dynamic subclassing (generated at runtime, not compile time)
        def _registry(cls : C) -> dict[str, C]:
            return { # Keep a registry of all charger implementations for convenience
                getattr(subclass, key_attr) : subclass
                    for subclass in cls.__subclasses__()
            }
        setattr(cls, reg_attr, _registry) # bind registry class property to target class. TODO : check for registry already present in class
        
        return cls # return back the modified class
    
    if cls is None: # null case (i.e. call without parens), return factory call
        return class_decorator
    return class_decorator(cls) # return literal class decorator call
