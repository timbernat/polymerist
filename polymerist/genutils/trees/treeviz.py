'''Wrappers for printing out tree-like data structures'''

__author__ = 'Timotej Bernat'
__email__ = 'timotej.bernat@colorado.edu'

from typing import Callable, Iterable, TypeVar, Union

from anytree import Node
from anytree.render import RenderTree
from anytree.render import AbstractStyle, AsciiStyle, ContStyle, ContRoundStyle, DoubleStyle

_render_style_aliases : dict[AbstractStyle, list[str]] = {
    AsciiStyle     : ['ASCII', 'Ascii', 'ascii'],
    ContStyle      : ['Cont', 'cont'],
    ContRoundStyle : ['Contround', 'contround', 'cont_round', 'cround', 'round'],
    DoubleStyle    : ['Doublestyle', 'Double', 'double', 'dub'],
}
RENDER_STYLE_MAP : dict[str, AbstractStyle] = {}
for StyleType in AbstractStyle.__subclasses__():
    RENDER_STYLE_MAP[StyleType.__name__] = StyleType
    for alias in _render_style_aliases[StyleType]: # register aliases for convenience
        RENDER_STYLE_MAP[alias] = StyleType


def treestr(root : Node, attr : str='name', style : Union[str, AbstractStyle]=ContStyle(), childiter : Callable[[tuple[Node]], Iterable[Node]]=list, maxlevel : int=None) -> str:
    '''Return a printable string representation of a tree from a root Node, reminiscent of GNU tree'''
    if isinstance(style, str):
        StyleType = RENDER_STYLE_MAP[style] # will raise KeyError if undefined style is provided
        style = StyleType()

    return RenderTree(
        node=root,
        style=style,
        childiter=childiter,
        maxlevel=maxlevel,
    ).by_attr(attr)