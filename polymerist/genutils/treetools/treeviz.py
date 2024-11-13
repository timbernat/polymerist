'''Wrappers for printing out tree-like data structures'''

from typing import Callable, Iterable, TypeVar, Union

from anytree import Node
from anytree.render import RenderTree
from anytree.render import AbstractStyle, AsciiStyle, ContStyle, ContRoundStyle, DoubleStyle

RENDER_STYLE_MAP = {
    StyleType.__name__ : StyleType
        for StyleType in AbstractStyle.__subclasses__()
}


def treestr(root : Node, style : Union[str, AbstractStyle]=ContStyle(), childiter : Callable[[tuple[Node]], Iterable[Node]]=list, maxlevel : int=None) -> str:
    '''Return a printable string representation of a tree from a root Node, reminiscent of GNU tree'''
    if isinstance(style, str):
        StyleType = RENDER_STYLE_MAP[style] # will raise KeyError if undefined style is provided
        style = StyleType()

    return RenderTree(
        node=root,
        style=style,
        childiter=childiter,
        maxlevel=maxlevel,
    )