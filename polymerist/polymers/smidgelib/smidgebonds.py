'''Representation classes for monomer graphs bonds which encode intermonomer Port information'''

__author__ = 'Timotej Bernat'
__email__ = 'timotej.bernat@colorado.edu'

import re
from rdkit import Chem

from typing import ClassVar, Union, Optional
from dataclasses import dataclass, field

from ...smileslib.primitives import BONDTYPE_BY_BOND_SMARTS, BOND_SMARTS_BY_BONDTYPE, BOND_PRIMITIVES_FOR_REGEX


# PROCESSING BOND TOKENS IN SMIDGE STRINGS
BOND_TOKEN_RE = re.compile(
    r'(?P<incoming_flavor>\d*)' \
    f'(?P<bondtype>{BOND_PRIMITIVES_FOR_REGEX})' \
    r'(?P<outgoing_flavor>\d*)'
)

@dataclass
class MonomerGraphBondInfo:
    '''Encapsulated information about an intermonomer bond in a monomer graph'''
    DEFAULT_BONDTYPE : ClassVar[Chem.BondType] = Chem.BondType.UNSPECIFIED
    DEFAULT_FLAVOR   : ClassVar[int] = 0

    incoming_flavor : int    = field(default=DEFAULT_FLAVOR)
    bondtype : Chem.BondType = field(default=DEFAULT_BONDTYPE)
    outgoing_flavor : int    = field(default=DEFAULT_FLAVOR)

    def __post_init__(self) -> None:
        '''Perform appropriate type conversions and apply defaults'''
        if self.incoming_flavor is None:
            self.incoming_flavor = self.DEFAULT_FLAVOR

        if self.bondtype is None:
            self.bondtype = self.DEFAULT_BONDTYPE

        if self.outgoing_flavor is None:
            self.outgoing_flavor = self.DEFAULT_FLAVOR

    @property
    def bond_str(self) -> str:
        '''SMARTS representation of the current bondtype (defaults to the default symbol if None or invalid bondtype is provided)'''
        return BOND_SMARTS_BY_BONDTYPE.get(self.bondtype, BOND_SMARTS_BY_BONDTYPE.get(self.DEFAULT_BONDTYPE))

    def __str__(self) -> str:
        return f'{self.incoming_flavor or ""}{self.bond_str}{self.outgoing_flavor or ""}'

    
    @staticmethod
    def parse_str_dict(str_dict : dict[str, str]) -> dict[str, Optional[Union[str, Chem.BondType]]]:
        '''Parse string-valued dict into dict of correct types and NoneType defaults for empty keys'''
        # 0) create copy of dict to avoid any in-place modification
        imb_info = {k : v for k, v in str_dict.items()}
        
        # 1) process bond type
        imb_info['bondtype'] = BONDTYPE_BY_BOND_SMARTS.get(imb_info.get('bondtype')) # set to None either if no bondtype is provided OR the provided type is not a registered primitive

        # 2) process port flavors
        for flavor_attr in ('incoming_flavor', 'outgoing_flavor'):
            if not (flavor_str := imb_info.get(flavor_attr)): # raised when the flavor_str is either not present or empty
                imb_info[flavor_attr] = None
            elif isinstance(flavor_str, str) and flavor_str.isdigit():
                imb_info[flavor_attr] = int(flavor_str) # convert parsed strings to ints where possible

        return imb_info

    @classmethod
    def from_dict(cls, str_dict : dict[str, str]) -> 'MonomerGraphBondInfo':
        '''Initialize from dictionary of values, after sanitizing'''
        return cls(**cls.parse_str_dict(str_dict))

    @classmethod
    def from_match(cls, match : Optional[re.Match]) -> 'MonomerGraphBondInfo':
        '''Initialize from groupdict of regex Match'''
        if match is None:
            # TODO: add logged warning
            return cls()
        return cls.from_dict(match.groupdict())