from dataclasses import dataclass

@dataclass(frozen=True)
class FnGroupSMARTSEntry:
    '''For encapuslating SMARTS group info from Daylight SMARTS registry'''
    category   : str
    category_desc : str

    group_type : str
    group_name : str

    SMARTS : str
    SMARTS_desc : str
