'''
Test that repeat unit SMARTS expansion complies with specification in
our publication (https://pubs.acs.org/doi/10.1021/acs.jcim.3c01691)
'''

__author__ = 'Timotej Bernat'
__email__ = 'timotej.bernat@colorado.edu'

import pytest

from polymerist.polymers.monomers.specification import compliant_mol_SMARTS


# NOTE: making hydrogens explicit here to enforce atom ordering (and implicitly check that spec converter preserves this order)
@pytest.mark.parametrize(
    'smiles,spec_smarts_expected',
    [
    # examples from Figure 4 of our JCIM paper
        ( ## PLA middle repeat unit
            '*OC(=O)C([H])([H])*',
            '[*:1]-[#8D2+0:2]-[#6D3+0:3](=[#8D1+0:4])-[#6D4+0:5](-[#1D1+0:6])(-[#1D1+0:7])-[*:8]'
        ), 
        ( ## terminal glycine (deprotonated)
            '*N([H])C([H])([H])C(=O)[O-]',
            '[*:1]-[#7D3+0:2](-[#1D1+0:3])-[#6D4+0:4](-[#1D1+0:5])(-[#1D1+0:6])-[#6D3+0:7](=[#8D1+0:8])-[#8D1-1:9]'
        ),
        ( ## amidoamide dendrimer core
            '*N(*)C([H])([H])C([H])([H])N(*)*',
            '[*:1]-[#7D3+0:2](-[*:3])-[#6D4+0:4](-[#1D1+0:5])(-[#1D1+0:6])-[#6D4+0:7](-[#1D1+0:8])(-[#1D1+0:9])-[#7D3+0:10](-[*:11])-[*:12]'
        ), 
    # other examples featuring particular edge cases
        ( ## polyacrylamide middle repeat unit
            '*C(C(=O)[NH2])C*',
            '[*:1]-[#6D4+0:2](-[#6D3+0:3](=[#8D1+0:4])-[#7D3+0:5](-[#1D1+0:9])-[#1D1+0:10])(-[#6D4+0:6](-[*:7])(-[#1D1+0:11])-[#1D1+0:12])-[#1D1+0:8]'
        ),
        # TODO: PAAm w/ explicit Hs?
        ( ## styrene monomer (requires kekulization)
            'c1([H])c([H])c([H])c([H])c([H])c1C([H])=C([H])([H])',
            '[#6D3+0:1]1(-[#1D1+0:2])=[#6D3+0:3](-[#1D1+0:4])-[#6D3+0:5](-[#1D1+0:6])=[#6D3+0:7](-[#1D1+0:8])-[#6D3+0:9](-[#1D1+0:10])=[#6D3+0:11]-1-[#6D3+0:12](-[#1D1+0:13])=[#6D3+0:14](-[#1D1+0:15])-[#1D1+0:16]'
        ),
        ( ## styrene monomer, without explicit hydrogens (atom labels will be different)
            'c1ccccc1C=C',
            '[#6D3+0:1]1(-[#1D1+0:9])=[#6D3+0:2](-[#1D1+0:10])-[#6D3+0:3](-[#1D1+0:11])=[#6D3+0:4](-[#1D1+0:12])-[#6D3+0:5](-[#1D1+0:13])=[#6D3+0:6]-1-[#6D3+0:7](=[#6D3+0:8](-[#1D1+0:15])-[#1D1+0:16])-[#1D1+0:14]'
        ),
        ( ## polythiophene middle repeat unit
            '*c1c([H])c([H])c(s1)*',
            '[*:1]-[#6D3+0:2]1=[#6D3+0:3](-[#1D1+0:4])-[#6D3+0:5](-[#1D1+0:6])=[#6D3+0:7](-[*:9])-[#16D2+0:8]-1'
        ),
        ( ## polythiophene middle repeat unit, without explicit hydrogens (atom labels will be different)
            '*c1ccc(s1)*',
            '[*:1]-[#6D3+0:2]1=[#6D3+0:3](-[#1D1+0:8])-[#6D3+0:4](-[#1D1+0:9])=[#6D3+0:5](-[*:7])-[#16D2+0:6]-1'
        ),
        ( ## ordinary isocyanate group (test higher-order linker)
            '[1*]N=C=O',
            '[1*:1]-[#7D2+0:2]=[#6D2+0:3]=[#8D1+0:4]'
        ),
        ( # hypervalent isocyanate group (test higher-order linker AND charges)
            '[1*]=[N+]=C=O',
            '[1*:1]=[#7D2+1:2]=[#6D2+0:3]=[#8D1+0:4]'
        ),
        ( ## PVC middle repeat unit with "weird" isotopes (including on linkers as labels)
            '[5*][13C]([H])([2H])C([37Cl])([H])[6*]',
            '[5*:1]-[13#6D4+0:2](-[#1D1+0:3])(-[2#1D1+0:4])-[#6D4+0:5](-[37#17D1+0:6])(-[#1D1+0:7])-[6*:8]'
        ),
        ( ## chlorinated variant of amino sulfobetaine (many charged atoms)
            '[Cl]C([N+]([H])([H])([H]))(C(=O)[O-])[S+]([O-])([O-])=O',
            '[#17D1+0:1]-[#6D4+0:2](-[#7D4+1:3](-[#1D1+0:4])(-[#1D1+0:5])-[#1D1+0:6])(-[#6D3+0:7](=[#8D1+0:8])-[#8D1-1:9])-[#16D4+1:10](-[#8D1-1:11])(-[#8D1-1:12])=[#8D1+0:13]'
        ), 
        ( ## calcium carbonate (has higher-order valence)
            'C([O-])([O-])=O.[Ca++]',
            '[#6D3+0:1](-[#8D1-1:2])(-[#8D1-1:3])=[#8D1+0:4].[#20D0+2:5]'
        ),
        ( ## calcium carbonate, with notational variant (has higher-order valence)
            'C([O-])([O-])=O.[Ca+2]',
            '[#6D3+0:1](-[#8D1-1:2])(-[#8D1-1:3])=[#8D1+0:4].[#20D0+2:5]'
        ),
        ( ## test that complex SMARTS queries (i.e. not in the SMILES subset of SMARTS) are also read correctly
            '[!$([#6]=O)]-O-[H]',
            '[*:1]-[#8D2+0:2]-[#1D1+0:3]'
        )
    ]
)
def test_spec_compliance(smiles : str, spec_smarts_expected : str) -> None:
    '''Test that "plain" SMILES strings are correctly upgraded to repeat unit template
    specification-compliant SMARTS strings with correct atom ordering and kekulization'''
    assert compliant_mol_SMARTS(smiles) == spec_smarts_expected
    
@pytest.mark.parametrize(
    'smiles',
    [
        'c1ccccc1C(=O)O',
        '[Cl]C([N+]([H])([H])([H]))(C(=O)[O-])[S+]([O-])([O-])=O',
        '*C(C)C(=O)O*',
    ]
)
def test_compliant_spec_SMARTS_idempotence(smiles : str) -> None:
    '''Test that compliant_mol_SMARTS() is idempotent, i.e. passing in an already-compliant SMARTS string returns the same string'''
    spec_smarts = compliant_mol_SMARTS(smiles)
    spec_smiles_secondary = compliant_mol_SMARTS(spec_smarts)
    assert spec_smiles_secondary == spec_smarts