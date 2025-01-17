'''Unit tests for `chemdbqueries` package'''

__author__ = 'Timotej Bernat'
__email__ = 'timotej.bernat@colorado.edu'

import pytest

from typing import Any
from dataclasses import dataclass, asdict

from requests import HTTPError

from polymerist.genutils.importutils.dependencies import modules_installed, MissingPrerequisitePackage
from polymerist.smileslib.chemdbqueries import (
    get_chemical_property,
    InvalidPropertyError,
    NullPropertyResponse,
    ChemicalDataQueryFailed,
    ChemDBServiceQueryStrategy,
    # NOTE: these strageies are implemented to be defined even if the packages in question aren't installed
    # it's just that instances will raise excception on most of their method calls
    NIHCACTUSQueryStrategy,
    PubChemQueryStrategy,
    get_chemical_property,
    
)

CHEMDB_STRATEGY_ONLINE : dict[str, bool] = {}
CHEMDB_STRATEGY_DEPENDENCIES_MET : dict[ChemDBServiceQueryStrategy, bool] = {}
for ChemDBStrategy in ChemDBServiceQueryStrategy.__subclasses__(): # dynamically determine criteria for which services should be tested
    CHEMDB_STRATEGY_ONLINE[          ChemDBStrategy] = ChemDBStrategy.is_online()
    CHEMDB_STRATEGY_DEPENDENCIES_MET[ChemDBStrategy] = modules_installed(*ChemDBStrategy.dependencies())

@pytest.mark.parametrize(
    'service_type', 
    [
        pytest.param(
            service_type,
            marks=pytest.mark.xfail(
                raises=MissingPrerequisitePackage,
                reason='Unsatisfied dependency needed for chemical database service to be imported',
                strict=True,
            )
        )
            for service_type, dependencies_met in CHEMDB_STRATEGY_DEPENDENCIES_MET.items()
                if not dependencies_met
    ]
)
def test_missing_dependency_xfail(service_type : ChemDBServiceQueryStrategy) -> None:
    '''Test whether checks for missing prerequisite dependencies are active'''
    _ = service_type.queryable_properties() # happen to know this requires respective dependencies at the time of writing


@dataclass
class ChemDBQueryExample:
    '''For encapsulating the many parameters passable to a chemical database service query'''
    prop_name : str
    representation : str
    namespace : str
    keep_first_only : bool
    allow_null_return : bool

# SERVICE_QUERY_EXAMPLES :

EXAMPLES_BY_SERVICE : list[tuple[ChemDBServiceQueryStrategy, ChemDBQueryExample, Any]] = [
    # for NIH CACTUS
    ( ## simple queries known to work for all services
        NIHCACTUSQueryStrategy,
        ChemDBQueryExample( 
            'iupac_name',
            'CCO',
            namespace='smiles',
            keep_first_only=True,
            allow_null_return=False,
        ),
        'ethanol'
    ),
    (
        NIHCACTUSQueryStrategy,
        ChemDBQueryExample( 
            'inchi',
            'N-methylformamide',
            namespace='name',
            keep_first_only=True,
            allow_null_return=False,
        ),
        'InChI=1/C2H5NO/c1-3-2-4/h2H,1H3,(H,3,4)/f/h3H',
    ),
    ( ## testing that different namespaces can be queried
        NIHCACTUSQueryStrategy,
        ChemDBQueryExample(
            'mw',
            'benzophenone',
            namespace='name', 
            keep_first_only=True,
            allow_null_return=False,
        ),
        '182.2214'
    ),
    ( ## testing that returns with multiple data values work
        NIHCACTUSQueryStrategy,
        ChemDBQueryExample(
            'names',
            'c1ccccc1-C(=S)S',
            namespace='smiles', 
            keep_first_only=False, 
            allow_null_return=False,
        ),
        ['Benzenecarbodithioic acid', '121-68-6', 'UPCMLD00WV-104', 'EINECS 204-491-4', 'Benzenecarbodithioic acid', 'Dithiobenzoic acid', 'NSC732246']
    ),
    ( ## testing that enabling and disabling None returns is handled properly in both cases
        NIHCACTUSQueryStrategy,
        ChemDBQueryExample(
            'inchi',
            'bogus-name', # this is obviously fake and should not return anything
            namespace='name', 
            keep_first_only=True,
            allow_null_return=True,
        ),
        None
    ),
    pytest.param( 
        NIHCACTUSQueryStrategy,
        ChemDBQueryExample(
            'inchi',
            'bogus-name', # this is obviously fake and should not return anything
            namespace='name', 
            keep_first_only=True,
            allow_null_return=False,
        ),
        None,
        marks=pytest.mark.xfail(
            raises=(NullPropertyResponse, ChemicalDataQueryFailed),
            reason='Did not allow response to be NoneType',
            strict=True,
        )
    ),
    pytest.param( ## testing that invalid property values are caught before attempting a query
        NIHCACTUSQueryStrategy,
        ChemDBQueryExample(
            'in_no_way_a_valid_property', # this should not even be considered a valid property
            'benophenone', 
            namespace='name', 
            keep_first_only=True,
            allow_null_return=False,
        ),
        None,
        marks=pytest.mark.xfail(
            raises=(InvalidPropertyError, ChemicalDataQueryFailed),
            reason='Tried to query a property that does not exist',
            strict=True,
        )
    ),
    
    # for PubChem
    ( ## simple queries known to work for all services
        PubChemQueryStrategy,
        ChemDBQueryExample(
            'iupac_name',
            'CCO',
            namespace='smiles',
            keep_first_only=True,
            allow_null_return=False
        ),
        'ethanol'
    ),
    (
        PubChemQueryStrategy,
        ChemDBQueryExample( 
            'inchi',
            'N-methylformamide',
            namespace='name',
            keep_first_only=True,
            allow_null_return=False,
        ),
        'InChI=1S/C2H5NO/c1-3-2-4/h2H,1H3,(H,3,4)',
    ),
    ( ## testing that different namespaces can be queried
        PubChemQueryStrategy,
        ChemDBQueryExample(
            'MolecularWeight',
            'InChI=1S/C2H5NO/c1-3-2-4/h2H,1H3,(H,3,4)',
            namespace='inchi', 
            keep_first_only=True,
            allow_null_return=False,
        ),
        '59.07',
    ),
    ( ## testing that returns with multiple data values work
        PubChemQueryStrategy,
        ChemDBQueryExample(
            'HeavyAtomCount',
            'CCO', 
            namespace='smiles', 
            keep_first_only=False,
            allow_null_return=False,
        ),
        [3], # note that this is wrapped in a list, as are all PubChem queries by deualt; I couldn't find a good example which returns more than one value like cirpy does
    ),
    pytest.param( ## testing sending malformed queries to PubChem
        PubChemQueryStrategy,
        ChemDBQueryExample(
            'inchi',
            'bogus-name', # this is obviously fake and should not return anything
            namespace='smiles', 
            keep_first_only=True,
            allow_null_return=True,
        ),
        None,
        marks=pytest.mark.xfail(
            raises=(HTTPError, ChemicalDataQueryFailed),
            reason='Invalid request sent to PubChem (queried a name as a SMILES string)',
            strict=True,
        )
    ),
    ( ## testing that enabling and disabling None returns is handled properly in both cases
        PubChemQueryStrategy,
        ChemDBQueryExample(
            'inchi',
            'bogus-name', # this is obviously fake and should not return anything
            namespace='name', 
            keep_first_only=True,
            allow_null_return=True,
        ),
        None,
    ),
    pytest.param(  
        PubChemQueryStrategy,
        ChemDBQueryExample(
            'inchi',
            'bogus-name', # this is obviously fake and should not return anything
            namespace='name', 
            keep_first_only=True,
            allow_null_return=False,
        ),
        None,
        marks=pytest.mark.xfail(
            raises=(NullPropertyResponse, ChemicalDataQueryFailed),
            reason='Did not allow response to be NoneType',
            strict=True,
        )
    ),
    pytest.param( ## testing that invalid property values are caught before attempting a query
        PubChemQueryStrategy,
        ChemDBQueryExample(
            'in_no_way_a_valid_property', # this should not even be considered a valid property
            'benophenone', 
            namespace='name', 
            keep_first_only=True,
            allow_null_return=False,
        ),
        None,
        marks=pytest.mark.xfail(
            raises=(InvalidPropertyError, ChemicalDataQueryFailed),
            reason='Tried to query a property that does not exist',
            strict=True,
        )
    ),
]
    
@pytest.mark.parametrize('service_type,query_example,expected_return', EXAMPLES_BY_SERVICE)
class TestChemicalDatabaseServiceQueries:
    def test_direct_service_property_query(self, service_type : ChemDBQueryExample, query_example : ChemDBQueryExample, expected_return : Any) -> None:
        '''Test if a chemical database query through a given service is executed completely and returns the expected result'''
        # prechecks to skip tests
        if not CHEMDB_STRATEGY_DEPENDENCIES_MET[service_type]:
            pytest.skip(f'{service_type.service_name} is missing Python dependencies')
            
        if not CHEMDB_STRATEGY_ONLINE[service_type]:
            pytest.skip(f'{service_type.service_name} cannot be connected to')
   
        # initialize and query service
        service = service_type()
        assert service.get_property(**asdict(query_example)) == expected_return
        
    def test_get_chemical_property_wrapper(self, service_type : ChemDBQueryExample, query_example : ChemDBQueryExample, expected_return : Any) -> None:
        '''Test that requests filtered through the get_chemical_properties() strategy wrapper are executed faithfully'''
        # prechecks to skip tests
        if not CHEMDB_STRATEGY_DEPENDENCIES_MET[service_type]:
            pytest.skip(f'{service_type.service_name} is missing Python dependencies')
            
        if not CHEMDB_STRATEGY_ONLINE[service_type]:
            pytest.skip(f'{service_type.service_name} cannot be connected to')
   
        # call get_chemical_property wrapper
        # try:
        print(asdict(query_example))
        assert get_chemical_property(**asdict(query_example), services=[service_type], fail_quietly=False) == expected_return # CRUCIAL that fail_quietly be False; rely on exceptions to match with xfails
        # except ChemicalDataQueryFailed:
            # print('foo')
    

# @pytest.mark.skipif(not CHEMDB_STRATEGIES_TO_TEST, reason=f'No chemical data server(s) is not available')
# def test_get_chemical_property() -> None:
#     a = 4