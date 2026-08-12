'''Unit tests for `chemdbqueries` package'''

__author__ = 'Timotej Bernat'
__email__ = 'timotej.bernat@colorado.edu'

import pytest
from _pytest.mark import ParameterSet

from typing import Any, Callable, Optional, TypeVar
from dataclasses import dataclass, asdict
T = TypeVar('T')

from requests import HTTPError

from polymerist.genutils.importutils.dependencies import modules_installed, MissingPrerequisitePackage
from polymerist.smileslib.chemdbqueries import (
    get_chemical_property,
    InvalidPropertyError,
    NullPropertyResponse,
    ChemicalDataQueryFailed,
    ChemDBServiceQueryStrategy,
    # NOTE: these strategies are defined even if their dependencies aren't installed
    # their object instances will just raise Exception on most of their method calls
    NIHCACTUSQueryStrategy,
    PubChemQueryStrategy,
    get_chemical_property,
    
)

TIMEOUT_ERR_CODES : set[int] = {
    500, # SERVER-SIDE PROBLEM
    502, # INVALID GATEWAY/PROXY RESPONSE UPSTREAM
    503, # MAINTENANCE OR OVERLOAD
    504, # REQUEST TIMEOUT
}
def skip_on_server_errors(func : Callable[[], T]) -> T:
    '''Boilerplate for skipping surrounding tests on server failures upon query'''
    try:
        output = func()
    except HTTPError as exc:
        http_err_code : Optional[int] = getattr(exc, 'code', None)
        if http_err_code in TIMEOUT_ERR_CODES:
            pytest.skip('Overloaded server or request denied; test will be inconclusive')
        else:
            raise exc
    else:
        return output

CHEMDB_STRATEGY_ONLINE : dict[type[ChemDBServiceQueryStrategy], bool] = {}
CHEMDB_STRATEGY_DEPENDENCIES_MET : dict[type[ChemDBServiceQueryStrategy], bool] = {}
for ChemDBStrategy in ChemDBServiceQueryStrategy.__subclasses__(): # dynamically determine criteria for which services should be tested
    CHEMDB_STRATEGY_ONLINE[          ChemDBStrategy] = ChemDBStrategy.is_online()
    CHEMDB_STRATEGY_DEPENDENCIES_MET[ChemDBStrategy] = modules_installed(*ChemDBStrategy.dependencies())
    
def skip_pytest_on_invalid_service(service_type : type[ChemDBServiceQueryStrategy]) -> None:
    '''
    Boilerplate function for skipping a test if the requested database
    query service is either missing local dependencies or is offline
    '''
    if not CHEMDB_STRATEGY_DEPENDENCIES_MET[service_type]:
        pytest.skip(f'{service_type.service_name} is missing Python dependencies')
        
    if not CHEMDB_STRATEGY_ONLINE[service_type]:
        pytest.skip(f'{service_type.service_name} cannot be connected to')

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
class ChemDBQueryParameters:
    '''For encapsulating the many parameters passable to a chemical database service query'''
    identifier : str
    namespace : str
    keep_first_only : bool
    allow_null_return : bool

def ethanol_params() -> ChemDBQueryParameters:
    '''Fixed chemical query inputs for ethanol'''
    return ChemDBQueryParameters(
        identifier='CCO',
        namespace='smiles',
        keep_first_only=True,
        allow_null_return=True, 
    )

def fixed_parameter_examples(chem_input : ChemDBQueryParameters) -> list[
    tuple[str, type[ChemDBServiceQueryStrategy], ChemDBQueryParameters]
]:
    '''Examples which test all queryable properties against a fixed chemical input'''
    return [
        (property_name, strategy_type, chem_input)
            for strategy_type, dependencies_met in CHEMDB_STRATEGY_DEPENDENCIES_MET.items()
                if dependencies_met
                    for property_name in strategy_type.queryable_properties()
    ]

# examples which test that many diverse inputs yield expected outputs
def varied_parameter_examples() -> list[
    tuple[str, type[ChemDBServiceQueryStrategy], ChemDBQueryParameters, Any] | ParameterSet
]:
    return [
        # for NIH CACTUS
        ( ## simple queries known to work for all services
            'iupac_name',
            NIHCACTUSQueryStrategy,
            ChemDBQueryParameters( 
                identifier='CCO',
                namespace='smiles',
                keep_first_only=True,
                allow_null_return=False,
            ),
            'ethanol'
        ),
        (
            'inchi',
            NIHCACTUSQueryStrategy,
            ChemDBQueryParameters( 
                identifier='N-methylformamide',
                namespace='name',
                keep_first_only=True,
                allow_null_return=False,
            ),
            'InChI=1/C2H5NO/c1-3-2-4/h2H,1H3,(H,3,4)/f/h3H',
        ),
        ( ## testing that different namespaces can be queried
            'mw',
            NIHCACTUSQueryStrategy,
            ChemDBQueryParameters(
                identifier='benzophenone',
                namespace='name', 
                keep_first_only=True,
                allow_null_return=False,
            ),
            '182.2214'
        ),
        ( ## testing that returns with multiple data values work
            'names',
            NIHCACTUSQueryStrategy,
            ChemDBQueryParameters(
                identifier='c1ccccc1-C(=S)S',
                namespace='smiles', 
                keep_first_only=False, 
                allow_null_return=False,
            ),
            ['Benzenecarbodithioic acid', '121-68-6', 'UPCMLD00WV-104', 'EINECS 204-491-4', 'Benzenecarbodithioic acid', 'Dithiobenzoic acid', 'NSC732246']
        ),
        ( ## testing that enabling and disabling None returns is handled properly in both cases
            'inchi',
            NIHCACTUSQueryStrategy,
            ChemDBQueryParameters(
                identifier='bogus-name', # this is obviously fake and should not return anything
                namespace='name', 
                keep_first_only=True,
                allow_null_return=True,
            ),
            None
        ),
        pytest.param( 
            'inchi',
            NIHCACTUSQueryStrategy,
            ChemDBQueryParameters(
                identifier='bogus-name', # this is obviously fake and should not return anything
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
            'in_no_way_a_valid_property', # this should not even be considered a valid property
            NIHCACTUSQueryStrategy,
            ChemDBQueryParameters(
                identifier='benophenone', 
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
            'iupac_name',
            PubChemQueryStrategy,
            ChemDBQueryParameters(
                identifier='CCO',
                namespace='smiles',
                keep_first_only=True,
                allow_null_return=False
            ),
            'ethanol'
        ),
        (
            'inchi',
            PubChemQueryStrategy,
            ChemDBQueryParameters( 
                identifier='N-methylformamide',
                namespace='name',
                keep_first_only=True,
                allow_null_return=False,
            ),
            'InChI=1S/C2H5NO/c1-3-2-4/h2H,1H3,(H,3,4)',
        ),
        ( ## testing that different namespaces can be queried
            'MolecularWeight',
            PubChemQueryStrategy,
            ChemDBQueryParameters(
                identifier='InChI=1S/C2H5NO/c1-3-2-4/h2H,1H3,(H,3,4)',
                namespace='inchi', 
                keep_first_only=True,
                allow_null_return=False,
            ),
            '59.07',
        ),
        ( ## testing that returns with multiple data values work
            'HeavyAtomCount',
            PubChemQueryStrategy,
            ChemDBQueryParameters(
                identifier='CCO', 
                namespace='smiles', 
                keep_first_only=False,
                allow_null_return=False,
            ),
            [3], # note that this is wrapped in a list, as are all PubChem queries by deualt; I couldn't find a good example which returns more than one value like cirpy does
        ),
        pytest.param( ## testing sending malformed queries to PubChem
            'inchi',
            PubChemQueryStrategy,
            ChemDBQueryParameters(
                identifier='bogus-name', # this is obviously fake and should not return anything
                namespace='smiles', 
                keep_first_only=True,
                allow_null_return=True,
            ),
            None,
            marks=pytest.mark.xfail(
                raises=(ChemicalDataQueryFailed, HTTPError), # DEVNOTE: HTTPError here is absolutely necessary! Request should return HTTP error 400
                reason='Invalid request sent to PubChem (queried a name as a SMILES string)',
                strict=True,
            )
        ),
        ( ## testing that enabling and disabling None returns is handled properly in both cases
            'inchi',
            PubChemQueryStrategy,
            ChemDBQueryParameters(
                identifier='bogus-name', # this is obviously fake and should not return anything
                namespace='name', 
                keep_first_only=True,
                allow_null_return=True,
            ),
            None,
        ),
        pytest.param(  
            'inchi',
            PubChemQueryStrategy,
            ChemDBQueryParameters(
                identifier='bogus-name', # this is obviously fake and should not return anything
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
            'in_no_way_a_valid_property', # this should not even be considered a valid property
            PubChemQueryStrategy,
            ChemDBQueryParameters(
                identifier='benophenone', 
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
    
class TestChemicalDatabaseServiceQueries:
    @pytest.mark.parametrize(
        'property_name,service_type,query_params',
        fixed_parameter_examples(ethanol_params()),
    )
    def test_queryable_properties(
        self,
        property_name : str,
        service_type : type[ChemDBServiceQueryStrategy],
        query_params : ChemDBQueryParameters,
    ) -> None:
        '''Test that the properties each service type lists as queryable do indeed return valid results'''
        skip_pytest_on_invalid_service(service_type=service_type)
        service = service_type()

        # no assert, simply checking that this doesn't raise uncaught Exception
        _ = skip_on_server_errors(
            lambda : service.get_property(
                property_name=property_name,
                **asdict(query_params)
            ) 
        )

    @pytest.mark.parametrize(
        'property_name,service_type,query_params,expected_property',
        varied_parameter_examples(),
    )
    def test_direct_service_property_query(
        self,
        property_name : str,
        service_type : type[ChemDBServiceQueryStrategy],
        query_params : ChemDBQueryParameters,
        expected_property : Any,
    ) -> None:
        '''Test if a chemical database query through a given service is executed completely and returns the expected result'''
        skip_pytest_on_invalid_service(service_type=service_type)
        service = service_type()

        actual_property = skip_on_server_errors(
            lambda : service.get_property(
                property_name=property_name,
                **asdict(query_params),
            )
        )
        assert actual_property == expected_property
        
    @pytest.mark.parametrize(
        'property_name,service_type,query_params,expected_property',
        varied_parameter_examples(),
    )
    def test_get_chemical_property_wrapper(
        self,
        property_name : str,
        service_type : type[ChemDBServiceQueryStrategy],
        query_params : ChemDBQueryParameters,
        expected_property : Any,
    ) -> None:
        '''Test that requests filtered through the get_chemical_properties() strategy wrapper are executed faithfully'''
        skip_pytest_on_invalid_service(service_type=service_type)

        actual_property = skip_on_server_errors(
            lambda : get_chemical_property(
                property_name,
                **asdict(query_params),
                services=[service_type],
                fail_quietly=False, # CRUCIAL that fail_quietly be False; rely on exceptions to match with xfails
            )
        )
        assert actual_property == expected_property
    