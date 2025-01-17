'''For querying chemical databases for information about molecules specified by SMILES string and other structures'''

__author__ = 'Timotej Bernat'
__email__ = 'timotej.bernat@colorado.edu'

import logging
LOGGER  = logging.getLogger(__name__)

from typing import Any, ClassVar, Container, Iterable, Optional, Sequence
from abc import ABC, abstractmethod

import requests

from ..genutils.decorators.classmod import register_abstract_class_attrs, register_subclasses
from ..genutils.importutils.dependencies import requires_modules, MissingPrerequisitePackage


# CUSTOM EXCEPTIONS
class InvalidPropertyError(Exception):
    '''Raised when attempting to query a property that a chemical database service cannot provide'''
    pass

class NullPropertyResponse(Exception):
    '''Raised when a chemical database query doesn't fail BUT returns a NoneType where not allowed'''
    pass

class ChemicalDataQueryFailed(Exception):
    '''Raised when a chemical data query is unfulfilled by a service'''
    pass

# STRATEGIES BASE FOR QUERYING CHEMICAL DATA
@register_subclasses(key_attr='service_name')
@register_abstract_class_attrs('service_name')
class ChemDBServiceQueryStrategy(ABC):
    '''Implementation of queries from a particular chemical database'''
    @abstractmethod
    def _get_property(self, prop_name : str, representation : str, **kwargs) -> Optional[Any]:
        ...
        
    @classmethod
    def dependencies(cls) -> Iterable[str]:
        '''For internals, allows dynamic checking for package dependencies (useful for automating unit test boilerplate)'''
        ...
        
    @classmethod
    @abstractmethod
    def is_online(cls) -> bool:
        '''Check if the service being queried is online and can accept requests'''
        ...
        
    @classmethod
    @abstractmethod
    def queryable_properties(cls) -> Container[str]:
        '''List which chemical property names can be queried from the service'''
        ...

    @classmethod
    @abstractmethod
    def queryable_namespaces(cls) -> Container[str]:
        '''List which chemical identification types can be searched through by the service'''
        ...
        
    def validate_property(self, prop_name : str) -> None:
        '''Pre-check to ensure that a property is queryable from a service before attempting HTTP query'''
        if prop_name not in self.queryable_properties():
            prop_options_str = '\n'.join(sorted(self.queryable_properties()))
            prop_error_msg = f'Cannot query property "{prop_name}" from {self.service_name}'
            LOGGER.error(prop_error_msg) # log briefer error message in cases where the ensuing ValueError is bypassed
            
            raise InvalidPropertyError(f'{prop_error_msg};\nChoose from one of the following property names:\n{prop_options_str}')
        
    def get_property(
            self, 
            prop_name : str, 
            representation : str, 
            namespace : Optional[str],
            keep_first_only : bool=True,
            allow_null_return : bool=False,
            **kwargs
        ) -> Optional[Any]:
        '''Fetch a property associated with a molecule from a chemical database query service'''
        LOGGER.info(f'Sent query request for property "{prop_name}" to {self.service_name}')
        self.validate_property(prop_name=prop_name)
        
        prop_val = self._get_property(prop_name=prop_name, representation=representation, namespace=namespace, **kwargs)
        if not prop_val:
            prop_val = None # cast empty lists, strings, etc to NoneType
        
        if isinstance(prop_val, Container) and not isinstance(prop_val, str) and keep_first_only: # avoid bug where first char of string response is returned
            prop_val = prop_val[0]
        
        if (prop_val is None) and (not allow_null_return): # NOTE: duplicated NoneType check is needed to catch empty containers which are cast to None above
            null_error_msg = f'{self.service_name} returned NoneType "{prop_name}", which is declared invalid by call signature'
            LOGGER.error(null_error_msg)
            
            raise NullPropertyResponse(null_error_msg)
        LOGGER.info(f'Successfully received property "{prop_name}" from {self.service_name}')
                
        return prop_val
    
# CONCRETE IMPLEMENTATIONS OF CHEMICAL DATABASE SERVICE QUERIES
## NIH CACTUS
cirpy_error = MissingPrerequisitePackage(
    importing_package_name=__spec__.name,
    use_case='Querying the NIH CACTUS Chemical Identifier Resolver (CIR)',
    install_link='https://cirpy.readthedocs.io/en/latest/guide/install.html',
    dependency_name='cirpy',
    dependency_name_formal='CIRpy',
)

class NIHCACTUSQueryStrategy(ChemDBServiceQueryStrategy):
    '''
    Implementation of chemical query requests to the NIH's CADD group 
    Cheminformatics Tools and User Services (CACTUS) Chemical Identifier Resolver (CIR)
    '''
    service_name : ClassVar[str] = 'NIH CACTUS CIR'
    
    @classmethod
    def dependencies(cls):
        return ['cirpy']
    
    @classmethod
    @requires_modules('cirpy', missing_module_error=cirpy_error)
    def queryable_properties(cls) -> set[str]:
        import cirpy 
        
        _CIR_PROPS = {  # see official docs for more info: https://cactus.nci.nih.gov/chemical/structure_documentation
            'smiles',
            'ficts',
            'ficus',
            'uuuuu',
            'hashisy',
            'names',
            'iupac_name',
            'cas',
            'chemspider_id',
            'image',
            'twirl',
            'mw',
            'formula',
            'h_bond_donor_count',
            'h_bond_acceptor_count',
            'h_bond_center_count',
            'rule_of_5_violation_count',
            'rotor_count',
            'effective_rotor_count',
            'ring_count',
            'ringsys_count',
            'inchi',
            'inchikey',
            # shortened aliases of InChI-related properties
            'stdinchi', 
            'stdinchikey',
            # these were not documented on CACTUS or by cirpy, but scraped from webchem: https://github.com/ropensci/webchem/blob/master/R/cir.R#L168-L174
            'deprotonable_group_count',
            'heavy_atom_count',
            'heteroatom_count',
            'hydrogen_atom_count',
            'monoisotopic_mass',
            'protonable_group_count',
            'xlogp2',
        }
        return _CIR_PROPS | cirpy.FILE_FORMATS
    
    @classmethod
    def is_online(cls):
        response = requests.head('https://cactus.nci.nih.gov/chemical/structure')
        return response.status_code < 500 # NOTE: could also be more stringent and check == 200 for OK; enough to just check server-side error for now
    
    @classmethod
    def queryable_namespaces(cls) -> set[str]:
        return { # obtained from https://cirpy.readthedocs.io/en/latest/guide/resolvers.html
            'smiles',
            'stdinchikey',
            'stdinchi',
            'ncicadd_identifier', # (for FICTS, FICuS, uuuuu)
            'hashisy',
            'cas_number',
            'name', # this is not documented but DOES work
            'name_by_opsin',
            'name_by_cir',
        }
    
    @requires_modules('cirpy', missing_module_error=cirpy_error)
    def _get_property(self, prop_name : str, representation : str, namespace : Optional[str]=None, **kwargs):
        import cirpy
        
        return cirpy.resolve(representation, prop_name, resolvers=[namespace], **kwargs)

## PubChem
pubchempy_error = MissingPrerequisitePackage(
    importing_package_name=__spec__.name,
    use_case='Querying the PubChem Compound database',
    install_link='https://pubchempy.readthedocs.io/en/latest/guide/install.html',
    dependency_name='pubchempy',
    dependency_name_formal='PubChemPy',
)

class PubChemQueryStrategy(ChemDBServiceQueryStrategy):
    '''
    Implementation of chemical query requests to PubChem via the
    PUG REST API (https://pubchem.ncbi.nlm.nih.gov/docs/pug-rest)
    '''
    service_name : ClassVar[str] = 'PubChem'
    
    @classmethod
    def dependencies(cls):
        return ['pubchempy']
    
    @classmethod
    def is_online(cls):
        response = requests.get('https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/aspirin/property/IUPACName/TXT') # sample query which is well-formatted
        return response.status_code < 500 # NOTE: enough to just check server-side error for now, but could be more stringent and check if ==200
        
    @classmethod
    @requires_modules('pubchempy', missing_module_error=pubchempy_error)
    def queryable_properties(cls) -> set[str]:
        from pubchempy import PROPERTY_MAP
        
        return set(PROPERTY_MAP.keys()) | set(PROPERTY_MAP.values()) | {'Fingerprint2D'} # also taken from webchem: https://github.com/ropensci/webchem/blob/master/R/pubchem.R#L377C21-L392C55
    
    @classmethod
    def queryable_namespaces(cls) -> set[str]:
        return { # obtained from https://pubchem.ncbi.nlm.nih.gov/docs/pug-rest#section=Input
            'cid',
            'name',
            'smiles',
            'inchi',
            'sdf',
            'inchikey',
            'formula',
            'listkey',
        }
    
    @requires_modules('pubchempy', missing_module_error=pubchempy_error)
    def _get_property(self, prop_name : str, representation : str, namespace : Optional[str]='smiles', **kwargs) -> Optional[Any]:
        from pubchempy import PROPERTY_MAP, get_properties, PubChemPyError
        
        official_prop_name = PROPERTY_MAP.get(prop_name, prop_name) # this is done internally, but needed here to extract the property value from the final return dict
        try:
            pubchem_result = get_properties(official_prop_name, identifier=representation, namespace=namespace, **kwargs)
        except PubChemPyError:
            raise requests.HTTPError # discards some information in return for making Strategy interface oblivious to pubchempy (i.e. in case it is not installed)
        else:
            if pubchem_result:
                pubchem_result = [
                    query_result[official_prop_name] # extract property value from extraneous CID (and any other) info
                        for query_result in pubchem_result
                            if official_prop_name in query_result # skip if return doesn't contain the info we specifically requested (happens occasionally for some reason)
                ] 
            return pubchem_result
        
# UTILITY FUNCTIONS EMPLOYING GENERIC STRATEG(Y/IES)
def get_chemical_property(
        prop_name : str, 
        representation : str, 
        namespace : str='smiles',
        keep_first_only : bool=True,
        allow_null_return : bool=False,
        fail_quietly : bool=False,
        services : Optional[Sequence['ChemDBServiceQueryStrategy']]=None,
        **kwargs,
    ) -> Optional[Any]:
    '''Attempt to fetch a molecular property from a variety of chemical database services, either
    provided manually (in the order they should be checked) or ALL implemented service queries by default
    
    Will return the first valid returned result or, if all services fail, raise Exception
    '''
    # determine services which should be queried
    if services is None:
        services = [chem_query_strat_type() for chem_query_strat_type in ChemDBServiceQueryStrategy.subclass_registry.values()]
    if not services: # check if "services" turns out to be an empty collection (either as-passed or because no subclasses are implemented when defaulting)
        raise IndexError('Must provide at least one chemical database querying strategy to "services"')
    n_services_to_try : int = len(services)
    
    # query services sequentially in order of appearance
    for i, service in enumerate(services, start=1):
        ## validate type of service strategies
        if isinstance(service, type):
            service = service() # allows ChemDBServiceQueryStrategy types to be passed in lieu of instances
        if not isinstance(service, ChemDBServiceQueryStrategy):
            raise TypeError(f'Services must be specified as {ChemDBServiceQueryStrategy.__name__} instances, not objects of type {type(service.__name)}')
        
        ## attempt to query result from service
        LOGGER.info(f'Attempting chemical property query to service {i}/{n_services_to_try} ("{service.service_name}"):')
        try:
            prop_val = service.get_property(
                prop_name,
                representation,
                namespace,
                keep_first_only=keep_first_only,
                allow_null_return=allow_null_return,
                **kwargs,
            )
            return prop_val
        except requests.HTTPError:
            LOGGER.error(f'Query to {service.service_name} failed, either due to connection timeout or invalid request')
            continue
        except (InvalidPropertyError, NullPropertyResponse): # skip over invalid property names (keep trying other services rather than failing)
            # log messages baken in to respective raises for these custom exceptions
            continue
    else: # take action when None of the provided services turn up fruitful
        fail_msg = 'Query could not be fulfilled by any of the provided chemical query services'
        if fail_quietly:
            LOGGER.error(f'{fail_msg}; returning NoneType')
            return None
        else: # fail vocally if none of the services can fulfill the property request
            raise ChemicalDataQueryFailed(fail_msg)