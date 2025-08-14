'''For handling serialization of OpenMM States'''

__author__ = 'Timotej Bernat'
__email__ = 'timotej.bernat@colorado.edu'

import logging
LOGGER = logging.getLogger(__name__)

from typing import Optional, Union

from pathlib import Path
from openmm import Context, State, XmlSerializer

from ....genutils.fileutils.pathutils import allow_string_paths


StateLike = Union[str, Path, State]
DEFAULT_STATE_PROPS : dict[str, bool] = {
    'getPositions'  : True,
    'getVelocities' : True,
    'getForces'     : True,
    'getEnergy'     : True,
    'getParameters' : True,
    'getParameterDerivatives' : False,
    'getIntegratorParameters' : False
}

def load_state_flexible(state : Optional[StateLike]=None) -> Optional[State]:
    '''Allows one to flexibly load an OpenMM state, either from a State object or file-like object'''
    if isinstance(state, State) or (state is None):
        state = state
    else:
        if isinstance(state, Path):
            state_path = state
        elif isinstance(state, str):
            state_path = Path(state)
        # TODO : add support for load from opened file
        else:
            raise TypeError('State can only be loaded from pathlike object') 
        
        try:
            with state_path.open('r') as state_file:
                LOGGER.info(f'Attempting to load State from file "{state_path}"')
                state = XmlSerializer.deserialize(state_file.read())
        except ValueError:
            state = None
    
    if state is None:
        LOGGER.warning('No valid State/State file provided, initializing State as None')
    else:
        LOGGER.info(f'Using successfully-initialized State {type(state)}')
    return state

@allow_string_paths
def serialize_state_from_context(
        state_path : Path,
        context : Context,
        state_params : dict[str, bool]=None,
    ) -> None:
    '''For saving State data within an existing OpenMM Context to file'''
    if state_params is None:
        state_params = DEFAULT_STATE_PROPS

    state = context.getState(**state_params)
    with state_path.open('w') as file:
        file.write(XmlSerializer.serialize(state))

def apply_state_to_context(context : Context, state : State) -> None: # TOSELF : this might be replaced with Context.getState() followed by context.reinitialize(preserveState=True)
    '''For applying saved State data to an existing OpenMM Context'''
    context.setPeriodicBoxVectors(*state.getPeriodicBoxVectors())
    context.setPositions(state.getPositions())
    context.setVelocities(state.getVelocities())
    context.setTime(state.getTime())

    context.reinitialize(preserveState=True)