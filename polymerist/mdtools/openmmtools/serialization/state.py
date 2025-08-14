'''For handling serialization of OpenMM States'''

__author__ = 'Timotej Bernat'
__email__ = 'timotej.bernat@colorado.edu'

from pathlib import Path
from openmm import Context, State, XmlSerializer

from ....genutils.fileutils.pathutils import assemble_path, allow_string_paths


DEFAULT_STATE_PROPS : dict[str, bool] = {
    'getPositions'  : True,
    'getVelocities' : True,
    'getForces'     : True,
    'getEnergy'     : True,
    'getParameters' : True,
    'getParameterDerivatives' : False,
    'getIntegratorParameters' : False
}

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