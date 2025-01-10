'''Tools for drawing and visulaizing RDKit molecules'''

__author__ = 'Timotej Bernat'
__email__ = 'timotej.bernat@colorado.edu'

from typing import Optional, Type, Union

import PIL
from PIL.Image import Image

import matplotlib.pyplot as plt
from matplotlib.colors import Normalize, Colormap

from rdkit.Chem.rdchem import Mol
from rdkit.Chem.Draw import rdMolDraw2D, SimilarityMaps, MolsToGridImage, IPythonConsole

from . import rdprops
from ..graphics import imageutils, plotutils
from ..graphics.named_colors import WHITE
           

# GLOBAL PREFERENCES
def set_rdkdraw_size(dim : int=300, aspect : float=3/2):
    '''Change image size and shape of RDKit Mol images'''
    IPythonConsole.molSize = (int(aspect*dim), dim) # Change IPython image display size
    
def enable_substruct_highlights() -> None:
    '''Turns on highlighting of found substructures when performing substructure matches'''
    IPythonConsole.highlightSubstructs = True

def disable_substruct_highlights() -> None:
    '''Turns off highlighting of found substructures when performing substructure matches'''
    IPythonConsole.highlightSubstructs = False

def enable_kekulized_drawing() -> None:
    '''Turns on automatic kekulization of aromatic bonds before drawing molecules in Jupyter Notebooks'''
    IPythonConsole.kekulizeStructures = True

def disable_kekulized_drawing() -> None:
    '''Turns off automatic kekulization of aromatic bonds before drawing molecules in Jupyter Notebooks'''
    IPythonConsole.kekulizeStructures = False


# SINGLE-MOLECULE DISPLAY OPTIONS
def clear_highlights(rdmol : Mol) -> None:
    '''Removes the highlighted atoms flags from an RDKit Mol if present'''
    if hasattr(rdmol, '__sssAtoms'):
        del rdmol.__sssAtoms


# PLOTTING
def tight_norm_for_rdmol_prop(rdmol : Mol, prop : str, prop_type : Union[Type[int], Type[float]]) -> tuple[Normalize, tuple[int, ...]]:
    '''Generate a matplotlib Normalize object with bounds matching the maxima, minima, and midpoint of a numeric Prop'''
    prop_vals = rdprops.aggregate_atom_prop(rdmol, prop, prop_type=prop_type) # explicitly ensure the property is interpreted as a numerical value
    vmin, vmax = min(prop_vals.values()), max(prop_vals.values())
    
    norm = Normalize(vmin, vmax)
    ticks = (norm.vmin, 0, norm.vmax)

    return norm, ticks


def rdmol_prop_heatmap(rdmol : Mol, prop : str, cmap : Colormap, norm : Optional[Normalize]=None, annotate : bool=False, annotate_precision : int=5, img_size : tuple[int, int]=(1_000, 1_000)) -> Image:
    '''Take a charged RDKit Mol and color atoms based on the magnitude of a particular atomwise property'''
    prop_vals = rdprops.aggregate_atom_prop(rdmol, prop, prop_type=float) # currently only support floats for plotting
    if norm is None:
        norm, ticks = tight_norm_for_rdmol_prop(rdmol, prop, prop_type=float)
    
    if annotate:
        rdprops.annotate_atom_prop(rdmol, prop, annotate_precision=annotate_precision, in_place=True)

    colors = {
        atom_num : cmap(norm(prop_val))
            for atom_num, prop_val in prop_vals.items()
    }

    # generate image of Mol
    draw = rdMolDraw2D.MolDraw2DCairo(*img_size) # or MolDraw2DCairo to get PNGs
    rdMolDraw2D.PrepareAndDrawMolecule(draw, rdmol, highlightAtoms=prop_vals.keys(), highlightAtomColors=colors)
    draw.FinishDrawing()
    img_bytes = draw.GetDrawingText()
    
    return imageutils.img_from_bytes(img_bytes)

def rdmol_prop_heatmap_colorscaled(rdmol : Mol, prop : str, cmap : Colormap=plt.get_cmap('turbo'), norm : Optional[Normalize]=None, ticks : Optional[tuple[float]]=None, cbar_label : str='', orient : Optional[str]=None, **heatmap_args) -> tuple[plt.Figure, plt.Axes]:
    '''Plot a labelled heatmap of the charge differences between 2 structurally identical RDKit Molecules with different partial charges'''
    if norm is None:
        norm, ticks = tight_norm_for_rdmol_prop(rdmol, prop, prop_type=float)

    image = rdmol_prop_heatmap(rdmol, prop=prop, cmap=cmap, norm=norm, **heatmap_args)
    image = imageutils.crop_borders(image, bg_color=WHITE)

    return plotutils.plot_image_with_colorbar(image, cmap, norm, ticks=ticks, label=cbar_label, orient=orient)
