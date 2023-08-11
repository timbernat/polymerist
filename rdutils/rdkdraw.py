'''Tools for drawing and visulaizing RDKit molecules'''

from typing import Optional
from . import rdprops
from .rdtypes import RDMol

import PIL
from PIL.Image import Image

import matplotlib.pyplot as plt
from matplotlib.colors import Normalize, Colormap

from rdkit.Chem.Draw import rdMolDraw2D, SimilarityMaps, MolsToGridImage, IPythonConsole

from ..graphics import imageutils, plotutils
from ..graphics.named_colors import WHITE
           

def set_rdkdraw_size(dim : int, aspect : float=1/1):
    '''Change image size and shape of RDMol images'''
    IPythonConsole.molSize = (int(aspect*dim), dim)   # Change image size

def rdmol_prop_heatmap(rdmol : RDMol, prop : str, cmap : Colormap, norm : Normalize, annotate : bool=False, precision : int=5, img_size : tuple[int, int]=(1_000, 1_000)) -> Image:
    '''Take a charged RDKit Mol and color atoms based on the magnitude of a particular atomwise property'''
    colors, prop_vals, atom_nums = {}, [], []
    for atom in rdmol.GetAtoms():
        atom_num, prop_val = atom.GetIdx(), atom.GetDoubleProp(prop)
        colors[atom_num] = cmap(norm(prop_val))
        atom_nums.append(atom_num)
        prop_vals.append(prop_val)

        if annotate:
            atom.SetProp('atomNote', str(round(prop_val, precision))) # need to convert to string, as double is susceptible to float round display errors (shows all decimal places regardless of rounding)

    draw = rdMolDraw2D.MolDraw2DCairo(*img_size) # or MolDraw2DCairo to get PNGs
    rdMolDraw2D.PrepareAndDrawMolecule(draw, rdmol, highlightAtoms=atom_nums, highlightAtomColors=colors)
    draw.FinishDrawing()
    img_bytes = draw.GetDrawingText()
    
    return imageutils.img_from_bytes(img_bytes)

def rdmol_prop_heatmap_colorscaled(rdmol : RDMol, prop : str, cmap : Colormap=plt.get_cmap('turbo'), cbar_label : str='', orient : Optional[str]=None, **heatmap_args) -> tuple[plt.Figure, plt.Axes]:
    '''Plot a labelled heatmap of the charge differences between 2 structurally identical RDKit Molecules with different partial charges'''
    prop_vals = rdprops.aggregate_atom_prop(rdmol, prop, prop_type=float) # explicitly ensure the property is interpreted as a numericla value
    vmin, vmax = min(prop_vals), max(prop_vals)
    norm = Normalize(vmin, vmax)
    ticks = (vmin, 0, vmax)

    image = rdmol_prop_heatmap(rdmol, prop=prop, cmap=cmap, norm=norm, **heatmap_args)
    image = imageutils.crop_borders(image, bg_color=WHITE)

    return plotutils.plot_image_with_colorbar(image, cmap, norm, label=cbar_label, ticks=ticks, orient=orient)
