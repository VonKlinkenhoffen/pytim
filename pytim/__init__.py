# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding: utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4

from   abc import ABCMeta, abstractmethod, abstractproperty
from   distutils.version import LooseVersion
import numpy as np
import MDAnalysis
from   MDAnalysis.topology import tables
from   difflib import get_close_matches



def PatchTrajectory(trajectory,interface):

    trajectory.interface = interface
    trajectory.original_read_next_timestep = trajectory._read_next_timestep

    class PatchedTrajectory(trajectory.__class__):

        def _read_next_timestep(self,ts=None):
            tmp = self.original_read_next_timestep(ts=ts)
            self.interface._assign_layers()
            return tmp

    oldname  = trajectory.__class__.__name__
    oldmodule= trajectory.__class__.__module__

    PatchedTrajectory.__name__   = oldname
    PatchedTrajectory.__module__ = oldmodule
    trajectory.__class__ = PatchedTrajectory


class PYTIM(object):
    """ The PYTIM metaclass
    """
    __metaclass__ = ABCMeta

    directions_dict={0:'x',1:'y',2:'z','x':'x','y':'y','z':'z','X':'x','Y':'y','Z:':'z'}
    symmetry_dict={'cylindrical':'cylindrical','spherical':'spherical','planar':'planar'}

    ALPHA_NEGATIVE = "parameter alpha must be positive"
    ALPHA_LARGE= "parameter alpha must be smaller than the smaller box side"
    MESH_NEGATIVE = "parameter mesh must be positive"
    MESH_LARGE= "parameter mesh must be smaller than the smaller box side"
    UNDEFINED_RADIUS= "one or more atoms do not have a corresponding radius in the default or provided dictionary"
    UNDEFINED_CLUSTER_SEARCH= "If extra_cluster_groups is defined, a cluster_cut should e provided"
    MISMATCH_CLUSTER_SEARCH= "cluster_cut should be either a scalar or an array matching the number of groups (including itim_group)"
    EMPTY_LAYER="One or more layers are empty"
    CLUSTER_FAILURE="Cluster algorithm failed: too small cluster cutoff provided?"
    UNDEFINED_LAYER="No layer defined: forgot to call _assign_layers() or not enough layers requested"
    WRONG_UNIVERSE="Wrong Universe passed to ITIM class"
    UNDEFINED_ITIM_GROUP="No itim_group defined"
    WRONG_DIRECTION="Wrong direction supplied. Use 'x','y','z' , 'X', 'Y', 'Z' or 0, 1, 2"
    CENTERING_FAILURE="Cannot center the group in the box. Wrong direction supplied?"


    def _basic_checks(self,universe):
        self._MDAversion=MDAnalysis.__version__
        assert LooseVersion(self._MDAversion) >  LooseVersion('0.15'), "Must use MDAnalysis  >= 0.15"

        assert isinstance(universe,MDAnalysis.core.universe.Universe), "You must pass an MDAnalysis Universe"

        if LooseVersion(self._MDAversion) >= LooseVersion('0.16'):  # new topology system
            if 'radii' not in dir(universe.atoms):
                from MDAnalysis.core.topologyattrs import Radii
                radii=np.zeros(len(universe.atoms))
                universe.add_TopologyAttr(Radii(radii))
            else:
                assert isinstance(universe.atoms.radii, property), "Internal error, 'radii' is not a property"
                
                
        

    def writepdb(self,filename='layers.pdb',centered='no',multiframe=True):
        """ Write the frame to a pdb file, marking the atoms belonging
            to the layers with different beta factor.

            :param filename:   string  -- the output file name
            :param centered:   string  -- 'origin', 'middle', or 'no'
            :param multiframe: boolean -- append to pdb file if True

            Example: save the positions (centering the interface in the cell) without appending

            >>> interface.writepdb('layers.pdb',multiframe=False)

            Example: save the positions without centering the interface. This will
                     not shift the atoms from the original position (still, they will
                     be put into the basic cell).
                     The :multiframe: option set to :False: will overwrite the file.

            >>> interface.writepdb('layers.pdb',centered='no')

        """
        center_options=['no','middle','origin',False,True]
        if centered not in center_options:
            centered='no'
        if centered == False: centered='no'
        if centered == True : centered='middle'
        try:
            if centered=='no':
                self.universe.atoms.positions=self.original_positions

            if centered=='middle':
                # NOTE: this assumes that all method relying on 'planar' symmetry must center the interface along the normal
                if self.symmetry=='planar':
                    translation = [0,0,0]
                    translation[self.normal] = self.universe.dimensions[self.normal]/2.
                    self.universe.atoms.positions+=np.array(translation)
                    self.universe.atoms.pack_into_box(self.universe.dimensions[:3])
            try:
                id(self.PDB[filename])>0 # it exists already, let's add information about the box, as MDAnalysis forgets to do so for successive frames. TODO MDA version check here!
                self.PDB[filename].CRYST1(self.PDB[filename].convert_dimensions_to_unitcell(self.universe.trajectory.ts))
            except:
                self.PDB[filename]=MDAnalysis.Writer(filename, multiframe=True, bonds=False,
                            n_atoms=self.universe.atoms.n_atoms)

            self.PDB[filename].write(self.universe.atoms)

        except:
            print("Error writing pdb file")

    def savepdb(self,filename='layers.pdb',centered='no',multiframe=True):
        """ An alias to :func:`writepdb`
        """
        self.writepdb(filename,centered,multiframe)

    def assign_radii(self,radii_dict):
        try:
            _groups = np.copy(self.extra_cluster_groups[:])
        except:
            _groups = []
        _groups.append(self.itim_group)
        for _g in _groups:
            # TODO: add a switch that allows to use the atom name instead of the type!
            if _g is not None:
                _types = np.copy(_g.types)
                if not np.any(np.equal(_g.radii, None)): # all radii already set
                    break
                if radii_dict is None : # some radii are not set  and no dict provided
                    _radii_dict = tables.vdwradii
                else: # use the provided dict.
                    _radii_dict = radii_dict

                _radii = np.zeros(len(_g.types))
                for _atype in np.unique(_types):
                    try:
                        matching_type         = get_close_matches(_atype, _radii_dict.keys() , n=1, cutoff=0.1)
                        _radii[_types==_atype]=_radii_dict[matching_type[0]]
                    except:
                        avg =  np.average(_radii_dict.values())
                        _radii[_types==_atype] = avg

                        print("!!                              WARNING")
                        print("!! No appropriate radius was found for the atomtype "+_atype)
                        print("!! Using the average radius ("+str(avg)+") as a fallback option...")
                        print("!! Pass a dictionary of radii (in Angstrom) with the option radii_dict")
                        print("!! for example: r={'"+_atype+"':1.2,...} ; inter=pytim.ITIM(u,radii_dict=r)")

                _g.radii=np.copy(_radii[:])

                assert not np.any(np.equal(_g.radii,None)) , self.UNDEFINED_RADIUS
                del _radii
                del _types

    def _assign_normal(self,normal):
        assert self.symmetry=='planar', "Error: wrong symmetry for normal assignement"
        assert self.itim_group is not None, self.UNDEFINED_ITIM_GROUP
        if normal=='guess':
            self.normal=utilities.guess_normal(self.universe,self.itim_group)
        else:
            assert normal in self.directions_dict, self.WRONG_DIRECTION
            self.normal = self.directions_dict[normal]


    def center(self, group, direction=None, halfbox_shift=True):
        """
        Centers the liquid slab in the simulation box.

        The algorithm tries to avoid problems with the definition
        of the center of mass. First, a rough density profile
        (10 bins) is computed. Then, the group is shifted
        and reboxed until the bins at the box boundaries have a
        density lower than a threshold delta


        In ITIM, the system along the normal direction is always
        centered at 0 (halfbox_shift==True). To center to the middle
        of the box along all directions, set halfbox_shift=False

        """
        if direction is None:
            direction = self.normal
        dim = group.universe.coord.dimensions
        total_shift=0

        assert direction in self.directions_dict, self.WRONG_DIRECTION
        _dir = self.directions_dict[direction]
        if _dir == 'x':
            direction = 0
            _pos_group =  utilities.get_x(group)
        if _dir == 'y':
            direction = 1
            _pos_group =  utilities.get_y(group)
        if _dir == 'z':
            direction = 2
            _pos_group =  utilities.get_z(group)

        shift=dim[direction]/100. ;


        _x = utilities.get_x(group.universe.atoms)
        _y = utilities.get_y(group.universe.atoms)
        _z = utilities.get_z(group.universe.atoms)

        if(halfbox_shift==True):
            _range=(-dim[direction]/2.,dim[direction]/2.)
        else:
            _range=(0.,dim[direction])

        histo,edges=np.histogram(_pos_group, bins=10, range=_range, density=True)

        max_val=np.amax(histo)
        min_val=np.amin(histo)
        delta=min_val+(max_val-min_val)/3. ;# TODO test different cases

        # let's first avoid crossing pbc with the liquid phase. This can fail:
        while(histo[0]>delta or histo[-1]> delta):
            total_shift+=shift
            assert total_shift<dim[direction], self.CENTERING_FAILURE
            _pos_group +=shift
            if _dir == 'x':
                utilities.centerbox(group.universe,x=_pos_group,center_direction=direction,halfbox_shift=halfbox_shift)
            if _dir == 'y':
                utilities.centerbox(group.universe,y=_pos_group,center_direction=direction,halfbox_shift=halfbox_shift)
            if _dir == 'z':
                utilities.centerbox(group.universe,z=_pos_group,center_direction=direction,halfbox_shift=halfbox_shift)
            histo,edges=np.histogram(_pos_group, bins=10, range=_range, density=True)

        _center=np.average(_pos_group)
        if(halfbox_shift==False):
            box_half = dim[direction]/2.
        else:
            box_half = 0.
        if _dir == 'x':
            _x += total_shift - _center + box_half
        if _dir == 'y':
            _y += total_shift - _center + box_half
        if _dir == 'z':
            _z += total_shift - _center + box_half
        # finally, we copy everything back
        group.universe.atoms.positions=np.column_stack((_x,_y,_z))

    @abstractmethod
    def _assign_layers(self):
        pass

    @abstractmethod
    def triangulate_layer(self):
        pass

    @abstractmethod
    def interpolate_surface(self):
        pass

    @abstractproperty
    def layers(self):
        pass

    def _define_groups(self):
        # we first make sure cluster_cut is either None, or an array
        if self.cluster_cut is not None and not isinstance(self.cluster_cut, (list, tuple, np.ndarray)):
            if type(self.cluster_cut) is int or type(self.cluster_cut) is float:
                self.cluster_cut = np.array([float(self.cluster_cut)])
        # same with extra_cluster_groups
        if self.extra_cluster_groups is not None and not isinstance(self.extra_cluster_groups, (list, tuple, np.ndarray)):
            self.extra_cluster_groups = [self.extra_cluster_groups]

        # fallback for itim_group
        if self.itim_group is None:
            self.itim_group = self.all_atoms


from pytim.itim import  ITIM
from pytim.gitim import GITIM
from pytim import observables, utilities
from pytim.datafiles import *

#__all__ = [ 'itim' , 'gitim' , 'observables', 'datafiles', 'utilities']


