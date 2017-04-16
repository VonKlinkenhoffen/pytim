# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding: utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
import MDAnalysis as mda
import pytim
from   pytim.datafiles import *
import numpy as np

u          = mda.Universe(WATER_GRO)
oxygens    = u.select_atoms("name OW")
g=oxygens
print g.atoms.indices
print g.atoms.ids
radii      = pytim_data.vdwradii(G43A1_TOP)
print dir(u._topology)
print u._topology.guessed_attributes
print u.atoms.names
print 'names' in dir(u)
print 'radii' not in dir(u.atoms)
print len(u.atoms)
interface  = pytim.ITIM(u,alpha=2.,max_layers=4,molecular=False)#,multiproc=True,radii_dict=radii,cluster_groups=oxygens,cluster_cut=3.5)
print np.unique(g.tempfactors)
print np.unique(u.atoms.tempfactors)
print u.atoms.radii
print oxygens.radii
print "----", len(u.atoms),interface.universe.atoms.n_atoms
tempf = getattr(u.atoms,'tempfactors')
print "->",np.unique(tempf)

layer      = interface.layers[0,0]  # first layer, upper side
print ("Interface computed. Upper layer:\n %s out of %s" % (layer,oxygens))
print interface._MDAversion
print mda.__version__
print len(interface.layers[0,0])
print len(interface.layers[0,1])

interface.writepdb('layers.pdb',centered=False)

