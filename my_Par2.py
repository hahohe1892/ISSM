import numpy as np
from paterson import *
from SetMarineIceSheetBC import *
from setmask import *
from calvingvonmises import *
from cuffey import *


## SMB

md.smb.mass_balance=np.zeros(md.mesh.numberofvertices)
md.calving=calvingvonmises()
#md.levelset.spclevelset=np.nan*np.ones(md.mesh.numberofvertices)
#md.frontalforcings.meltingrate=np.zeros(md.mesh.numberofvertices)

## Transient

md.transient.isgroundingline=1
md.transient.isthermal=1
md.transient.ismovingfront=1

## Constants
md.initialization.pressure=md.materials.rho_ice*md.constants.g*md.geometry.thickness

md.initialization.vx=200*np.ones(md.mesh.numberofvertices)
md.initialization.vy=np.ones(md.mesh.numberofvertices)

## time stepping

md.timestepping.start_time=0

## Materials
md.initialization.temperature=(273.15-5.)*np.ones((md.mesh.numberofvertices))
#md.materials.rheology_B=cuffey(md.initialization.temperature)
md.materials.rheology_n=3*np.ones(md.mesh.numberofelements)

## Boundary conditions
#md=SetMarineIceSheetBC(md) 
#pos=np.where(md.mesh.x==0)
#md.stressbalance.spcvx=np.nan*np.ones(md.mesh.numberofvertices)
#md.stressbalance.spcvy=np.nan*np.ones(md.mesh.numberofvertices)                      
#md.stressbalance.spcvx[pos]=400
#md.stressbalance.spcvy[pos]=0



## Friction

#md.friction.coefficient=90*np.ones(md.mesh.numberofvertices)
#md.friction.coefficient=(1./md.geometry.thickness)*300+90
md.friction.p=np.ones(md.mesh.numberofelements) 
md.friction.q=np.ones(md.mesh.numberofelements) 


## Basal melting

md.basalforcings.groundedice_melting_rate=np.zeros(md.mesh.numberofvertices)
#md.basalforcings.floatingice_melting_rate=np.zeros(md.mesh.numberofvertices)

