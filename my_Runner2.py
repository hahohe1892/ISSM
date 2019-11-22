import numpy as np
from model import *
from bamg import *
from transient import *
from parameterize import *
from asymetric_geometries2 import Geometry
from solve import *
from InterpFromMeshToMesh2d import *
from setflowequation import *
from setmask import *
from SetMarineIceSheetBC import *
from generate_exp2 import *
from os import path
from cuffey import *
from transientrestart2 import *
from loadmodel import *
from calvingvonmises import *


def SpinUp(params, run_name, load_name):
    y_dim, x_dim, slope, dc, gap_halfwidth, step = standardvalues()
    start_icefront=params['start_icefront']
    slab_thickness=params['slab_thickness']
    steepness=params['steepness']
    null_level=params['null_level']
    x_dim=params['x_dim']


    generateEXP(params, x_dim, y_dim) 
    md=bamg(model(), 'domain','./exp_file2.exp', 'hmax',100, 'hmin', 100) 
    #md=bamg(model(), 'domain', './Square2.exp', 'hmax', 100, 'hmin', 100) 
    md.geometry.bed=Geometry(md, y_dim+20000, x_dim, slope, params['bump_spread'], params['bump_height'], params['bump_pos'], steepness, gap_halfwidth, dc, params['bay_spread1'], params['bay_height1'], params['bay_pos1'], params['bay_spread2'], params['bay_height2'], params['bay_pos2'], step)
    old_mesh_elements=md.mesh.elements
    old_mesh_x=md.mesh.x
    old_mesh_y=md.mesh.y
    old_mesh_geometry=md.geometry.bed
    
    h=np.nan*np.ones(md.mesh.numberofvertices)
    h[np.where(np.logical_and(np.logical_and(md.mesh.y<17800, md.mesh.y>12200), np.logical_and(md.mesh.x<70000, md.mesh.x>18000)))]=100 #155000,120000

    md=bamg(md, 'field', old_mesh_geometry, 'hmax', 2000, 'hmin', params['hmin'], 'gradation', 2, 'hVertices', h)
    md.miscellaneous.name=run_name

    md.geometry.bed=InterpFromMeshToMesh2d(old_mesh_elements, old_mesh_x, old_mesh_y, old_mesh_geometry, md.mesh.x, md.mesh.y)[0][:,0]+null_level  
    md.timestepping.final_time=params['final_time']
    md.timestepping.time_step=params['timestepping']
    md.settings.output_frequency=params['output_frequency']

    md.geometry.surface=md.geometry.bed+1
    surface_pos=np.where(md.mesh.x<=start_icefront)
    #md.geometry.surface[surface_pos]=slab_thickness+null_level+md.mesh.x[surface_pos]*slope
    md.geometry.surface[surface_pos]=np.sqrt(35*(start_icefront-md.mesh.x[surface_pos]))-null_level+150+md.mesh.x[surface_pos]*slope

    below_pos=np.where(md.geometry.surface<md.geometry.bed)
    md.geometry.surface[below_pos]=md.geometry.bed[below_pos]+1
    
    md.geometry.thickness=md.geometry.surface-md.geometry.bed
    md.geometry.base=md.geometry.surface-md.geometry.thickness

    md=setmask(md, '','')
    mask_pos=np.where(md.geometry.surface-md.geometry.bed>params['min_thickness_mask'])
    md.mask.ice_levelset=np.ones(md.mesh.numberofvertices)
    md.mask.ice_levelset[mask_pos]=-1

    #flat_bed=md.geometry.bed-(null_level+slope*md.mesh.x)
    #y_dim_friction=scale(-md.mesh.x)*5
    #md.friction.coefficient=flat_bed*0.03+params['friction']+y_dim_friction
    md.friction.coefficient=params['friction']*np.ones(md.mesh.numberofvertices)

    md=parameterize(md, 'my_Par2.py')
    md=setflowequation(md, 'SSA', 'all')

    md=SetMarineIceSheetBC(md)

    md.levelset.spclevelset=np.nan*np.ones(md.mesh.numberofvertices)
    md.levelset.spclevelset[np.where(md.geometry.bed>0)]=-1
    md.mask.ice_levelset[np.where(md.levelset.spclevelset==-1)]=-1
    
    md.stressbalance.spcvx=np.nan*np.ones(md.mesh.numberofvertices)
    md.stressbalance.spcvx[np.where(md.mesh.x==0)]=params['spcvx']
    md.stressbalance.spcvy=np.nan*np.ones(md.mesh.numberofvertices) 
    md.stressbalance.spcvy[md.mesh.vertexonboundary]=0
    md.stressbalance.spcvy[np.where(md.mesh.x==0)]=np.nan
    md.stressbalance.spcvy[np.where(np.logical_and(md.mesh.x<=x_dim, md.mesh.x>=x_dim-100))]=np.nan
    md.initialization.vx=params['spcvx']*np.ones(md.mesh.numberofvertices)
    md.masstransport.min_thickness=1

    md.frontalforcings.meltingrate=params['frontal_melt']*np.ones(md.mesh.numberofvertices)
    md.basalforcings.floatingice_melting_rate=params['floating_melt']*np.ones(md.mesh.numberofvertices)

    #rheology_B=cuffey(md.initialization.temperature)*1.1
    #md.materials.rheology_B=rheology_B*((md.mesh.y/(np.mean(md.mesh.y)))**0.85)
    #upper_half=np.where(md.mesh.y-10000>0.5*y_dim)
    #md.materials.rheology_B[upper_half]=rheology_B[upper_half]*(((-md.mesh.y[upper_half]+10000-params['bay_height1']+max(md.mesh.y))/(np.mean(-md.mesh.y+10000-params['bay_height1']+max(md.mesh.y))))**0.85)

    md.materials.rheology_B=cuffey(md.initialization.temperature)
    
    md.materials.rheology_law='Cuffey'

    md.calving.stress_threshold_groundedice=params['max_stress']
    md.calving.stress_threshold_floatingice=params['max_stress_floating'] 

    thk_dif=(dc+params['influx_height']+params['null_level'])*np.ones(len(md.mesh.x[np.where(md.mesh.x==0)]))-md.geometry.base[np.where(md.mesh.x==0)]
    md.masstransport.spcthickness[np.where(md.mesh.x==0)]=thk_dif
    md.geometry.thickness[np.where(md.mesh.x==0)]=md.masstransport.spcthickness[np.where(md.mesh.x==0)]
    md.geometry.surface[np.where(md.mesh.x==0)]=md.geometry.base[np.where(md.mesh.x==0)]+md.geometry.thickness[np.where(md.mesh.x==0)]

    #md.transient.isgroundingline=0
    #md.transient.ismovingfront=0
    #md.calving.min_thickness=10
    return md



def SpinUp_load(params, run_name, load_name):

    restart_time=100
    y_dim, x_dim, slope, dc, gap_halfwidth, step = standardvalues()
    x_dim=params['x_dim']
    
    md=loadmodel(load_name) 
    md=transientrestart(md, md, restart_time)
    md.miscellaneous.name=run_name
    md.mask.ice_levelset[np.where(md.results.TransientSolution[restart_time].Thickness<=2)[0]]=1
    md.mask.ice_levelset[np.where(md.levelset.spclevelset==-1)]=-1
    md.timestepping.final_time=params['final_time']
    md.timestepping.time_step=params['timestepping']
    md.settings.output_frequency=params['output_frequency']

    md.frontalforcings.meltingrate=params['frontal_melt']*np.ones(md.mesh.numberofvertices)
    md.basalforcings.floatingice_melting_rate=params['floating_melt']*np.ones(md.mesh.numberofvertices)
    md.stressbalance.spcvx[np.where(md.mesh.x==0)]=params['spcvx']
    md.calving.stress_threshold_groundedice=params['max_stress']
    md.calving.stress_threshold_floatingice=params['max_stress_floating'] 
    md.friction.coefficient=params['friction']*np.ones(md.mesh.numberofvertices)
    thk_dif=(dc+params['influx_height']+params['null_level'])*np.ones(len(md.mesh.x[np.where(md.mesh.x==0)]))-md.geometry.base[np.where(md.mesh.x==0)]
    md.masstransport.spcthickness[np.where(md.mesh.x==0)]=thk_dif
    md.geometry.thickness[np.where(md.mesh.x==0)]=md.masstransport.spcthickness[np.where(md.mesh.x==0)]
    md.geometry.surface[np.where(md.mesh.x==0)]=md.geometry.base[np.where(md.mesh.x==0)]+md.geometry.thickness[np.where(md.mesh.x==0)]
    
    return md

def extenddomain(params, run_name, load_name):
    restart_time=150
    y_dim, x_dim, slope, dc, gap_halfwidth, step = standardvalues()
    start_icefront=params['start_icefront']
    slab_thickness=params['slab_thickness']
    steepness=params['steepness']
    null_level=params['null_level']

    md2=loadmodel(load_name)

    x_dim=params['x_dim']
    generateEXP(params, x_dim, y_dim)
    md=bamg(model(), 'domain','exp_file2.exp', 'hmax',100, 'hmin', 100)

    md.geometry.bed=Geometry(md, y_dim+20000, x_dim, slope, params['bump_spread'], params['bump_height'], params['bump_pos'], steepness, gap_halfwidth, dc, params['bay_spread1'], params['bay_height1'], params['bay_pos1'], params['bay_spread2'], params['bay_height2'], params['bay_pos2'], step)
    old_mesh_elements=md.mesh.elements
    old_mesh_x=md.mesh.x
    old_mesh_y=md.mesh.y
    old_mesh_geometry=md.geometry.bed
    h=np.nan*np.ones(md.mesh.numberofvertices)
    h[np.where(np.logical_and(np.logical_and(md.mesh.y<17800, md.mesh.y>12200), np.logical_and(md.mesh.x<100000, md.mesh.x>52000)))]=100

    ## 70000 und 42000 für spcvx 600 und 3 retreived, new: 90000, 48000, new: 80000, 47000 retreived, new: 50000,75000 or 52000 and 100000
    ## 65000 und 40000 für spcvx 500 und 3 retreived, new: 90000, 43000, from TS17, new: 90000, 50000
    ## 80000 und 45000 für spcvx 500 und 4 retreived, new: 90000, 52000, new: 105000, 55000 retreived, new: spinUp_load
    
    md=bamg(md, 'field', old_mesh_geometry, 'hmax', 1000, 'hmin', params['hmin'], 'gradation', 1.7, 'hVertices', h)
    md.miscellaneous.name=run_name

    md.geometry.bed=InterpFromMeshToMesh2d(old_mesh_elements, old_mesh_x, old_mesh_y, old_mesh_geometry, md.mesh.x, md.mesh.y)[0][:,0]+null_level
    
    md.geometry.thickness=InterpFromMeshToMesh2d(md2.mesh.elements, md2.mesh.x, md2.mesh.y, md2.results.TransientSolution[restart_time].Thickness, md.mesh.x, md.mesh.y)[0][:,0]
    md.mask.groundedice_levelset=InterpFromMeshToMesh2d(md2.mesh.elements, md2.mesh.x, md2.mesh.y, md2.results.TransientSolution[restart_time].MaskGroundediceLevelset, md.mesh.x, md.mesh.y)[0][:,0]
    md.initialization.pressure=InterpFromMeshToMesh2d(md2.mesh.elements, md2.mesh.x, md2.mesh.y, md2.results.TransientSolution[restart_time].Pressure, md.mesh.x, md.mesh.y)[0][:,0]
    md.geometry.base=InterpFromMeshToMesh2d(md2.mesh.elements, md2.mesh.x, md2.mesh.y, md2.results.TransientSolution[restart_time].Base, md.mesh.x, md.mesh.y)[0][:,0]

    deep=np.where(md.geometry.base<md.geometry.bed)
    md.geometry.base[deep]=md.geometry.bed[deep]
    
    
    md.initialization.vx=InterpFromMeshToMesh2d(md2.mesh.elements, md2.mesh.x, md2.mesh.y, md2.results.TransientSolution[restart_time].Vx, md.mesh.x, md.mesh.y)[0][:,0]
    md.initialization.vy=InterpFromMeshToMesh2d(md2.mesh.elements, md2.mesh.x, md2.mesh.y, md2.results.TransientSolution[restart_time].Vy, md.mesh.x, md.mesh.y)[0][:,0]


    ## Parameterization
    md.smb.mass_balance=np.zeros(md.mesh.numberofvertices)
    md.calving=calvingvonmises()

    md.transient.isgroundingline=1
    md.transient.isthermal=1
    md.transient.ismovingfront=1

    md.timestepping.start_time=0
    md.initialization.temperature=(273.15-5.)*np.ones((md.mesh.numberofvertices))
    md.materials.rheology_n=3*np.ones(md.mesh.numberofelements)
    md.friction.p=np.ones(md.mesh.numberofelements) 
    md.friction.q=np.ones(md.mesh.numberofelements) 
    md.basalforcings.groundedice_melting_rate=np.zeros(md.mesh.numberofvertices)

    grounded_mask=np.where(md.mask.groundedice_levelset>0)
    md.geometry.base[grounded_mask]=md.geometry.bed[grounded_mask]
    md.geometry.surface=md.geometry.base+md.geometry.thickness

    
    mask_pos=np.where(md.geometry.surface>1)
    md.mask.ice_levelset=np.ones(md.mesh.numberofvertices)
    md.mask.ice_levelset[mask_pos]=-1
    md.mask.ice_levelset[np.where(md.geometry.thickness>1)]=-1

    md.friction.coefficient=params['friction']*np.ones(md.mesh.numberofvertices)
    
    md.timestepping.final_time=params['final_time']
    md.timestepping.time_step=params['timestepping']
    md.settings.output_frequency=params['output_frequency']

    md=setflowequation(md, 'SSA', 'all')

    md=SetMarineIceSheetBC(md)
    
    md.materials.rheology_B=cuffey(md.initialization.temperature)
    md.materials.rheology_law='Cuffey'

    md.levelset.spclevelset=np.nan*np.ones(md.mesh.numberofvertices)
    md.levelset.spclevelset[np.where(md.geometry.bed>0)]=-1
    md.mask.ice_levelset[np.where(md.levelset.spclevelset==-1)]=-1
    
    md.stressbalance.spcvx=np.nan*np.ones(md.mesh.numberofvertices)
    md.stressbalance.spcvx[np.where(md.mesh.x==0)]=params['spcvx']
    md.stressbalance.spcvy=np.nan*np.ones(md.mesh.numberofvertices) 
    md.stressbalance.spcvy[md.mesh.vertexonboundary]=0
    md.stressbalance.spcvy[np.where(md.mesh.x==0)]=np.nan
    md.stressbalance.spcvy[np.where(np.logical_and(md.mesh.x<=x_dim, md.mesh.x>=x_dim-100))]=np.nan

    md.masstransport.min_thickness=1

    md.frontalforcings.meltingrate=params['frontal_melt']*np.ones(md.mesh.numberofvertices)
    md.basalforcings.floatingice_melting_rate=params['floating_melt']*np.ones(md.mesh.numberofvertices)
    md.cluster=md2.cluster
    
    md.calving.stress_threshold_groundedice=params['max_stress']
    md.calving.stress_threshold_floatingice=params['max_stress_floating'] 
    #md.mask.ice_levelset[np.where(md.geometry.thickness<=2)]=1
    #md.mask.ice_levelset[np.where(md.levelset.spclevelset==-1)]=-1
    #md.initialization.vx[np.where(md.geometry.thickness<=2)]=0
    #md.initialization.vy[np.where(md.geometry.thickness<=2)]=0
    
    thk_dif=(dc+params['influx_height']+params['null_level'])*np.ones(len(md.mesh.x[np.where(md.mesh.x==0)]))-md.geometry.base[np.where(md.mesh.x==0)]
    md.masstransport.spcthickness[np.where(md.mesh.x==0)]=thk_dif
    md.geometry.thickness[np.where(md.mesh.x==0)]=md.masstransport.spcthickness[np.where(md.mesh.x==0)]
    md.geometry.surface[np.where(md.mesh.x==0)]=md.geometry.base[np.where(md.mesh.x==0)]+md.geometry.thickness[np.where(md.mesh.x==0)]

    
    return md


def remesh(params, run_name, load_name):
    restart_time=75
    y_dim, x_dim, slope, dc, gap_halfwidth, step = standardvalues()
    start_icefront=params['start_icefront']
    slab_thickness=params['slab_thickness']
    steepness=params['steepness']
    null_level=params['null_level']

    md2=loadmodel(load_name)

    x_dim=200000
    generateEXP(params, x_dim, y_dim)
    md=bamg(model(), 'domain','exp_file2.exp', 'hmax',100, 'hmin', 100)

    md.geometry.bed=Geometry(md, y_dim+20000, x_dim, slope, params['bump_spread'], params['bump_height'], params['bump_pos'], steepness, gap_halfwidth, dc, params['bay_spread1'], params['bay_height1'], params['bay_pos1'], params['bay_spread2'], params['bay_height2'], params['bay_pos2'], step)
    old_mesh_elements=md.mesh.elements
    old_mesh_x=md.mesh.x
    old_mesh_y=md.mesh.y
    old_mesh_geometry=md.geometry.bed
    
    md=bamg(md, 'field', md2.results.TransientSolution[restart_time].Vel, 'hmax', 1000, 'hmin', params['hmin'], 'gradation', 1.7)
    md.miscellaneous.name=run_name

    md.geometry.bed=InterpFromMeshToMesh2d(old_mesh_elements, old_mesh_x, old_mesh_y, old_mesh_geometry, md.mesh.x, md.mesh.y)[0][:,0]+null_level
    
    md.geometry.thickness=InterpFromMeshToMesh2d(md2.mesh.elements, md2.mesh.x, md2.mesh.y, md2.results.TransientSolution[restart_time].Thickness, md.mesh.x, md.mesh.y)[0][:,0]
    #md.mask.ice_levelset=InterpFromMeshToMesh2d(md2.mesh.elements, md2.mesh.x, md2.mesh.y, md2.results.TransientSolution[120].MaskIceLevelset, md.mesh.x, md.mesh.y)[0][:,0]
    md.mask.groundedice_levelset=InterpFromMeshToMesh2d(md2.mesh.elements, md2.mesh.x, md2.mesh.y, md2.results.TransientSolution[restart_time].MaskGroundediceLevelset, md.mesh.x, md.mesh.y)[0][:,0]
    md.initialization.pressure=InterpFromMeshToMesh2d(md2.mesh.elements, md2.mesh.x, md2.mesh.y, md2.results.TransientSolution[restart_time].Pressure, md.mesh.x, md.mesh.y)[0][:,0]
    md.geometry.base=InterpFromMeshToMesh2d(md2.mesh.elements, md2.mesh.x, md2.mesh.y, md2.results.TransientSolution[restart_time].Base, md.mesh.x, md.mesh.y)[0][:,0]

    deep=np.where(md.geometry.base<md.geometry.bed)
    md.geometry.base[deep]=md.geometry.bed[deep]
    
    
    md.initialization.vx=InterpFromMeshToMesh2d(md2.mesh.elements, md2.mesh.x, md2.mesh.y, md2.results.TransientSolution[restart_time].Vx, md.mesh.x, md.mesh.y)[0][:,0]
    md.initialization.vy=InterpFromMeshToMesh2d(md2.mesh.elements, md2.mesh.x, md2.mesh.y, md2.results.TransientSolution[restart_time].Vy, md.mesh.x, md.mesh.y)[0][:,0]


    ## Parameterization
    #pressure_undefined=np.nonzero([not isinstance(x, float) for x in md.initialization.pressure]) 
    #md.initialization.pressure[pressure_undefined]=md.materials.rho_ice*md.constants.g*md.geometry.thickness
    md.smb.mass_balance=np.zeros(md.mesh.numberofvertices)
    md.calving=calvingvonmises()

    md.transient.isgroundingline=1
    md.transient.isthermal=1
    md.transient.ismovingfront=1

    #Vel_undefined=np.nonzero([not isinstance(x, float) for x in md.initialization.vx])
    #md.initialization.vx[Vel_undefined]=0
    #md.initialization.vy[Vel_undefined]=0
    md.timestepping.start_time=0
    md.initialization.temperature=(273.15-5.)*np.ones((md.mesh.numberofvertices))
    md.materials.rheology_n=3*np.ones(md.mesh.numberofelements)
    md.friction.p=np.ones(md.mesh.numberofelements) 
    md.friction.q=np.ones(md.mesh.numberofelements) 
    md.basalforcings.groundedice_melting_rate=np.zeros(md.mesh.numberofvertices)

    #groundedice_undefined=np.nonzero([not isinstance(x, float) for x in md.mask.groundedice_levelset]) 
    #md.mask.groundedice_levelset[groundedice_undefined]=1

    grounded_mask=np.where(md.mask.groundedice_levelset>0)
    md.geometry.base[grounded_mask]=md.geometry.bed[grounded_mask]
    md.geometry.surface=md.geometry.base+md.geometry.thickness

    
    mask_pos=np.where(md.geometry.surface>1)
    md.mask.ice_levelset=np.ones(md.mesh.numberofvertices)
    md.mask.ice_levelset[mask_pos]=-1
    md.mask.ice_levelset[np.where(md.geometry.thickness>1)]=-1
    md.mask.ice_levelset[np.where(md.levelset.spclevelset==-1)]=-1

    #flat_bed=md.geometry.bed-(null_level+slope*md.mesh.x)
    #y_dim_friction=scale(-md.mesh.x)*5
    #md.friction.coefficient=flat_bed*0.03+params['friction']+y_dim_friction
    md.friction.coefficient=params['friction']*np.ones(md.mesh.numberofvertices)

    md.timestepping.final_time=params['final_time']
    md.timestepping.time_step=params['timestepping']
    md.settings.output_frequency=params['output_frequency']

    md=setflowequation(md, 'SSA', 'all')

    md=SetMarineIceSheetBC(md)
    
    md.materials.rheology_B=cuffey(md.initialization.temperature)
    #md.materials.rheology_B=rheology_B*((md.mesh.y/(np.mean(md.mesh.y)))**0.85)
    #upper_half=np.where(md.mesh.y-10000>0.5*y_dim)
    #md.materials.rheology_B[upper_half]=rheology_B[upper_half]*(((-md.mesh.y[upper_half]+10000-params['bay_height1']+max(md.mesh.y))/(np.mean(-md.mesh.y+10000-params['bay_height1']+max(md.mesh.y))))**0.85)


    md.levelset.spclevelset=np.nan*np.ones(md.mesh.numberofvertices)
    md.levelset.spclevelset[np.where(md.geometry.bed>0)]=-1
    md.mask.ice_levelset[np.where(md.levelset.spclevelset==-1)]=-1
    
    md.stressbalance.spcvx=np.nan*np.ones(md.mesh.numberofvertices)
    md.stressbalance.spcvx[np.where(md.mesh.x==0)]=params['spcvx']
    md.stressbalance.spcvy=np.nan*np.ones(md.mesh.numberofvertices) 
    md.stressbalance.spcvy[md.mesh.vertexonboundary]=0
    md.stressbalance.spcvy[np.where(md.mesh.x==0)]=np.nan
    md.stressbalance.spcvy[np.where(np.logical_and(md.mesh.x<=x_dim, md.mesh.x>=x_dim-100))]=np.nan
    md.mask.ice_levelset[np.where(md.levelset.spclevelset==-1)]=-1 

    md.masstransport.min_thickness=1

    md.frontalforcings.meltingrate=params['frontal_melt']*np.ones(md.mesh.numberofvertices)
    md.basalforcings.floatingice_melting_rate=params['floating_melt']*np.ones(md.mesh.numberofvertices)
    md.cluster=md2.cluster
    md.mask.ice_levelset[np.where(md.results.TransientSolution[-1].Thickness<=50)[0]]=1
    md.mask.ice_levelset[np.where(md.levelset.spclevelset==-1)]=-1

    return md
