import numpy as np
### md is model that should get results, md2 is model that contains results, can also be the same, n is time step to be set as new initial value
def transientrestart(md, md2, n):
   # md.results.TransientSolution[n].Vel=np.squeeze(md.results.TransientSolution[n].Vel)
    md2.results.TransientSolution[n].Vx=np.squeeze(md2.results.TransientSolution[n].Vx)
    md2.results.TransientSolution[n].Vy=np.squeeze(md2.results.TransientSolution[n].Vy)
    md2.results.TransientSolution[n].Pressure=np.squeeze(md2.results.TransientSolution[n].Pressure)
    #md.results.TransientSolution[n].Temperature=np.squeeze(md.results.TransientSolution[n].Temperature)

    md2.results.TransientSolution[n].MaskIceLevelset=np.squeeze(md2.results.TransientSolution[n].MaskIceLevelset)
    md2.results.TransientSolution[n].MaskGroundediceLevelset=np.squeeze(md2.results.TransientSolution[n].MaskGroundediceLevelset)

    #md.initialization.vel=md.results.TransientSolution[n].Vel
    md.initialization.vx=md2.results.TransientSolution[n].Vx
    md.initialization.vy=md2.results.TransientSolution[n].Vy
    md.initialization.pressure=md2.results.TransientSolution[n].Pressure
    md.mask.ice_levelset=md2.results.TransientSolution[n].MaskIceLevelset
    md.mask.groundedice_levelset=md2.results.TransientSolution[n].MaskGroundediceLevelset

    md2.results.TransientSolution[n].Base=np.squeeze(md2.results.TransientSolution[n].Base)
    md2.results.TransientSolution[n].Thickness=np.squeeze(md2.results.TransientSolution[n].Thickness)

    deep=np.where(md2.results.TransientSolution[n].Base<md2.geometry.bed)
    md.geometry.base=md2.results.TransientSolution[n].Base
    md.geometry.base[deep]=md2.geometry.bed[deep]
    md.geometry.thickness=md2.results.TransientSolution[n].Thickness
    md.geometry.surface=md.geometry.base+md.geometry.thickness

    return md
