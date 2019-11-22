import glob
from matplotlib.backends.backend_pdf import PdfPages
from loadmodel import *
import numpy as np
from plotmodel import *
from GetAreas import *

def generateEXP(params, x_dim, y_dim):
    e=open('./exp_file2.exp', 'w')
    
    x_coordinates=[0, int(params['bay_pos1']-(params['bay_spread1'])/2),int(params['bay_pos1']-(params['bay_spread1'])/2), int(params['bay_pos1']+(params['bay_spread1'])/2), int(params['bay_pos1']+(params['bay_spread1'])/2), x_dim, x_dim, int(params['bay_pos2']+(params['bay_spread2'])/2), int(params['bay_pos2']+(params['bay_spread2'])/2), int(params['bay_pos2']-(params['bay_spread2'])/2), int(params['bay_pos2']-(params['bay_spread2'])/2), 0, 0]
    y_coordinates=[10000, 10000, 10000-params['bay_height1'], 10000-params['bay_height1'], 10000, 10000, 10000+y_dim, 10000+y_dim, 10000+y_dim+params['bay_height2'], 10000+y_dim+params['bay_height2'], 10000+y_dim, 10000+y_dim, 10000]

    p=1
    tuple_list=[]
    for i in range(0, len(x_coordinates)):
        check=(x_coordinates[i], y_coordinates[i])
        if (check in tuple_list)==True:
            continue
        else:
            tuple_list.append(check)
            p+=1

    e.write('## Name:Square\n## Icon:0\n# Points Count  Value\n{} 1.\n# X pos Y pos\n'.format(p))

    tuple_list=[]
    for i in range(0, len(x_coordinates)):
        check=(x_coordinates[i], y_coordinates[i])
        if (check in tuple_list)==True:
            continue
        else:
            tuple_list.append(check)
            e.write('{} {} \n'.format(x_coordinates[i], y_coordinates[i]))
    e.write('0 10000')
    e.close()

    
 
def writelog(params, run_name, x_dim, y_dim, dc, gap_halfwidth, slope, steepness, step, clusterident, output_file, status):
    f=open(output_file, 'a')
    fill_string='%s,'*28+'%s\n'
    f.write(fill_string%(run_name,' ', str(params['spcvx']), str(params['null_level']), str(params['final_time']), str(params['timestepping']), str(params['output_frequency']), str(params['hmin']), str(params['frontal_melt']), str(params['floating_melt']),str(params['friction']),str(params['bump_height']), str(params['bump_pos']), str(params['bump_spread']), str(params['bay_height1']), str(params['bay_pos1']), str(params['bay_spread1']), str(params['bay_height2']), str(params['bay_pos2']), str(params['bay_spread2']), str(x_dim),str(y_dim), str(dc), str(gap_halfwidth), str(slope), str(round(params['steepness'],3)), str(step), clusterident, status))
    f.close()



def standardvalues():
    y_dim=10000
    x_dim=100000
    slope=-2./1000
    dc=2000
    gap_halfwidth=2500
    step=100

    return y_dim, x_dim, slope, dc, gap_halfwidth, step

def plotPDF(pattern, retrieve, name):
    all_model=glob.glob('./Models/'+pattern)
    pattern=pattern[1:].split('_')[0]
    pp=PdfPages('runplots_%s.pdf'%(name))
    for i in range(0,len(all_model)):
        md=loadmodel(all_model[i])
        r=np.linspace(0,len(md.results.TransientSolution),10)
        r = [ int(x) for x in r ]
        on_off={'on': md.miscellaneous.name, 'off': all_model[i]}
        if max(r)<=10:
            pp.savefig(plotmodel(md, 'data', md.results.TransientSolution[0].Vel, 'data', md.results.TransientSolution[-1].Vel, 'colorbar', 'off', 'title',  on_off[retrieve].split('TS')[0]+'\n'+on_off[retrieve].split('TS')[1]+'_Vel'))
            pp.savefig(plotmodel(md, 'data', md.results.TransientSolution[0].Thickness, 'data', md.results.TransientSolution[-1].Thickness, 'colorbar', 'off', 'title', on_off[retrieve].split('TS')[0]+'\n'+on_off[retrieve].split('TS')[1] +'_Thk'))
            pp.savefig(plotmodel(md, 'data', md.results.TransientSolution[0].IceMaskNodeActivation, 'data', md.results.TransientSolution[-1].IceMaskNodeActivation, 'colorbar', 'off', 'title', on_off[retrieve].split('TS')[0]+'\n'+on_off[retrieve].split('TS')[1] +'_IceMask'))
            pp.savefig(plotmodel(md, 'data', md.results.TransientSolution[0].MaskGroundediceLevelset, 'data', md.results.TransientSolution[-1].MaskGroundediceLevelset,'log',10,'log',10,'log',10,'log',10, 'log',10, 'log',10,'log',10,'log',10,'log', 10, 'colorbar', 'off', 'title', on_off[retrieve].split('TS')[0]+'\n'+on_off[retrieve].split('TS')[1] +'_GroundMask'))
        else:
            
            pp.savefig(plotmodel(md, 'data', md.results.TransientSolution[r[0]].Vel,'data', md.results.TransientSolution[r[1]].Vel, 'data',md.results.TransientSolution[r[3]].Vel,'data',md.results.TransientSolution[r[4]].Vel, 'data',md.results.TransientSolution[r[5]].Vel, 'data',md.results.TransientSolution[r[6]].Vel,'data',md.results.TransientSolution[r[7]].Vel, 'data',md.results.TransientSolution[r[8]].Vel, 'data',md.results.TransientSolution[-1].Vel, 'colorbar', 'off', 'title#2',  on_off[retrieve].split('TS')[0]+'\n'+on_off[retrieve].split('TS')[1]+'\n_Vel(%s%s%s%s%s%s%s)'%(r[1], r[3], r[4],r[5], r[6], r[7],r[8])))
            pp.savefig(plotmodel(md, 'data', md.results.TransientSolution[0].Thickness, 'data',md.results.TransientSolution[r[1]].Thickness, 'data',md.results.TransientSolution[r[3]].Thickness,'data',md.results.TransientSolution[r[4]].Thickness, 'data',md.results.TransientSolution[r[5]].Thickness, 'data',md.results.TransientSolution[r[6]].Thickness, 'data',md.results.TransientSolution[r[7]].Thickness, 'data',md.results.TransientSolution[r[8]].Thickness,'data',md.results.TransientSolution[-1].Thickness, 'colorbar', 'off', 'title#2', on_off[retrieve].split('TS')[0]+'\n'+on_off[retrieve].split('TS')[1]+'\n_Thk(%s%s%s%s%s%s%s)'%(r[1], r[3], r[4],r[5], r[6], r[7],r[8])))
            pp.savefig(plotmodel(md, 'data', md.results.TransientSolution[0].IceMaskNodeActivation, 'data',md.results.TransientSolution[r[1]].IceMaskNodeActivation, 'data',md.results.TransientSolution[r[3]].IceMaskNodeActivation, 'data',md.results.TransientSolution[r[4]].IceMaskNodeActivation, 'data',md.results.TransientSolution[r[5]].IceMaskNodeActivation, 'data',md.results.TransientSolution[r[6]].IceMaskNodeActivation, 'data',md.results.TransientSolution[r[7]].IceMaskNodeActivation, 'data',md.results.TransientSolution[r[8]].IceMaskNodeActivation,'data',md.results.TransientSolution[-1].IceMaskNodeActivation, 'colorbar', 'off', 'title#2', on_off[retrieve].split('TS')[0]+'\n'+on_off[retrieve].split('TS')[1]+'\n_IceMask(%s%s%s%s%s%s%s)'%(r[1], r[3], r[4],r[5], r[6], r[7],r[8])))
            pp.savefig(plotmodel(md, 'data', md.results.TransientSolution[0].MaskGroundediceLevelset, 'data',md.results.TransientSolution[r[1]].MaskGroundediceLevelset, 'data',md.results.TransientSolution[r[3]].MaskGroundediceLevelset, 'data',md.results.TransientSolution[r[4]].MaskGroundediceLevelset, 'data',md.results.TransientSolution[r[5]].MaskGroundediceLevelset, 'data',md.results.TransientSolution[r[6]].MaskGroundediceLevelset, 'data',md.results.TransientSolution[r[7]].MaskGroundediceLevelset, 'data',md.results.TransientSolution[r[8]].MaskGroundediceLevelset,'data',md.results.TransientSolution[-2].MaskGroundediceLevelset,'log',10,'log',10,'log',10, 'log',10, 'log',10,'log',10,'log',10,'log', 10,'log',10, 'colorbar', 'off', 'title#2', on_off[retrieve].split('TS')[0]+'\n'+on_off[retrieve].split('TS')[1]+'\n_GroundMask(%s%s%s%s%s%s%s)'%(r[1], r[3], r[4],r[5], r[6], r[7],r[8])))

        fig=plt.figure(10)
        for t in range(0, len(md2.results.TransientSolution)):
            along(md2,md2.results.TransientSolution[t].Surface)
            along(md2, md2.results.TransientSolution[t].Base)
            along(md2, md2.results.TransientSolution[t].Vel)
        pp.savefig(plt.figure(10))
        plt.close(10)

        mod=md3
        ice_Volume=[]
        floating_area=[]
        grounded_area=[]
        calving=[]
        fV=[]
        for q in range(0, (len(mod.results.TransientSolution)-1)):
            ice_Volume.append(mod.results.TransientSolution[q].IceVolume)
            grounded_area.append(mod.results.TransientSolution[q].GroundedArea)
            floating_area.append(mod.results.TransientSolution[q].FloatingArea)
            if hasattr(mod.results.TransientSolution[0], 'TotalCalvingFluxLevelset'):
                calving.append(mod.results.TransientSolution[q].TotalCalvingFluxLevelset)
            fV.append(mod.results.TransientSolution[q].IceVolume-mod.results.TransientSolution[q].IceVolumeAboveFloatation)
                

        mval=[]
        mvel=[]
        mthk=[]
        mvy=[]
        for i in range(0, len(mod.results.TransientSolution)):
            mval.append(np.max(mod.mesh.x[np.where(mod.results.TransientSolution[i].Thickness>30)[0]]))
            mvel.append(np.max(mod.results.TransientSolution[i].Vel))
            mthk.append(np.max(mod.results.TransientSolution[i].Thickness))
            mvy.append(np.max(abs(mod.results.TransientSolution[i].Vy)))
            
        fig=plt.figure(1, dpi=20)
        plt.plot(ice_Volume, label=i)
        fig.suptitle('Ice Volume')
        fig=plt.figure(4, dpi=20)
        plt.plot(ice_Volume, label=i)
        fig.suptitle('Ice Volume')
        fig=plt.figure(2, dpi=20)
        plt.plot(grounded_area, label=i)
        fig.suptitle('Grounded Area')
        fig=plt.figure(3, dpi=20)
        plt.plot(floating_area, label=i)
        fig.suptitle('Floating Area')
        fig=plt.figure(5, dpi=20)
        plt.plot(calving, label=i)
        fig.suptitle('Calving Flux')

        
    plt.figure(1)
    plt.legend()
    plt.figure(2)
    plt.legend()
    plt.figure(3)
    plt.legend()
    plt.figure(4)
    plt.legend()
    plt.figure(5)
    plt.legend()
    pp.savefig(plt.figure(1))
    pp.savefig(plt.figure(2))
    pp.savefig(plt.figure(3))
    pp.savefig(plt.figure(4))
    pp.savefig(plt.figure(5))
    plt.close('all')
    pp.close()




def scale(data):
    output=(data-np.mean(data))/np.std(data)

    return output


def across(md, parameter, limit1, limit2):
    across=np.where(np.logical_and(md.mesh.x<limit2, md.mesh.x>limit1))
    array=np.array((md.mesh.y[across], np.squeeze(parameter[across])))
    ind=np.argsort(array[0])
    array=array[:,ind]
    plt.plot(array[0],array[1])

def along(md, parameter, limit1=15001, limit2=14999):
    along=np.where(np.logical_and(md.mesh.y<limit1, md.mesh.y>limit2))
    array=np.array((md.mesh.x[along], np.squeeze(parameter[along])))
    ind=np.argsort(array[0])
    array=array[:,ind]
    plt.plot(array[0], array[1])



def getcontrib(md):
    areas=GetAreas(md.mesh.elements, md.mesh.x, md.mesh.y)
    
