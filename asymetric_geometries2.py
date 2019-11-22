#    from model import *
#    from bamg import *
#    from plotmodel import *

#    x_dim=100000 #length of setup
#    y_dim=40000 #width of setup
#    x_mesh=1000 #resolution in x direction
#    y_mesh=400 #resolution in y direction
#    x_res=x_dim/x_mesh
#    y_res=y_dim/y_mesh
#
#    md=bamg(model(), 'domain', '../../Exp/Square.exp', 'hmax', 500)

    ## Topography in x-direction

#    slope=-1/1000

#    bump_spread=20000
#    bump_height=500
#    bump_pos=60000

import numpy as np

def Geometry(md, y_dim, x_dim, slope, bump_spread, bump_height, bump_pos, steepness, gap_halfwidth, dc, bay_spread1, bay_height1, bay_pos1, bay_spread2, bay_height2, bay_pos2, step):

    #slope_area1=np.where(np.logical_and(md.mesh.x>=0,md.mesh.x<bump_pos-bump_spread/2))
    #slope_area2=np.where(md.mesh.x>bump_pos+bump_spread/2)
    bump_area=np.where(np.logical_and(md.mesh.x>bump_pos-bump_spread/2, md.mesh.x<bump_pos+bump_spread/2))

    Bx=np.ones(md.mesh.numberofvertices)*md.mesh.x*slope
    Bx[bump_area]=-np.sin((md.mesh.x[bump_area]-bump_spread/4-bump_pos)*2*np.pi/bump_spread)*.5*bump_height+.5*bump_height+md.mesh.x[bump_area]*slope



    #steepness=1/200 # steepness of the fjord walls
    #gap_halfwidth=5000 # this times 2 yields width at half of fjord depth, for asymetric fjords create gap1 and gap2 and put that into equation
    #dc=800 #m, fjord depth

    #bay_spread1=50000 #extent of bay/bottleneck
    #bay_height1=7000 #value that is added/subtracted from gap_halfwidth to create embayment/bottleneck 
    #bay_pos1=60000

    #bay_spread2=50000
    #bay_height2=7000
    #bay_pos2=60000




    By=np.ones(md.mesh.numberofvertices)*(dc/(1+np.exp(steepness*(md.mesh.y-y_dim/2+gap_halfwidth)))+dc/(1+np.exp(-steepness*(md.mesh.y-y_dim/2-gap_halfwidth))))
    bay_area1=np.where(np.logical_and(md.mesh.x>bay_pos1-bay_spread1/2, md.mesh.x<bay_pos1+bay_spread1/2))
    bay_area2=np.where(np.logical_and(md.mesh.x>bay_pos2-bay_spread2/2, md.mesh.x<bay_pos2+bay_spread2/2))
    bay1=-np.sin((md.mesh.x[bay_area1]-bay_spread1/4-bay_pos1)*2*np.pi/bay_spread1)*.5*bay_height1+.5*bay_height1
    bay2=-np.sin((md.mesh.x[bay_area2]-bay_spread2/4-bay_pos2)*2*np.pi/bay_spread2)*.5*bay_height2+.5*bay_height2

    overlap=np.intersect1d(bay_area1, bay_area2)

    if len(overlap)==0:
        By[bay_area1]=dc/(1+np.exp(steepness*(md.mesh.y[bay_area1]-y_dim/2+gap_halfwidth+bay1)))+dc/(1+np.exp(-steepness*(md.mesh.y[bay_area1]-y_dim/2-gap_halfwidth)))
        By[bay_area2]=dc/(1+np.exp(steepness*(md.mesh.y[bay_area2]-y_dim/2+gap_halfwidth)))+dc/(1+np.exp(-steepness*(md.mesh.y[bay_area2]-y_dim/2-gap_halfwidth-bay2)))

    else:
        By[bay_area1]=dc/(1+np.exp(steepness*(md.mesh.y[bay_area1]-y_dim/2+gap_halfwidth+bay1)))
        By[bay_area2]=By[bay_area2]+dc/(1+np.exp(-steepness*(md.mesh.y[bay_area2]-y_dim/2-gap_halfwidth-bay2)))


    #step=100
    #canyon=np.where(By<1)[0]
    #test=np.zeros(md.mesh.numberofvertices)
    iterator=np.arange(step,x_dim+step,step)
    only_bump=Bx+-slope*md.mesh.x

    #for r in np.arange(0,len(canyon)):
       # x_canyon=md.mesh.x[canyon[r]]
       # buff=[x_canyon-100, x_canyon+100]
      #  vertices_in_buffer=np.where(np.logical_and(md.mesh.x>buff[0],md.mesh.x<buff[1]))
     #   By[np.where(only_bump[vertices_in_buffer]>By[canyon[r]])]=only_bump[np.where(only_bump[vertices_in_buffer]>By[canyon[r]])]
     #   if len(np.where(only_bump[vertices_in_buffer]>By[canyon[r]])[0])>0:
    #           print(r)


    if bump_height>=0:

        for i in iterator:
            uplift_area=np.where(np.logical_and(md.mesh.x<i, md.mesh.x>i-step))[0]
    #    canyon_reference=np.intersect1d(uplift_area, canyon)[0]
            too_low=np.where(only_bump[uplift_area]>By[uplift_area])

            By[uplift_area[too_low]]=only_bump[uplift_area[too_low]]
            md.geometry.bed=By+md.mesh.x*slope


    else:
        for i in iterator:
            uplift_area=np.where(np.logical_and(md.mesh.x<i, md.mesh.x>i-step))[0]
    #    canyon_reference=np.intersect1d(uplift_area, canyon)[0]
            too_high=np.where(-only_bump[uplift_area]>By[uplift_area])

            By[uplift_area[too_high]]=only_bump[uplift_area[too_high]]
            md.geometry.bed=By+md.mesh.x*slope


#        without_bump=By+md.mesh.x*slope
#        plateau=np.where(without_bump>dc/2-x_dim*slope)
#        md.geometry.bed=Bx+By
#        md.geometry.bed[plateau]=without_bump[plateau]




    return md.geometry.bed

#plotmodel(md, 'data', md.geometry.bed)






























#    CreBump1=CreateBump(bump_spread, bump_height, bump_pos, md, x_dim, y_dim, x_mesh, y_mesh)
#    CreBump2=CreateBump(bump_spread, bump_height, bump_pos, md, x_dim, y_dim, x_mesh, y_mesh)

#    bump1=CreBump1[0]
#   insert_left1=CreBump1[1]
#    insert_right1=CreBump1[2]
#    insert_bump1=CreBump1[3]


#    Bx[insert_left1:insert_right1]=bump1[0:insert_bump1]+Bx[insert_left1:insert_right1]


#    steepness=1/200 # steepness of the fjord walls
#    gap_halfwidth=5000 # this times 2 yields width at half of fjord depth, for asymetric fjords create gap1 and gap2 and put that into equation
#    dc=800 #m, fjord depth

#    bay_spread1=50000 #extent of bay/bottleneck
#    bay_height1=-7000 #value that is added/subtracted from gap_halfwidth to create embayment/bottleneck 
#    bay_pos1=60000

#    bay_spread2=50000
#    bay_height2=7000
#    bay_pos2=60000

#    CreBay1=CreateBump(bay_spread1,bay_height1,bay_pos1,md,x_dim, y_dim,x_mesh,y_mesh)

#    bay1=CreBay1[0]
#    insert_left_bay1=CreBay1[1]
#    insert_right_bay1=CreBay1[2]
#    insert_bay1=CreBay1[3]

#    CreBay2=CreateBump(bay_spread2, bay_height2, bay_pos2, md, x_dim,y_dim,x_mesh,y_mesh)

#    bay2=CreBay2[0]
#    insert_left_bay2=CreBay2[1]
#    insert_right_bay2=CreBay2[2]
#    insert_bay2=CreBay2[3]

#    By=dc/(1+np.exp(steepness*(md.mesh.y-y_dim/2+gap_halfwidth)))+dc/(1+np.exp(-steepness*(md.mesh.y-y_dim/2-gap_halfwidth)))


#    set1=set(np.array(range(insert_left_bay1,insert_right_bay1)))
#    set2=set(np.array(range(insert_left_bay2,insert_right_bay2)))
#    overlap=list(set1.intersection(set2))

#    if len(overlap)==0:  ## partial overlap not supported, either complete overlap or no overlap

#        By[insert_left_bay1:insert_right_bay1]=dc/(1+np.exp(steepness*(md.mesh.y[insert_left_bay1:insert_right_bay1]-y_dim/2+(gap_halfwidth+bay1[0:insert_bay1]))))+dc/(1+np.exp(-steepness*(md.mesh.y[insert_left_bay1:insert_right_bay1]-y_dim/2-gap_halfwidth)))

#        By[insert_left_bay2:insert_right_bay2]=dc/(1+np.exp(steepness*(md.mesh.y[insert_left_bay2:insert_right_bay2]-y_dim/2+(gap_halfwidth))))+dc/(1+np.exp(-steepness*(md.mesh.y[insert_left_bay2:insert_right_bay2]-y_dim/2-gap_halfwidth-bay2[0:insert_bay2])))

#    else:

#        By[insert_left_bay1:insert_right_bay1]=dc/(1+np.exp(steepness*(md.mesh.y[insert_left_bay1:insert_right_bay1]-y_dim/2+(gap_halfwidth+bay1[0:insert_bay1]))))

#        By[insert_left_bay2:insert_right_bay2]=dc/(1+np.exp(-steepness*(md.mesh.y[insert_left_bay2:insert_right_bay2]-y_dim/2-gap_halfwidth-bay2[0:insert_bay2])))+By[insert_left_bay2:insert_right_bay2]



#    if bump_height>=0:
    ## for each 'line' in y direction, select those vertices that are lower than the centerline with bump and for each such vertice, add the corresponding centerline value to a list
#        without_bump=By+md.mesh.x*slope
#        new_elevation=[]
#        adjust_locations=[]
#        extract_lowpoint=[]
#        for t in np.array(range(1,x_mesh+1)):
#            row_range=np.array(range(t*y_mesh-y_mesh, t*y_mesh))
#            extract_lowpoint.append(np.where(without_bump[row_range]==min(without_bump[row_range]))[0][0]+(t-1)*y_mesh)
#        for t in np.array(range(1,x_mesh+1)):
#            bump_c=Bx[extract_lowpoint]
#            for q in np.array(range(t*y_mesh-y_mesh, t*y_mesh)):
#               q_mesh=(slope)*md.mesh.x[q]+By[q]
#                if q_mesh <bump_c[t-1]:
#                    adjust_locations.append(q)
#                    new_elevation.append(bump_c[t-1])
#        md.geometry.bed=without_bump  ## computate the bed without bump
#        for r in range(1,len(adjust_locations)-1):  ## replace old values with new elevation of bump
#            md.geometry.bed[adjust_locations[r]]=new_elevation[r]
#
#    else:
#    ## calculate geometry with and without depression; all values higher than half the fjord depth are set to their value without depression, all values below to their value with depression
#        without_bump=By+md.mesh.x*slope
#        plateau=np.where(without_bump>dc/2-x_dim*slope)
#        md.geometry.bed=Bx+By
#        md.geometry.bed[plateau]=without_bump[plateau]

#    plotmodel(md, 'data', md.geometry.bed)
