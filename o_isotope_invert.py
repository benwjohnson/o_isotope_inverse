#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 18 14:22:01 2018

@author: benjohnson
"""

##----Function that takes column of oxygen isotope, temperature, and spatial data
# (need horizontal (x) and vertical parts (z) of a cross section) and outputs fluid vector, 
# incoming and outgoing fluid composition

def o_isotope_invert(del18O,temp,x,z,d18Oinit):
    import sys
    sys.path.insert(0, '/Users/benjohnson/Science/python_scripts/')
    sys.path.insert(0, '/Users/benjohnson/Science/projects/geochemical_inverse_modeling')
    import numpy as np
    from scipy.interpolate import griddata
    import matplotlib.pyplot as plt
    from matrix_average import matrix_average 
    #basic stuff
    #some constants
    Rstd=0.00205 #ratio of 18O/16O in SMOW
    Orock=8.6e6 #moles O per m^3 rock
#    Orock=4*Orock #moles O in plag per m^3 in rock
    d18Of = np.nanmean(del18O)
    rho_rock=3 #density in g/cm^3

    ##inversion
    numsteps=10

    Zmax = max(z); Zmin=min(z)
    Xmax = max(x); Xmin=min(x)

    xi = np.arange(Xmin,Xmax+1,(Xmax-Xmin)/numsteps)
    zi = np.arange(Zmin,Zmax+1,(Zmax-Zmin)/numsteps)
    xi,zi = np.meshgrid(xi,zi)  

    #fill value sets values outside grid area to 5.5 permil, or average basalt d18O initial    
    iso_grid = griddata((x,z),del18O,(xi,zi),method='linear')
    
    #find rows with only nans
    num_rows = iso_grid.shape[0]
    rows2del = []
    for irow in range(0,num_rows):
        rownan = np.argwhere(np.isnan(iso_grid[irow]))
        if len(rownan) == (iso_grid.shape[1]):
            rows2del.append(irow)
    iso_grid = np.delete(iso_grid,rows2del,0)     
    from nan_helper import nan_helper
   
    num_rows = iso_grid.shape[0]
    for irow in range(0,num_rows):
        temp_row = iso_grid[irow]
        nans, tempx = nan_helper(iso_grid[irow])
        temp_row[nans] = np.interp(tempx(nans), tempx(~nans), temp_row[~nans])
        iso_grid[irow] = temp_row
     
    iso_grid_trim = iso_grid #uncomment if trimming nans
    x_grid_trim = np.delete(xi,rows2del,0)
    z_grid_trim = np.delete(zi,rows2del,0)
    
        
    depth = np.subtract(max(z),z)  
    ktemp = 20000
    #temp = np.add(5,np.multiply(350,np.divide(depth,np.add(depth,ktemp))))
   # temp = np.multiply(np.divide(1,del18O),1000)
    geothermx,geothermy = np.meshgrid(depth,temp) 
    #same steps for temperature
    temp_grid = griddata((x,z),temp,(xi,zi),method='linear') #comment if trimming nans below
    temp_grid = np.delete(temp_grid,rows2del,0)    
    
    num_rows = temp_grid.shape[0]
    for irow in range(0,num_rows):
        temp_row = temp_grid[irow]
        nans, tempx = nan_helper(temp_grid[irow])
        temp_row[nans] = np.interp(tempx(nans), tempx(~nans), temp_row[~nans])
        temp_grid[irow] = temp_row
    
   
    temp_grid_trim = temp_grid
#    
#    
#    fig73=plt.figure(num=73); plt.clf(); #plt.title('grid')
#    plt.contourf(x_grid_trim,z_grid_trim,temp_grid_trim,cmap=plt.cm.BrBG)
#    cm = plt.cm.get_cmap('RdYlBu')
#    cbar=plt.colorbar()
#    cbar.ax.get_yaxis().labelpad = 15
#    cbar.ax.set_ylabel('Temp (degC)',rotation=270)
#    plt.xlim([45000,375000]); plt.ylim([-1000,110000])
#    
#    plt.plot(x,z,'w^',markeredgecolor='k',markersize=12)
#     
    #
    #find first column with no nans
#    num_cols = temp_grid.shape[1]
#    cols2del=[]
#    for first in range(num_cols):
#        [idx] = np.where(np.isfinite(temp_grid[first]))
#        cols2del.append(first)
#        if idx.size>0:
#            break
        
    #loop to interpolate out to edges of grid to make it totally square...
#    from nan_helper import nan_helper
#    
#    for icol in range(first,num_cols):
#        temp_col = temp_grid[:,icol]
#        nans, tempx = nan_helper(temp_col)
#        temp_col[nans] = np.interp(tempx(nans), tempx(~nans), temp_col[~nans])
#        temp_grid[:,icol] = temp_col
    #    
    ##delete all columns with all nans
    #temp_grid_trim = np.delete(temp_grid,cols2del,1)
    #x_grid_trim_temp=np.delete(xi,cols2del,1)
    #z_grid_trim_temp=np.delete(zi,cols2del,1)
    
    #### Now, processing and inverting
    #say grid is regularly spaced, with 100m sides, and is 1cm thick, use units of m
    
    length=x_grid_trim[0,1]-x_grid_trim[0,0]; #all in m
    height=z_grid_trim[1,0]-z_grid_trim[0,0];
    width=0.01
    volume=length*height*width
    
    moles_O_box = Orock*volume
    moles18O_init_box = (5.5/1000+1)*Rstd*Orock*volume
    
    [auto_test,col_avgs_test,row_avgs_test] = matrix_average(iso_grid_trim)
    
    grid_dim = x_grid_trim.shape
    num_rows = grid_dim[0]
    num_cols = grid_dim[1]
    total_raw_entries = auto_test.size
    total_avg_entries = (num_cols-1)*num_rows + (num_rows-1)*num_cols# total number of averages
    
    #use matrix_average to get temp 
    [temp_test,temp_cols,temp_rows] = matrix_average(temp_grid_trim)
    
    temp_test_kelvin = temp_test[0]+273
    
    ##calculate water del18O in equilibrium with rock at temps
    epsilon_list_test=[]
    epsilon_r_w = []
    AA = -3.9e3; BB = 3.42e6  #1000ln alpha(rock-water) = AA/T + BB/T^2  Cole et al., 1987
#    AA = 2.15e6 ; BB = 0.692; #1000ln alpha(r-w) = AAx10^6/T^2+BB Cole, 1980 and Cathles, 1983
#    AA = 2.53e6 ; BB = 3.61; #1000ln alpha(r-w) = AAx10^6/T^2+BB Norton and Taylor, 1979
    for ientry in range(len(temp_test_kelvin)):
        aa = AA/temp_test_kelvin[ientry]; bb = BB/(temp_test_kelvin[ientry]**2) #1000ln alpha(rock-water) = AA/T + BB/T^2 
#        aa = AA/(temp_test_kelvin[ientry]**2); bb=BB  #1000ln alpha(r-w) = AAx10^6/T^2+BB
        alpha = np.exp((aa+bb)*1e-3)
        alphaw = 1/alpha
        epsilon_list_test.append((alphaw-1)*1000)
        epsilon_r_w.append((alpha-1)*1000)
        
    water_isotope_test = np.add(auto_test,epsilon_list_test)
    #
    ##setting up Cq=delC to solve for q, fluid flow vector
    ##Basic properties
    num_boxes =  (num_rows-1)*(num_cols-1)   
    
    
    #C
    conc_matrix_test = np.zeros([num_boxes*2,total_avg_entries])
    idx=[]  
    
    #find indices for each box, loop puts out list where each row is flow and concentration indices for each box
    start_idx = -(num_cols*2-1)
    for irow in range(num_rows-1):
        start_idx = start_idx + (num_cols*2-1)
        for ibox in range(num_cols-1):
            ww = ibox+start_idx
            xx = ww+num_cols-1
            yy = xx+1
            zz = yy+num_cols-1
            idx.append([ww,xx,yy,zz])
         
    #populate C with concentration, first two entries are + second two are -, second
            #row is a 1 or -1 for conservation of fluid
    
    for irow in range(0,len(idx)):
        conc_row = irow*2
        aa =  water_isotope_test[0][idx[irow][0:2]]#first two entries, positive
        bb = -water_isotope_test[0][idx[irow][2:4]] #second two entries, negative
        cc = np.concatenate((aa,bb)) #concatenated 
        conc_matrix_test[conc_row,idx[irow]]=cc
        conc_matrix_test[conc_row+1,idx[irow]]=[1,1,-1,-1]
    #    
    #delC
    change_matrix = np.zeros([num_boxes*2,1])
    box_avgs=[]
    for irow in range(num_rows-1):
        for icol in range(num_cols-1):
            aa=[irow,icol] ; bb=[irow,icol+1]
            cc =[irow+1,icol]; dd=[irow+1,icol+1]        
            box_avgs.append(np.mean([iso_grid_trim[aa[0],aa[1]],iso_grid_trim[bb[0],bb[1]],iso_grid_trim[cc[0],cc[1]],iso_grid_trim[dd[0],dd[1]]]))
    box_avgs = np.round(box_avgs,decimals=2)
    #
    change_matrix = np.zeros([num_boxes*2,])
    #change_matrix = np.add(change_matrix,4)
    counter = -2
    for ibox in range(num_boxes):
        counter = counter + 2 #put in change in d18O every other 
        change_matrix[counter] = d18Oinit-box_avgs[ibox]
    
    ## Place boundary condition: edges of hydro cell set to 0 flow, so first column/last column and bottom row
    matrix_shape = conc_matrix_test.shape
    isotope_rows = list(range(0,matrix_shape[0],2)) #row indices where avg isotope values are 
    fluid_rows = list(range(1,matrix_shape[0],2))#row indices where fluid direction/balance are    
    
    top_row_idx = list(range(0,num_cols-1,1))# find the indicies of the top row
    left_edge_idx = list(range(num_cols-1,total_avg_entries-num_rows+2,2*num_cols-1)) # find indices of left side of area
    right_edge_idx = list(np.add(left_edge_idx,(num_cols-1)))
    bottom_row_idx = list(range(total_avg_entries-num_cols+1,total_avg_entries,1)) 
    
    ### Uncomment to say zero flow around edges and bottom of cell ######
    for irow in range(len(fluid_rows)):
        temprow = fluid_rows[irow]
        tempisorow = isotope_rows[irow]
#        conc_matrix_test[temprow,left_edge_idx] = 0 
#        conc_matrix_test[temprow,right_edge_idx] = 0 
        conc_matrix_test[temprow,bottom_row_idx] = 0 
#        conc_matrix_test[tempisorow,left_edge_idx] = 0 
#        conc_matrix_test[tempisorow,right_edge_idx] = 0 
        conc_matrix_test[tempisorow,bottom_row_idx] = 0 
    
    ##use least squares to solve
    #
    #multiply change_matrix by moles O in rock
    change_matrix = np.multiply(change_matrix,moles_O_box)
    #add row on bottom to saw flow in left top is equal to flow out top right
    #conc_matrix_test[-1,:] = 0
    #out_idx = top_row_idx[0:np.int(len(top_row_idx)/2)]
    #in_idx = top_row_idx[np.int(len(top_row_idx)/2):np.int(len(top_row_idx))]
    #conc_matrix_test[-1,out_idx] = 1
    #conc_matrix_test[-1,in_idx] = -1
    #
    #change_matrix[-1,] = 0
    
    
    #solution
    qsolved_test = np.linalg.lstsq(conc_matrix_test,change_matrix)[0]  #least squares
    #qsolved_test = qsolved_test[0]
    ##resid = np.subtract(change_matrix,np.dot(conc_matrix_test,qsolved_test))
    ##tapaered least squares....  
    alpha2=2.05e-4 # 2.75e-4 3e-4 to match EPR (-0.2 permil), 1.45e-4 for no side or bottom flow, 3e-4 for no bottom flow but flow on sides, 3.75e-4 for no flow restritcions 
    conctranspose = np.transpose(conc_matrix_test)
    concdot = np.dot(conctranspose,conc_matrix_test)
    iden_mat = np.identity(len(concdot))
    addinside = np.add(concdot,np.multiply(alpha2,iden_mat))
    invert = np.linalg.inv(addinside)
    aa=np.dot(invert,conctranspose)
    qsolved_taper = np.dot(aa,change_matrix)
    
    
    taper_resid = np.subtract(change_matrix,np.dot(conc_matrix_test,qsolved_taper))
    
    qsolved_test = qsolved_taper
    
#    
#    plt.figure(num=114); plt.clf()
#    plt.plot(qsolved_test,qsolved_taper,'o')
#    plt.xlabel('Least squares standard'); plt.ylabel('Tapered least squares')
#    
#    plt.figure(num=115); plt.clf()
#    plt.plot(change_matrix,np.divide(taper_resid,change_matrix),'o')
#    plt.xlabel('$\Delta$C'); plt.ylabel('Residual %')
#    #plt.ylim([min(np.divide(taper_resid,change_matrix)),max(np.divide(taper_resid,change_matrix))])
#    #plot 
#    
    [xfluid_points,xfluid_cols,xfluid_rows]=matrix_average(x_grid_trim)
    [zfluid_points,zfluid_cols,zfluid_rows]=matrix_average(z_grid_trim)
    
    #get indices of horizontal and vertical vectors 
    vert_ind = []; vert_counter = -(2*num_cols-1)
    for ivert in range(num_rows):
        vert_counter = vert_counter+(2*num_cols-1)
        vert_ind.extend(range(vert_counter,vert_counter+num_cols-1))
    
    vector_components = np.zeros([2,len(qsolved_test)])
    vector_components[0,vert_ind] = np.reshape(qsolved_test[vert_ind],[1,len(vert_ind)]) #add vertical components     
    
    aa = np.where(vector_components[0,:]==0)
    horz_ind=[]
    horz_ind.extend(aa[0])
    vector_components[1,horz_ind] = np.reshape(qsolved_test[horz_ind],[1,len(horz_ind)])
    
#    
#    plt.figure(num=222); plt.clf(); #plt.title('Autoplotted') 
#    #plt.plot(xfl,test_y,'go')
#    
#    plt.contourf(x_grid_trim,z_grid_trim,iso_grid_trim,cmap=plt.cm.BrBG)
#    #cm = plt.cm.get_cmap('RdYlBu')
#    cbar=plt.colorbar()
#    cbar.ax.get_yaxis().labelpad = 15
#    cbar.ax.set_ylabel('$\delta^{18}$O',rotation=270)
#    
#    plt.quiver(xfluid_points,zfluid_points,vector_components[1],vector_components[0])
#    plt.plot(x,z,'w^',markeredgecolor='k',markersize=12)
#    #plt.plot(x_grid_trim,y_grid_trim_temp,'k*')
#    #plt.ylim([-500, 3000])
#    ax = plt.gca() # grab the current axis
#    #ax.set_xticks([1,2,3]) # choose which x locations to have ticks
#    #labels = [0,500,1000,1500,2000,2500,3000,3500]
#    #ax.set_xticklabels(labels) # set the labels to display at those ticks
#    plt.xlabel('Horizontal (m)'); plt.ylabel('Vertical (m)')
#    
    #back out total moles fluid
    
    moles_fluid = np.linalg.norm(qsolved_test)/2
    W_R_predicted = moles_fluid/(num_boxes*moles_O_box)
    moles_fluid_taper = np.linalg.norm(qsolved_taper)
    moles_O_rock = num_boxes*moles_O_box
    W_R_taper = (np.linalg.norm(qsolved_taper)/2)/(num_boxes*moles_O_box)
    water_initial_d18O = (d18Of-d18Oinit)/(W_R_predicted) + d18Of - np.nanmean(epsilon_r_w)
#    water_initial_d18O = (d18Of-d18Oinit)/W_R_predicted + d18Of - (-np.nanmean(epsilon_list_test))
    water_initial_d18O_taper = (d18Of-d18Oinit)/(W_R_taper) + d18Of - (-np.nanmean(epsilon_list_test))
    water_outgoing_d18O = (water_initial_d18O_taper*moles_fluid + moles_O_box*num_boxes*d18Oinit)/moles_fluid
    #Water-rock ratio for each box
    W_R_box = np.zeros([num_boxes])
    for iidx in range(len(idx)):
        fluidsum = np.linalg.norm(qsolved_test[idx[iidx]])/2
        wrtemp = fluidsum/(moles_O_box)
        W_R_box[iidx] = wrtemp
        
    
    # add up vectors for all boxes
    summed_horz_vectors = []
    summed_vert_vectors = []
    for irow in range(len(idx)):
        tempidx=idx[irow]
        summed_vert_vectors.append(np.add(qsolved_test[tempidx[0]],qsolved_test[tempidx[3]]))
        summed_horz_vectors.append(np.add(qsolved_test[tempidx[1]],qsolved_test[tempidx[2]]))
        
    summed_vert_vectors = np.multiply(summed_vert_vectors,0.5)
    summed_horz_vectors = np.multiply(summed_horz_vectors,0.5)
    
    horz_vec_normalized = []
    vert_vec_normalized = []
    vertnorm = []
    horznorm=[]
    for ivec in range(len(summed_horz_vectors)):
        horznorm.append(np.linalg.norm(summed_horz_vectors[ivec]))
        vertnorm.append(np.linalg.norm(summed_vert_vectors[ivec]))
        horz_vec_normalized.append(summed_horz_vectors[ivec]/(horznorm[ivec]**0.5))
        vert_vec_normalized.append(summed_vert_vectors[ivec]/(vertnorm[ivec]**0.5)) 
    #get center point coordinates of all boxes
    x_center_points = np.zeros([num_rows-1,num_cols-1])
    for irow in range(num_rows-1): 
        for icol in range(num_cols-1):
            x_center_points[irow,icol] = np.mean([x_grid_trim[irow,icol],x_grid_trim[irow,icol+1]])
    
    y_center_points = np.zeros([num_rows-1,num_cols-1])
    
    for irow in range(num_rows-1): 
        
        #y_center_points[irow]=irow
        y_center_points[irow] = np.mean([z_grid_trim[num_rows-irow-1],z_grid_trim[num_rows-irow-2]])
        #print(y_center_points)    
        
    #plt.figure(num=76); plt.clf()
    #plt.plot(x_grid_trim,y_grid_trim,'ko')
    #plt.plot(x_center_points,y_center_points,'r^')    
    x_center_list = x_center_points.reshape([len(W_R_box)])
    y_center_list = y_center_points.reshape([len(W_R_box)])
    
    W_R_grid = griddata((x_center_list,y_center_list),W_R_box,(xi,zi),method='linear')
    #
    for irow in range(W_R_grid.shape[0]):
        W_R_grid[irow][np.isnan(W_R_grid[irow])]=np.nanmean(W_R_grid[irow])
    
    W_R_grid[0] = W_R_grid[1]  
    W_R_grid[-1]  = W_R_grid[-2]  
   

    return [iso_grid_trim,x_grid_trim,z_grid_trim,\
            conc_matrix_test,change_matrix,qsolved_taper,moles_fluid_taper,W_R_taper,\
            water_initial_d18O, water_outgoing_d18O,temp_test,temp_grid,moles_O_rock,epsilon_r_w]
