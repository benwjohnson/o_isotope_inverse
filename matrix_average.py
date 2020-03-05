#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 20 10:20:13 2018

@author: benjohnson
"""

#function to take a matrix input and:
# 1. Calculate average value in between all entries 
# Example, for matrix [3 4 5; 1 2 3; 3 4 5], this returns a matrix (col_avgs) 
# for averages between columns for all rows [3.5, 4.5; 1.5, 2.5; 3.5, 4.5] a matrix
# (row_avgs) for values between rows for all colums [2 3 4; 2 6 4]
# 2. packages col_avgs and row_avgs as a vector, where rows of col_avgs alternate with rows of row_avgs
# User puts matrix to be processed in matrix_average(), then defines the name of the output vectors and col/row avgs

def matrix_average(A):
    import numpy as np
    
    grid_dim = A.shape
    num_rows = grid_dim[0]
    num_cols = grid_dim[1]
    total_raw_entries = A.size
    total_avg_entries = (num_cols-1)*num_rows + (num_rows-1)*num_cols# total number of averages


    col_avgs = []
    # calculate means between columns in all rows
    for irow in range(grid_dim[0]):
        for icol in range(grid_dim[1]-1):
            col_avgs.append(np.mean([A[irow,icol],A[irow,icol+1]]))
        
    col_avgs = np.round(col_avgs,decimals=2)
    col_avgs = col_avgs.reshape([num_rows,num_cols-1])

    #calculate means between row entries in all columns
    row_avgs = []  
    for irow in range(grid_dim[0]-1):  
        for icol in range(grid_dim[1]):    
            row_avgs.append(np.mean([A[irow,icol],A[irow+1,icol]]))
            
    row_avgs = np.round(row_avgs,decimals=2)
    row_avgs = row_avgs.reshape([num_rows-1,num_cols])

    #put means together into a single vector
  
    avg_list_auto = np.zeros([1,total_avg_entries]) 
   
    counter_col1 = -(row_avgs.shape[1]+col_avgs.shape[1])

    for icol in range(col_avgs.shape[0]):
        counter_col1 = counter_col1+row_avgs.shape[1]+col_avgs.shape[1]
        counter_col2 = counter_col1+col_avgs.shape[1] 
        avg_list_auto[0][counter_col1:counter_col2] = col_avgs[icol]
   
    counter_row1 = -row_avgs.shape[1]
    for irow in range(row_avgs.shape[0]):
        counter_row1 = counter_row1 + row_avgs.shape[1]+col_avgs.shape[1] 
        counter_row2 = counter_row1 + row_avgs.shape[1]
        avg_list_auto[0][counter_row1:counter_row2] =  np.transpose(row_avgs)[:,irow] 
    
    
    output_name = avg_list_auto  
    
    return [output_name, col_avgs, row_avgs]