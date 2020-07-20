# -*- coding: utf-8 -*-
"""
Created on Tue Jun 16 16:42:58 2020

@author: Dell
"""

import numpy as np
import pandas as pd

###############################################################################
"""                   importing data from csv file                          """

heterogeneous_reservoir_data = pd.read_csv('./heterogeneous_reservoir.csv')
print(heterogeneous_reservoir_data.count())
print("\n")
print(heterogeneous_reservoir_data)
k = np.array(heterogeneous_reservoir_data['permeability(md)'])

###############################################################################
L = 10000
print("\n\tThe lenght of the reservoir is ", str(L))
number_nodes = 4
print("\n\tBlock node in the reservoir is ", str(number_nodes))
P0 = 1000
print("\n\tThe intial pressure of the reservoir is ", str(P0))
P_left = 2000
print("\n\tThe pressure at the left boundary of the reservoir is ", str(P_left))
dx = L/number_nodes
print("\n\tThe lenght of Blocks in the reservoir is ", str(dx))

porosity = 0.2
print("\n\tthe porosity value of the reservoir is ", str(porosity))

permeability  = 50
print("\n\tthe permeability(md) value of the reservoir is ", str(permeability))

viscosity = 1
print("\n\tthe viscosity value is ", str(viscosity))

area = 200000
print("\n\tCross sectional area of the reservoir ", str(area))

compressibility = 1*10**(-6)
print("\n\tcompressibility of the reservoir is ", str(compressibility))

Bw = 1
print("\n\t water formation volume factor is " +str(Bw)+ " rb/stb" )


###############################################################################
###############################################################################
#################### final time for simulation is  ############################

t_final = 3
print("\n\t the reservoir simulation time in days is  " +  str(t_final) + "days")

#################### time step  ###############################################

dt_increment = 1
print("\n\t the reservoir simulation incremental time step in days is "+ str(dt_increment)+ "day")

###############################################################################
###############################################################################
###############################################################################
############            pressure and  boundary condition          #############

pressure_previous = np.ones([number_nodes,1])*P0
print("\n############## pressure distribution ################\n")
print("pressure distribution at day 0 is\n", str(pressure_previous))

boundary_condition_array =  np.zeros([number_nodes,1])
print("\n boundary_condition_array at day 0 is\n", str(boundary_condition_array))

###############################################################################
###############################################################################

print("########################")

def permeability(i,k):
    perm = 2*(( 1/k[i] + 1/k[i+1] )**(-1))
    return perm

transmissibility_1 = np.ones([5,1])
print(transmissibility_1)
print("########################")
      
for i in range(0,3):
    print(i)
#    print(permeability(i,k))
    k_permeability = permeability(i,k)
    print(k_permeability)
    transmissibility_1[0] =  2*( k[0]*area )/( Bw*viscosity*dx )
    transmissibility_1[i+1] = ( k_permeability*area )/( Bw*viscosity*dx )
    transmissibility_1[4] =  0
    #print(transmissibility_1[i])
        
print(transmissibility_1)
print("########################")
      
print("\n number = ", number_nodes)
print("\n ")

transmissibility_2 = np.ones([number_nodes+1,1])
print(transmissibility_2)

print("########################")
print(k)      

for i in range(0,number_nodes-1):    
    print(i)
    print(k[i])
    transmissibility_2[0] = 2*( k[0]*area )/( Bw*viscosity*dx )
    print("transmissibility of block " + str ( i +1/2 ) + "hai")
#    k_permeability = permeability(i,k)
#    print(k_permeability)
#    transmissibility_2[i+1] = ( k_permeability*area )/( Bw*viscosity*dx )
    transmissibility_2[i+1] = ( permeability(i,k)*area )/( Bw*viscosity*dx )
transmissibility_2[number_nodes] = 1000    

print("\n transibillity matrix = \n",transmissibility_2)
print("########################")
    

      
    




        
###############################################################################
###############################################################################
###############################################################################
print("\n############## transmisibility_matrix ################\n")

transmisibility_matrix = np.zeros([number_nodes , number_nodes]) 
print("\n transmisibility_matrix is\n", str(transmisibility_matrix))

for i in range(1,number_nodes,1):
    print("i value" , i)
#    print(transmissibility_2[i])

    transmisibility_matrix[i][i] = transmissibility_2[i] + transmissibility_2[i+1]  
    transmisibility_matrix[i][i-1] = - transmissibility_2[i]
    transmisibility_matrix[i-1][i] = - transmissibility_2[i]

transmisibility_matrix[0][0] =  transmissibility_2[0] +  transmissibility_2[1] 
transmisibility_matrix[number_nodes-1][number_nodes-1] = transmissibility_2[number_nodes-1]
print("\n transmisibility_matrix is\n", str(transmisibility_matrix))









###############################################################################
print("\n############## B_matrix ################\n")
    
B_matrix = np.zeros([number_nodes , number_nodes])
print("\n B_matrix is\n", str(B_matrix))

B = porosity*compressibility*area*dx

for i in range(1,number_nodes,1):
    #print("i value" , i)
    B_matrix[i][i] = B
    B_matrix[0][0]= B
    B_matrix[number_nodes-1][number_nodes-1]= B
print("\n B_matrix is\n", str(B_matrix))

###############################################################################
print("\n############## Q_matrix ################\n")
    
Q_matrix = np.zeros([number_nodes , 1])
print("\n Q_matrix is\n", str(Q_matrix))

Q = 2000*1600*6.33*10**(-3)

Q_matrix[0][0]= Q
#Q_matrix[number_nodes-1][number_nodes-1]= B
print("\n Q_matrix is\n", str(Q_matrix))

###############################################################################
print("\n############## N_plus_1_matrix ################\n")

inverse_dt =  ( 1 / dt_increment )

N_plus_1_matrix = ( 6.33*10**(-3)*transmisibility_matrix + inverse_dt*B_matrix )

print("\n N_plus_1_matrix is\n", str(N_plus_1_matrix))

###############################################################################
print("\n############## N__matrix ################\n")

inverse_dt =  ( 1 / dt_increment )

N__matrix = (inverse_dt*B_matrix )

print("\n N__matrix is\n", str(N__matrix))


print("\n",pressure_previous)
for k in range(0, t_final, 1):
    
    #print("\n time step value",k)
#    boundary_condition_array[0][0] = 2*neta*P_left
    pressure_previous = np.dot(np.linalg.inv(N_plus_1_matrix), (np.dot(N__matrix , pressure_previous) + Q_matrix))
    print("\n",pressure_previous)

print("\nfinal pressure \n", pressure_previous )










