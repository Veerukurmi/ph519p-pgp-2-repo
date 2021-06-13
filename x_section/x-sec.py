# This program is used for calculation of the phase-shift 
# Code by Apoorav Singh Deo
# Github: https://github.com/apoorav-singh




# Importing pandas to export data from text file 
# Importing Matplotlib.pyplot to draw graphs
# Importing numpy to to make faster array like structure to store data
import pandas as pd
import matplotlib.pyplot as plt 
import numpy as np
import math

# Loading the text file
df_txt_0 = pd.read_csv('ph_sh_0.txt',delimiter='\s+')
df_txt_1 = pd.read_csv('ph_sh_1.txt',delimiter='\s+')
df_txt_2 = pd.read_csv('ph_sh_2.txt',delimiter='\s+')
df_txt_3 = pd.read_csv('ph_sh_3.txt',delimiter='\s+')
df_txt_4 = pd.read_csv('ph_sh_4.txt',delimiter='\s+')
df_txt_5 = pd.read_csv('ph_sh_5.txt',delimiter='\s+')
df_txt_6 = pd.read_csv('ph_sh_6.txt',delimiter='\s+')
df_txt_7 = pd.read_csv('ph_sh_7.txt',delimiter='\s+')
df_txt_8 = pd.read_csv('ph_sh_8.txt',delimiter='\s+')
df_txt_9 = pd.read_csv('ph_sh_9.txt',delimiter='\s+')
df_txt_10 = pd.read_csv('literature.txt',delimiter='\s+')

#
# Printing the data file out to check if data is correctly loaded or not
# This can be later removed
# print(df_txt)
E_lit = np.array(df_txt_10.iloc[0:130,0])/27.211324570273
X_lit = np.array(df_txt_10.iloc[0:130,1])*100
#print(X_lit)

grid = 50 #1000000
#print(grid)

E = np.array(df_txt_0.iloc[0:grid,0])
phi_0 = np.array(df_txt_0.iloc[0:grid,1])
phi_1 = np.array(df_txt_1.iloc[0:grid,1])
phi_2 = np.array(df_txt_2.iloc[0:grid,1])
phi_3 = np.array(df_txt_3.iloc[0:grid,1])
phi_4 = np.array(df_txt_4.iloc[0:grid,1])
phi_5 = np.array(df_txt_5.iloc[0:grid,1])
phi_6 = np.array(df_txt_6.iloc[0:grid,1])
phi_7 = np.array(df_txt_7.iloc[0:grid,1])
phi_8 = np.array(df_txt_8.iloc[0:grid,1])
phi_9 = np.array(df_txt_9.iloc[0:grid,1])


k = [0]*(grid-1)

i=0
#print(k)

while i<(grid-1):
    k[i] = math.sqrt(2*E[i])
    #print(k[i])
    i+=1

X = [0]*(grid)

j=0
while j<(grid-1):
    #X[j] = ((4*(22/7))/k[j]**2)*(1*(math.sin(phi_0[j]))**2 + 3*(math.sin(phi_1[j]))**2 + 5*(math.sin(phi_2[j]))**2 + 7*(math.sin(phi_3[j]))**2 + 9*(math.sin(phi_4[j]))**2 + 11*(math.sin(phi_5[j]))**2 + 13*(math.sin(phi_6[j]))**2 + 15*(math.sin(phi_7[j]))**2)
    X[j] = ((4*(22/7))/k[j]**2)*(15*(math.sin(phi_7[j]))**2)
    #print(X[j])
    j+=1


plt.figure(figsize=(8,5), dpi=300)
plt.plot(E[0:(grid-1)],X[0:(grid-1)],"b->")
#plt.plot(E_lit,X_lit,"r->")
#plt.legend(["Computed Value","Literature Value"])
plt.legend(["l=7"])
plt.xlabel('E (a.u.)')
#plt.ylabel(r'$\delta_l$')
plt.ylabel("Cross-section (a.u.)")
plt.grid()
plt.title("Cross-Section")
plt.show()
#plt.text(1,0.2, r'step size $(n) = 100000$')
plt.savefig('x-sec_7.png', dpi=300)

# Saving the csv file for matlab computation
# Creating a dictionary to save data
dict_1 = {'x':E[0:(grid-1)],'psi':X[0:(grid-1)]}
df_1 = pd.DataFrame(dict_1)

# Saving the csv file
df_1.to_csv('cross_section_7.csv')