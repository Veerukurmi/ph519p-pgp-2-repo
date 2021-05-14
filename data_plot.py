# Importing pandas to export data from text file 
# Importing Matplotlib.pyplot to draw graphs
# Importing numpy to to make faster array like structure to store data
import pandas as pd
import matplotlib.pyplot as plt 
import numpy as np

# Loading the text file
df_txt = pd.read_csv('wave_data.txt',delimiter='\s+')
#
# Printing the data file out to check if data is correctly loaded or not
# This can be later removed
# print(df_txt)
grid = 1000000

psi = np.array(df_txt.iloc[0:grid,1])
x = np.array(df_txt.iloc[0:grid,0])


plt.figure(figsize=(8,5), dpi=300)
plt.plot(x,psi)
plt.xlabel('x-position')
plt.ylabel(r'$\psi(r)$')
plt.grid()
plt.title("Wavefunction's radial part only")
plt.text(1,0.2, r'step size $(n) = 100000$')
plt.savefig('wave_check_2.png', dpi=300)

# Saving the csv file for matlab computation
# Creating a dictionary to save data
dict_1 = {'x':x,'psi':psi}
df_1 = pd.DataFrame(dict_1)

# Saving the csv file
#df_1.to_csv('wavefunction.csv')