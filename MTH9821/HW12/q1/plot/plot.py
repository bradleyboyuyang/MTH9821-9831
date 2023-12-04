import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# Read the data from the csv file
exact = 4.3756
data = pd.read_csv('output.csv', names=['N', 'Value', 'Error']).set_index('N')

# get the values of index N when the error first decrease and then increase

# Find the local minima indices
local_min_indices = np.where((data['Error'].shift(1) > data['Error']) & (data['Error'].shift(-1) > data['Error']))[0]

# Print the local minima indices
print("Local Minima Indices: ", local_min_indices+data.index[0])
# print value of the local minima indices
print("Local Minima Values: ", data.iloc[local_min_indices]['Error'])

# plot the Error with respect to N
# generate ax object with figsize 10x10
data['Error'].plot(figsize=(12, 6), title='Error vs. N', color='b')
# plot the exact value as a horizontal line
plt.axhline(0, color='r', linestyle='--')
plt.savefig('error.png',dpi=200) 

# scatterplot the Error with respect to N = 23, 51, 92, 144, 207, 282, 369, 467, 577, 698, 831, 975
data.loc[[23, 51, 92, 144, 207, 282, 369, 467, 577, 698, 831, 975], 'Error'].plot(style='o', title='Error vs. N', color='Brown')
plt.savefig('error2.png',dpi=200)    

# data.loc[[23, 51, 92, 144, 207, 282, 369, 467, 577, 698, 831, 975], 'Error'].plot(figsize=(12, 6), title='Error vs. N (optimal)', color='Brown')
# plt.savefig('error3.png', dpi=200)