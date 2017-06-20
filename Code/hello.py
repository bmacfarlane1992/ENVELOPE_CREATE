# hello.py


import numpy as np          # Import Python system modules
import math                 # When names are long and convoluted, can be called "as"
import os
import matplotlib.pyplot as plt

import bye                  # Import user modules, which should be in cwd, when they aren't in PYTHONPATH

arch_dir = os.getcwd()         # get cwd, and output to variable printed to terminal
print arch_dir

# variables

pi = math.pi            # Ensure all module specified routines are called AFTER importing
print pi

x = np.array([0,1,2,3])         # Define array - useful to use numpy (as np) arrays for multi-D analysis

for i in range(0, len(x)):
    print x[i]
                                # Two ways of using for loop - second way removes len() requirement
for i in x:                     # but I think the x[i] notation is more intuitive
    print i

print "Hello world"             # Output something to terminal - can be either print "blah" or print("blah")

x[2], pi = bye.byebye(x[2],pi)      # Call bye module, specifically the byebye function which reads 2 variables
                                    # and modifies them in some way, returning the values back
print x[2], pi

print "Still here mofo"             # Printed when exit() condition isn't invoked in bye module

y = np.array([0,1,2,3])

plt.scatter(x,y,color='r')                      # With imported matplotlib.pyplot as plt, plot x and y arrays
plt.savefig(arch_dir+'/helloyoucunt.png')
plt.clf()                                       # Make sure you use plt.clf() to prevent masking in any future plots
                                                # (Though a useful thing when you want to overplot stuff!)
x2 = [1.3,1.5,1.7]
y2 = [1.5, 1.7, 1.9]

plt.scatter(x2,y2,color='g')
plt.scatter(x,y, color = 'r')                   # Demonstration of overplotting applied
plt.savefig(arch_dir+'/harambe.eps')
plt.clf()

bye.byebye2()                       # Call byebye2() function in bye module (needs no variables)

exit()                          # Not explicitly needed, but useful to exit Python
