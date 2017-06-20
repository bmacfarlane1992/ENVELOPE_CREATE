import numpy as np

x = [1,2,3,4,[5.,4,3,2],6]
val = x[4][0]
val = val/2
print val, x[0]

x = np.array([[1],[2],[3],[4],[5.,4,3,2],[6]])
val = x[4][0]
val = val/2
print val,x, x[4][0]
