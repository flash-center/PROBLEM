import numpy as np

a = np.zeros((5,5))
b = np.zeros((5,5))
b[:,:] = 1

with open('test.txt', 'wb') as f:
    f.write(bytearray('# Hello World', 'utf-8'))
    np.savetxt(f, a, delimiter=',')
    np.savetxt(f, b, delimiter=',')

