#!/usr/bin/python

import struct

prop_main  = 3500.0
prop_layer = 2000.0

# size of the domain
sx = 1000.0
sy = 1000.0
sz = 1000.0

# coordinates of the layer
Y0 = 300.0
Y1 = 500.0

# grid
nx = 50
ny = 50
nz = 50

tol = 1e-8

hx = sx / nx; print 'hx = ', hx
hy = sy / ny; print 'hy = ', hy
hz = sz / nz; print 'hz = ', hz


prop_values = []
for iz in range(nz):
  for iy in range(ny):
    y = (iy+0.5)*hy
    for ix in range(nx):
      if y > Y0 and y < Y1:
        prop_values.append(prop_layer)
      else:
        prop_values.append(prop_main)


# write the properties into a binary file with a single precision
out = open('properties.bin', 'wb')
s = struct.pack('f'*len(prop_values), *prop_values)
out.write(s)
out.close()

