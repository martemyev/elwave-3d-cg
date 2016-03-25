#!/usr/bin/python

import struct

prop_main  = 2000.0
prop_layer = 1500.0

# size of the domain
sx = 1000.0
sy = 1000.0
sz = 1000.0

# coordinates of the layer
Y0 = 400.0
Y1 = 600.0

# grid
nx = 100
ny = 100
nz = 100

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

