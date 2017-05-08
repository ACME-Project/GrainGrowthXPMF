from fipy import *
import random, math, time, sys, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

def gradGrid(mesh, nx, ny, bins, vx = [], vy = [], scale = 1.0):
	num_bins = len(bins)
	grains = sum(bins)

	if not vx:
		for bin_num in range(num_bins):
			start = bin_num*(nx/num_bins)
			end = (bin_num + 1)*(nx/num_bins)
			for i in range(bins[bin_num]):
				vx.append(random.randrange(start,end))
				vy.append(random.randrange(ny))
	else:
		for i in range(len(vx)):
			vx[i] = int(float(vx[i]*scale))
		for i in range(len(vy)):
			vy[i] = int(float(vy[i]*scale))

	phi_field = [CellVariable(name='Grain-'+repr(n),mesh=mesh) for n in range(grains)]
	diagram = np.zeros((nx, ny), dtype=np.int)
	for x in range(nx):
		for y in range(ny):
			dmin = math.hypot(nx-1, ny-1)
			j = -1
			for i in range(grains):
				d = math.hypot(vx[i]-x, vy[i]-y)
				if d <= dmin:
					dmin = d
					j = i
				diagram[x][y] = j
	listed = []
	for y in range(ny):
		for x in range(nx):
			listed.append(diagram[x][y])
	rav = np.array(listed)
	filename = str(str(nx)+'x'+str(ny)+'.txt')
	with open(filename, 'w') as f:
		for x in range(len(rav)):
			f.write(str(rav[x]) + '\n')
#	curr_id = -1
	with open(str('data'+filename), 'w') as f:
		for x in range(grains):
			f.write(str(x) + '\n')
	for n in range(grains):
		phi_field[n].setValue(0.)
		phi_field[n].setValue(1., where= rav==n)
	points = zip(vx, vy)
	print points
	return phi_field, points

min_size = 700

aspect_ratio = 3.0 #(height/width) of individual bins

bin_min = 50
bias = 5
num_bins = 9
bins = []
for i in range(num_bins):
#	bins.append(int((grains/num_bins) + bias*(num_bins/2 - i)))
	bins.append(bin_min + bias*i)

bins = [110,110,110,110,110,110,110,110,110] #bias = 0
#bins = [90,95,100,105,110,115,120,125,130] #bias = 5
#bins = [70,80,90,100,110,120,130,140,150] #bias = 10
#bins = [50,65,80,95,110,125,140,155,170] #bias = 15
#bins = [30,50,70,90,110,130,150,170,190] #bias = 20, below threshold of 700
bin_area = 170*min_size

bin_x = int(np.sqrt(bin_area/aspect_ratio))
bin_y = int(bin_x*aspect_ratio)

nx = len(bins)*bin_x
ny = bin_y

print "Bins of",bins,"initializing on grid of",nx,"by",ny

dx = dy = 1

mesh = Grid2D(nx=nx, dx=dx, ny=ny, dy=dy)
phi_field = []

phi_field, points = gradGrid(mesh, nx, ny, bins)

grains = len(phi_field)

viewer_field = CellVariable(name='Grain-ID', mesh=mesh)
viewer_field.value = 0.0
for n in range(grains):
    viewer_field.value += (phi_field[n] > 0.99)*(n+1)
outputname = str(nx) + 'x' + str(ny) + '_init'
viewer = Matplotlib2DViewer(vars=viewer_field)
viewer.plot(filename=outputname + '.png')
