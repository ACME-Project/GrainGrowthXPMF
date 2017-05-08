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
	#	for y in range(ny):
	#		for x in range(nx):
	#			bin_num = int(x / (nx/num_bins))
	#			weight = bins[bin_nnum]
	#			diagram[x][y] = int(0 + 1*(x >= nx/2))
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
#	rav = diagram.ravel()
#	print type(rav)
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

#bins = [12,14,16,18,20,22,24,26,28]
#bins = [20,20,20,20,20,20,20,20,20]
#bins = [16,17,18,19,20,21,22,23,24]
#bin_area = max(bins)*min_size

bins = [110,110,110,110,110,110,110,110,110]
#bins = [90,95,100,105,110,115,120,125,130]
#bins = [70,80,90,100,110,120,130,140,150]
#bins = [50,65,80,95,110,125,140,155,170]
#bins = [30,50,70,90,110,130,150,170,190] #below threshold of 700
bin_area = 170*min_size

bin_x = int(np.sqrt(bin_area/aspect_ratio))
bin_y = int(bin_x*aspect_ratio)

nx = len(bins)*bin_x
ny = bin_y

print "Bins of",bins,"initializing on grid of",nx,"by",ny

#nx = 4000
#ny = 1000
dx = dy = 1

mesh = Grid2D(nx=nx, dx=dx, ny=ny, dy=dy)
phi_field = []

#scale = float(657. / 612.)

#vx = [0, 15, 60, 29, 41, 6, 3, 2, 58, 65, 17, 50, 14, 32, 1, 16, 28, 31, 49, 46, 129, 127, 80, 84, 99, 132, 73, 127, 70, 102, 99, 125, 92, 84, 100, 86, 107, 83, 84, 76, 177, 202, 138, 187, 143, 138, 157, 202, 196, 195, 181, 148, 189, 176, 181, 186, 175, 194, 156, 184, 264, 233, 225, 260, 212, 206, 260, 267, 241, 215, 250, 240, 220, 252, 242, 235, 235, 240, 268, 233, 331, 319, 327, 325, 272, 337, 332, 307, 272, 299, 291, 331, 293, 311, 313, 331, 301, 291, 323, 303, 343, 400, 393, 375, 357, 376, 402, 402, 391, 360, 358, 383, 374, 372, 391, 362, 343, 345, 359, 355, 419, 443, 459, 421, 437, 424, 426, 460, 444, 416, 464, 420, 439, 465, 434, 438, 473, 413, 436, 409, 512, 526, 501, 484, 531, 526, 517, 487, 514, 502, 542, 537, 492, 537, 521, 504, 525, 529, 522, 526, 604, 563, 578, 586, 564, 558, 596, 604, 566, 547, 577, 573, 587, 610, 573, 607, 562, 591, 609, 561]

#vy = [166, 47, 52, 94, 20, 33, 137, 98, 59, 141, 139, 122, 12, 155, 150, 174, 174, 26, 94, 73, 52, 202, 131, 43, 156, 92, 150, 90, 106, 28, 202, 143, 178, 52, 198, 166, 36, 119, 201, 162, 101, 19, 150, 71, 180, 96, 177, 120, 161, 169, 181, 48, 93, 175, 157, 123, 30, 75, 34, 0, 118, 40, 110, 66, 85, 91, 2, 46, 87, 165, 94, 161, 44, 103, 190, 199, 155, 160, 84, 178, 102, 153, 88, 38, 36, 76, 91, 65, 179, 168, 157, 111, 189, 174, 198, 9, 202, 186, 156, 167, 53, 56, 170, 173, 116, 162, 136, 58, 53, 29, 11, 169, 62, 15, 20, 143, 198, 194, 102, 171, 104, 180, 160, 186, 202, 7, 37, 170, 149, 23, 86, 30, 17, 130, 136, 41, 173, 45, 97, 1, 177, 30, 197, 202, 53, 74, 170, 77, 135, 6, 17, 117, 9, 93, 173, 129, 97, 3, 9, 32, 103, 59, 19, 2, 81, 5, 25, 6, 23, 145, 65, 159, 83, 55, 149, 137, 90, 172, 61, 16]
#print scale
#phi_field, points = gradGrid(mesh, nx, ny, bins, vx, vy, scale)

phi_field, points = gradGrid(mesh, nx, ny, bins)

#vx = []
#vy = []
#for i in range(len(points)):
#	vx.append(points[i][0])
#	vy.append(points[i][1])

#print vx
#print vy
grains = len(phi_field)

viewer_field = CellVariable(name='Grain-ID', mesh=mesh)
viewer_field.value = 0.0
for n in range(grains):
    viewer_field.value += (phi_field[n] > 0.99)*(n+1)
outputname = str(nx) + 'x' + str(ny) + '_init'
#TSVViewer(vars=phi_field).plot(filename=outputname + '.txt')
viewer = Matplotlib2DViewer(vars=viewer_field)
viewer.plot(filename=outputname + '.png')
