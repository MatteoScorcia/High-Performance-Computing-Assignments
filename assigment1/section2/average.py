#!/usr/bin/env python3
import numpy
from sys import argv

numpy.set_printoptions(suppress=True)

# arg1 = directory
# arg2 = network

files = []

basename = argv[1] + "/result"
for i in range(0, 50):
	name = basename + "_{}_".format(i) + argv[2]
	files.append(name)

data = numpy.genfromtxt(files[0], delimiter=',')
for f in files[1:]:
    data += numpy.genfromtxt(f, delimiter=',')

# row is: bytes, reps, time, bw

data /= len(files)
for row in data:
	for num in row:
		string = "%.2f"%num
		print(string, end=',')
	print()
