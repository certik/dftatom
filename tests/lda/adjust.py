from numpy import loadtxt, savetxt

data = loadtxt("data8.dat", dtype="int")
data[:, 1] += 1000
savetxt("data8.dat.new", data, fmt="%d")
