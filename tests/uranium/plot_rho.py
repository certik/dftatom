from pylab import plot, show, legend
from numpy import array
from h5py import File


data = File("data.h5")
iter = 2
R = array(data["/%04d/R" % iter])
rho = array(data["/%04d/rho" % iter])
Vtot = array(data["/%04d/V_tot" % iter])

Zeff = -Vtot * R

#for i in range(1, 19):
#    P = array(data["/%04d/P%04d" % (iter, i)])
#    plot(R, P, label="P%04d" % i)
for i in range(1, 11):
    iter = i
    R = array(data["/%04d/R" % iter])
    rho = array(data["/%04d/rho" % iter])
    plot(R, rho*R**2, label="iter=%d" % iter)
legend()
show()
