from pylab import plot, savefig, legend, grid, xlabel, ylabel, title
from numpy import loadtxt

for p in [8, 7, 6, 5, 4, 3]:
    data = loadtxt("data%d.dat" % p, dtype="int")
    Z = data[:, 0]
    N = data[:, 1]
    plot(Z, N, label="1e-%d" % p)
legend(loc="upper left")
xlabel("Z")
ylabel("N")
title("Lowest N needed for the given accuracy for each atom")
grid()
savefig("plot_n.pdf")
