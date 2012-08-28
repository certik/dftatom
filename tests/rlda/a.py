from numpy import fromstring
from pylab import plot, savefig, clf, semilogx

f = open("log.txt")
i = 1
while 1:
    l = f.readline()
    if not l: break
    R = fromstring(l, sep=" ")
    P = fromstring(f.readline(), sep=" ")
    Q = fromstring(f.readline(), sep=" ")
    idx = int(f.readline())
#    clf()
#    semilogx(R, P)
#    savefig("/tmp/a%04d.png" % i)
#    i += 1

clf()
plot(R[:idx], P[:idx], "ko-")
savefig("b.png")
