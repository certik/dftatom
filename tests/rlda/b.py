from numpy import fromstring, select
from pylab import plot, savefig, clf, semilogx, legend

def ctp(V, E):
    for i in range(len(V)-1, 1, -1):
        if V[i] -E >= 0: continue
        return i
    stop

f = open("log.txt")
R = fromstring(f.readline(), sep=" ")
E = float(f.readline())
l = int(f.readline())
V = fromstring(f.readline(), sep=" ")
Veff = fromstring(f.readline(), sep=" ")
#Veff = l*(l+1) + 2*(V-E)*R**2
Z = 66
V = -Z/R
Veff = l*(l+1)/(2*R**2) + V
print ctp(Veff, E)
stop
Veff = Veff - E
print "potential min:", min(Veff), "E=", E
Veff = select([Veff < 0], [Veff])
semilogx(R, Veff, "ko-", label="from V")
legend()
savefig("c.png")
