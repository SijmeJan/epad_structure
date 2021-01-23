import numpy as np
from scipy.integrate import ode
import matplotlib.pyplot as plt

class PolyPlanet():
    def __init__(self, rho0, c, n):
        self.rho0 = rho0
        self.c = c
        self.n = n
        self.G = 6.67408e-11

    def _rho(self, P):
        if P < 0:
            return self.rho0

        return self.rho0 + self.c*np.power(P, self.n)

    def _f(self, t, y):
        M = y[0]
        P = y[1]
        rho = self._rho(P)
        f1 = 4.0*np.pi*t*t*rho   # dM/dr
        f2 = -self.G*M*rho/t/t   # dP/dr

        #print(t, M, P)
        return [f1, f2]

    def _solout(self, t, y):
        if y[1] < 1.0e-2:
            return -1
        return 0

    def solve(self, Pc):
        r = ode(self._f).set_integrator('dopri5')
        r.set_initial_value([0.0, Pc], 1.0e-10)
        r.set_solout(self._solout)

        res = r.integrate(1.0e13)
        return r.t, res

    def single_solution(self, Pc, radii):
        r = ode(self._f).set_integrator('dopri5')
        r.set_initial_value([0.0, Pc], 1.0e-10)
        r.set_solout(self._solout)

        res = 0.0*radii
        for i in range(0, len(radii)):
            res[i] = self._rho(r.integrate(radii[i])[1])
        return res


    def manifold(self, Pc):
        mass = 0.0*Pc
        radius = 0.0*Pc

        for i in range(0, len(Pc)):
            t, y = self.solve(Pc[i])
            radius[i] = t
            mass[i] = y[0]

        return mass, radius


## p = PolyPlanet(rho0=1460.0, c=0.00311, n=0.513)
## pressure = np.logspace(8, 17, 100)
## mass, radius = p.manifold(pressure)

#fig, ax = plt.subplots()
#print(fig.canvas)

#plt.xscale('log')
#plt.yscale('log')
#plt.plot(mass/me, radius/re)

#p = PolyPlanet(rho0=4260.0, c=0.00127, n=0.549)
#pressure = np.logspace(8, 17, 100)
#mass, radius = p.manifold(pressure)
#plt.plot(mass/me, radius/re)


# python 3.x
#import xml.etree.ElementTree as ET, urllib.request, gzip, io
#url = "https://github.com/OpenExoplanetCatalogue/oec_gzip/raw/master/systems.xml.gz"
#oec = ET.parse(gzip.GzipFile(fileobj=io.BytesIO(urllib.request.urlopen(url).read())))

#oec_masses = []
#oec_radii = []
#oec_names = []
# Output mass and radius of all planets
#for planet in oec.findall(".//planet"):
#    if (planet.findtext("mass") is not None and
#        planet.findtext("radius") is not None and
#        planet.findtext("mass") != '' and
#        planet.findtext("radius") != ''):

#        m = float(planet.findtext("mass"))*mjup
#        r = float(planet.findtext("radius"))*rjup

        #plt.plot([m/me], [r/re], marker='o', color='b')

#        oec_masses.append(m/me)
#        oec_radii.append(r/re)
#        oec_names.append(planet.findtext("name"))

#sp = ScatterPlot(fig, oec_masses, oec_radii, oec_names)

#plt.show()
