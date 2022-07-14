"""
Debris disk class
@author: Eve J. Lee
Feb. 15th 2016
"""

from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as pl
import fn_draw as fd
import betadistrb as bd
import numpy as np
import numpy.random as nr
import fn_freq as ff
import read_disk_input
import prob_position_orbit as ppo
import consts
import pdb
import math


class DebrisDisk:
    def __init__(self, inputfile, freeelem=False, freeelemtxt=''):
        self.inputdata = read_disk_input.ReadInput(inputfile)
        self.Mstar = self.inputdata["Mstar"] * consts.Msun
        self.age = self.inputdata["Age"] * 1e6 * consts.s2yr
        self.Mp = []
        self.ap = []
        self.ep = []
        self.Ip = []
        self.Omegap = []
        self.omegap = []
        self.freeelem = freeelem
        self.freeelemtxt = freeelemtxt
        self.a_all = []
        self.e_all = []
        self.I_all = []
        self.Omega_all = []
        self.omega_all = []
        self.beta_all = []
        self.beta_bounded = False
        self.verbose = False
        self.print_every_x_dust = 0
        if "verbose" in self.inputdata:
            self.verbose = self.inputdata["verbose"]
            if "print_every_x_dust" in self.inputdata:
                self.print_every_x_dust = self.inputdata["print_every_x_dust"]

    def AddSinglePlanet(self, manual=False, Mp=4., ap=5., ep=0.25, Ip=0., Omegap=0., omegap=0.):
        # Add a planet
        # Mp: planet mass in Mearth
        # ap: planet semi-major axis in AU
        # ep: planet eccentricity
        # Ip: initial mutual inclination IN DEGREES
        # Omegap: nodal angle IN DEGREES
        # omegap: argument of periapse IN DEGREES
        if manual:
            self.Mp.append(Mp)
            self.ap.append(ap)
            self.ep.append(ep)
            self.Ip.append(Ip)
            self.Omegap.append(Omegap)
            self.omegap.append(omegap)
        else:
            self.Mp.append(self.inputdata["Mp"])
            self.ap.append(self.inputdata["ap"])
            self.ep.append(self.inputdata["ep"])
            self.Ip.append(self.inputdata["Ip"] * np.pi / 180.)
            self.Omegap.append(self.inputdata["Omegap"] * np.pi / 180.)
            self.omegap.append(self.inputdata["omegap"] * np.pi / 180.)

        self.MakePlanetArray()

    def MakePlanetArray(self):
        self.Mp = np.array(self.Mp)
        self.ap = np.array(self.ap)
        self.ep = np.array(self.ep)
        self.Ip = np.array(self.Ip)
        self.Omegap = np.array(self.Omegap)
        self.omegap = np.array(self.omegap)

    def ComputeParentSingle(self, manual=False, Nparticles=500, amin=4., amax=6., I0=5. * (np.pi / 180.), e0=0.2,
                            Omega0=0., omega0=0., random=False):
        # Compute p, q, h, k of parent bodies
        # Single planet only
        # Nparticles: number of parent bodies
        # amin, amax: min, max semi-major axes of parent bodies
        # I0, e0: initial mutual inclination and eccentricities
        # Omega0, omega0: initial nodal angle and argument of periapse
        print("Computing Parent Orbits (single planet)...")
        if self.freeelem:
            self.a, self.I0, self.e0, self.Omega0, self.omega0 = np.loadtxt(self.freeelemtxt, unpack=True)
        elif manual:
            self.a = np.linspace(amin, amax, int(self.inputdata["Nparticles"]))
            self.I0 = I0
            self.e0 = e0
            self.Omega0 = Omega0
            self.omega0 = omega0
        else:
            if self.inputdata["InnerPlanet"] == 1:
                amin = self.ap * (1 + 2.0 * (self.Mp * consts.Mearth / self.Mstar) ** (2. / 7.))
                if self.inputdata["OuterPlanet"] == 0:
                    # amax = amin/(1-self.inputdata["awidth"])
                    amax = amin * (1. + self.inputdata["awidth"])
            if self.inputdata["OuterPlanet"] == 1:
                amax = self.ap * (1 - 2.0 * (self.Mp * consts.Mearth / self.Mstar) ** (2. / 7.))
                if self.inputdata["InnerPlanet"] == 0:
                    amin = amax * (1 - self.inputdata["awidth"])

            self.a = np.linspace(amin, amax, int(self.inputdata["Nparticles"]))

            if not (self.inputdata["Random"]):
                self.I0 = self.inputdata["I0"] * np.pi / 180.
                self.e0 = self.inputdata["e0"]
                self.Omega0 = self.inputdata["Omega0"] * np.pi / 180.
                self.omega0 = self.inputdata["omega0"] * np.pi / 180.
            elif not (self.inputdata["EndState"]):
                if not (self.inputdata["SingleCollision"]):
                    self.I0 = nr.uniform(0, self.inputdata["I0"] * np.pi / 180., int(self.inputdata["Nparticles"]))
                    self.e0 = nr.uniform(0, self.inputdata["e0"], int(self.inputdata["Nparticles"]))
                    self.Omega0 = nr.uniform(0, 2 * np.pi, int(self.inputdata["Nparticles"]))
                    self.omega0 = nr.uniform(0, 2 * np.pi, int(self.inputdata["Nparticles"]))
                else:
                    self.I0 = nr.uniform(
                        np.max((0, self.inputdata["Icent"] * np.pi / 180. - self.inputdata["I0"] * np.pi / 180.)), \
                        self.inputdata["Icent"] * np.pi / 180. + self.inputdata["I0"] * np.pi / 180.,
                        int(self.inputdata["Nparticles"]))
                    self.e0 = nr.uniform(np.max((0, self.inputdata["ecent"] - self.inputdata["e0"])), \
                                         self.inputdata["ecent"] + self.inputdata["e0"],
                                         int(self.inputdata["Nparticles"]))
                    self.Omega0 = nr.uniform(self.inputdata["Omcent"] - self.inputdata["Om0"], \
                                             self.inputdata["Omcent"] + self.inputdata["Om0"],
                                             int(self.inputdata["Nparticles"]))
                    self.omega0 = nr.uniform(self.inputdata["omcent"] - self.inputdata["om0"], \
                                             self.inputdata["omcent"] + self.inputdata["om0"],
                                             int(self.inputdata["Nparticles"]))

        self.A = ff.A(self.a, self.Mp[0], self.ap[0], Mstar=self.Mstar)
        self.Aj = ff.Aj(self.a, self.Mp[0], self.ap[0], Mstar=self.Mstar)
        self.B = ff.B(self.a, self.Mp[0], self.ap[0], Mstar=self.Mstar)
        self.Bj = ff.Bj(self.a, self.Mp[0], self.ap[0], Mstar=self.Mstar)

        # forced component
        self.h0 = -(self.Aj / self.A) * self.ep[0] * np.sin(self.Omegap[0] + self.omegap[0])
        self.k0 = -(self.Aj / self.A) * self.ep[0] * np.cos(self.Omegap[0] + self.omegap[0])
        self.p0 = -(self.Bj / self.B) * self.Ip[0] * np.sin(self.Omegap[0])
        self.q0 = -(self.Bj / self.B) * self.Ip[0] * np.cos(self.Omegap[0])

        if not (self.inputdata["EndState"]):
            Ifreesin = self.I0 * np.sin(self.Omega0)
            Ifreecos = self.I0 * np.cos(self.Omega0)
            efreesin = self.e0 * np.sin(self.Omega0 + self.omega0)
            efreecos = self.e0 * np.cos(self.Omega0 + self.omega0)

            Ifree = np.sqrt(Ifreesin ** 2 + Ifreecos ** 2)
            gamma = np.arctan2(Ifreesin, Ifreecos)
            efree = np.sqrt(efreesin ** 2 + efreecos ** 2)
            beta = np.arctan2(efreesin, efreecos)

            self.h = efree * np.sin(self.A * self.age + beta) + self.h0
            self.k = efree * np.cos(self.A * self.age + beta) + self.k0
            self.p = Ifree * np.sin(self.B * self.age + gamma) + self.p0
            self.q = Ifree * np.cos(self.B * self.age + gamma) + self.q0
        else:
            efree = nr.uniform(0, self.inputdata["e0"], int(self.inputdata["Nparticles"]))
            Ifree = nr.uniform(
                np.max((0, self.inputdata["Icent"] * np.pi / 180. - self.inputdata["I0"] * np.pi / 180.)), \
                self.inputdata["I0"] * np.pi / 180. + self.inputdata["Icent"] * np.pi / 180.,
                int(self.inputdata["Nparticles"]))
            if self.inputdata["launchstyle"] >= 4: #if parent body orbits should be aligned
                omega = np.zeros(int(self.inputdata["Nparticles"]))
                Omega = np.zeros(int(self.inputdata["Nparticles"]))
            else:
                omega = nr.uniform(0, 2 * np.pi, int(self.inputdata["Nparticles"]))
                Omega = nr.uniform(0, 2 * np.pi, int(self.inputdata["Nparticles"]))
            self.h = efree * np.sin(omega + Omega) + self.h0
            self.k = efree * np.cos(omega + Omega) + self.k0
            self.p = Ifree * np.sin(Omega) + self.p0
            self.q = Ifree * np.cos(Omega) + self.q0

    def ComputeParentOrbital(self):
        # Compute orbital parameters of parent bodies
        self.Omega = np.arctan2(self.p, self.q)
        pomega = np.arctan2(self.h, self.k)
        self.omega = pomega - self.Omega
        self.I = np.sqrt(self.p ** 2 + self.q ** 2)
        self.e = np.sqrt(self.h ** 2 + self.k ** 2)

    def ComputeDustGrains(self, manual=False, beta=0.3, Nlaunch=10):
        # Compute orbital parameters of launched dust grains
        # beta = Prad/Pgrav
        # Nlaunch = launch points per parent body orbit
        print("Computing Dust Grain Orbits...")
        if manual:
            self.beta = beta
            Nlaunch = Nlaunch
        elif self.inputdata["betadistrb"] == 0:
            self.beta = self.inputdata["beta"]
            Nlaunch = int(self.inputdata["Nlaunch"])
            self.beta_dust = np.ones((len(self.h), Nlaunch)) * self.beta
        else:
            Nlaunch = int(self.inputdata["Nlaunch"])
            betapow = self.inputdata["betadistrb"]
            betamin, betamax = self.inputdata["betamin"], self.inputdata["betamax"]
            self.beta_dust = np.zeros((len(self.h), Nlaunch))

        self.beta_bounded = False
        if "beta_bounded" in self.inputdata and self.inputdata["beta_bounded"] == 1:
            self.beta_bounded = True
        self.a_dust = np.zeros((len(self.h), Nlaunch))  # len(self.h) is NParticles - this is Nparticles x Nlaunch array
        self.e_dust = np.zeros((len(self.h), Nlaunch))
        self.I_dust = np.zeros((len(self.h), Nlaunch))
        self.Omega_dust = np.zeros((len(self.h), Nlaunch))
        self.omega_dust = np.zeros((len(self.h), Nlaunch))

        lps = np.zeros((len(self.h), Nlaunch))
        mu = consts.G * self.Mstar

        for i in range(len(self.h)):  # for each parent body
            print("%i/%i parent body" % (i + 1, len(self.h)))
            if self.inputdata["launchstyle"] == 1:
                # uniform in f
                fp = nr.uniform(0, 2*np.pi, Nlaunch)
                lps[i] = fp
                cosfp = np.cos(fp)
                sinfp = np.sin(fp)
            elif self.inputdata["launchstyle"] == 2:
                # uniform in cosf
                cosfp = nr.uniform(-1, 1) * np.ones(Nlaunch)
                sinfp = nr.uniform(-1, 1) * np.ones(Nlaunch)
            elif self.inputdata["launchstyle"] == 3:
                # uniform in M
                f, fweight = ppo.OutputPosition(self.e[i], Npts=Nlaunch)  # ??
                lps[i] = f
                cosfp = np.cos(f)
                sinfp = np.sin(f)
            elif self.inputdata["launchstyle"] == 4:
                # all at peri
                cosfp = np.ones(Nlaunch)
                sinfp = np.zeros(Nlaunch)
            elif self.inputdata["launchstyle"] == 5:
                # all from quadrature
                cosfp = np.zeros(Nlaunch)
                sinfp = np.ones(Nlaunch)
            if self.inputdata["betadistrb"] != 0:
                for j in range(len(cosfp)):  # for each dust grain

                    if self.inputdata["sizedistrb"] == 1:
                        self.beta_dust[i, j] = bd.Donhanyi(self.e[i], cosfp[j], betapow=betapow, betamin=betamin,
                                                           betamax=betamax, Ndust=1, beta_bounded=self.beta_bounded,
                                                           a=self.a[i], Mstar=self.Mstar, Tage=self.age)
                    elif self.inputdata["sizedistrb"] == 2:
                        self.beta_dust[i, j] = bd.OrbTimeCorr(self.e[i], cosfp[j], betapow=betapow, betamin=betamin,
                                                              betamax=betamax, Ndust=1, beta_bounded=self.beta_bounded,
                                                              a=self.a[i], Mstar=self.Mstar, Tage=self.age)
                    elif self.inputdata["sizedistrb"] == 3:
                        self.beta_dust[i, j] = bd.OrbTimeCorrCirc(betapow=betapow, betamin=betamin, Ndust=1,
                                                                  beta_bounded=self.beta_bounded, a=self.a[i],
                                                                  Mstar=self.Mstar, Tage=self.age)
            if "ejected" in self.inputdata and self.inputdata["ejected"] == 1: #including ejecta velocity
                dv_ratio = self.inputdata["dv_ratio"]
                P1 = self.get_p1_matrix(w=self.omega[i])
                P2 = self.get_p2_matrix(I=self.I[i])
                P3 = self.get_p3_matrix(om=self.Omega[i])
                M = P3 @ P2 @ P1 # M matrix translates coordinates from orbital plane to equatorial system
                #numpy optimizations:
                # radius = self.a[i] * (1 - self.e[i] ** 2) / ( 1 + self.e[i] * cosfp)
                # velocity = np.sqrt(consts.G * self.Mstar * (2 / radius - 1 / self.a[i]))
                # drdt = self.a[i] * (1 - self.e[i] ** 2) * self.e[i] * sinfp / ((1 + self.e[i] * cosfp) ** 2)
                # def grSign(num):
                #   if num > 0:
                #       return 1
                #   return -1
                # temp_arr1 = np.where(cosfp > 0, velocity, -velocity)
                # choose_array = np.array([np.zeros(len(cosfp)), temp_arr1, np.zeros(len(cosfp))]).T
                # np.where(abs(-1 * radius * sinfp[j] + drdt * cosfp[j]) > 1e-40, np.zeros()
                # dydx = np.where(abs(-1 * radius * sinfp + drdt * cosfp) > 1e-40, (radius * cosfp - drdt * sinfp) / (-1 * radius * sinfp + drdt * cosfp), np.zeros(len(cosfp)))
                #
                # num = radius * cosfp - drdt * sinfp
                # div_zero = np.where(abs(-1 * radius * sinfp + drdt * cosfp) < 1e-40)
                # denom =
                #numpy optimizations end
                for j in range(len(self.a_dust[i])): #for each dust particle launched from this parent body
                    if self.verbose and j % self.print_every_x_dust == 0:
                        print("Computing jth dust grain...")
                    radius = self.a[i] * (1 - self.e[i] ** 2) / (1 + self.e[i] * cosfp[j])
                    velocity = np.sqrt(consts.G * self.Mstar * (2 / radius - 1 / self.a[i]))
                    drdt = self.a[i] * (1 - self.e[i] ** 2) * self.e[i] * sinfp[j] / ((1 + self.e[i] * cosfp[j]) ** 2)
                    if abs(-1 * radius * sinfp[j] + drdt * cosfp[j]) < 1e-40:
                        if cosfp[j] > 0:
                            velocity_orbplane = [0, velocity, 0]
                        else:
                            velocity_orbplane = [0, -velocity, 0]
                    else:
                        dydx = (radius * cosfp[j] - drdt * sinfp[j]) / (-1 * radius * sinfp[j] + drdt * cosfp[j])
                        velocity_orbplane = [1, dydx, 0]
                        velocity_orbplane = velocity/np.linalg.norm(velocity_orbplane) * velocity_orbplane
                    velocity_eq = M @ velocity_orbplane
                    coords_orbplane = [radius * cosfp[j], radius * sinfp[j], 0]
                    coords_eq = M @ coords_orbplane

                    a, e, I, O, w, f = self.get_orbital_elements_rand_dv(coords_eq, velocity_eq, dv_ratio, mu, 1e-40)
                    cosf = np.cos(f)
                    sinf = np.sin(f)
                    self.a_dust[i][j] = (1 - self.beta_dust[i][j]) * a * (1 - e ** 2) / (
                                1 - e ** 2 - 2 * self.beta_dust[i][j] * (1 + e * cosf))
                    self.e_dust[i][j] = np.sqrt(
                        e ** 2 + 2 * self.beta_dust[i][j] * e * cosf + self.beta_dust[i][j] ** 2) / (
                                                    1 - self.beta_dust[i][j])
                    self.omega_dust[i][j] = w + np.arctan2(self.beta_dust[i][j] * sinf,
                                                           e + self.beta_dust[i][j] * cosf)
                    self.I_dust[i][j] = I
                    self.Omega_dust[i][j] = O
                    cosfp[i] = cosf
                    sinfp[i] = sinf
            else:
                self.a_dust[i, :] = (1 - self.beta_dust[i, :]) * self.a[i] * (1 - self.e[i] ** 2) / \
                                    (1 - self.e[i] ** 2 - 2 * self.beta_dust[i, :] * (1 + self.e[i] * cosfp))
                self.e_dust[i, :] = np.sqrt(
                    self.e[i] ** 2 + 2 * self.beta_dust[i, :] * self.e[i] * cosfp + self.beta_dust[i, :] ** 2) / (
                                                1 - self.beta_dust[i, :])
                self.omega_dust[i,:] = self.omega[i]+np.arctan2(self.beta_dust[i,:]*sinfp, self.e[i]+self.beta_dust[i,:]*cosfp)

                self.I_dust[i, :] = self.I[i]
                self.Omega_dust[i, :] = self.Omega[i]
            uboundi = np.where(self.a_dust[i, :] < 0)[0] 
            # if len(uboundi) > 0 and self.inputdata["betadistrb"] != 0:
            #     pdb.set_trace()

        lps = lps.flatten()
        np.savetxt('launchpoints.txt', lps)

        self.a_dust = self.a_dust.flatten()
        self.e_dust = self.e_dust.flatten()
        self.I_dust = self.I_dust.flatten()
        self.Omega_dust = self.Omega_dust.flatten()
        self.omega_dust = self.omega_dust.flatten()
        self.beta_dust = self.beta_dust.flatten()

        uboundi = np.where(self.a_dust < 0)[0]
        if len(uboundi) == len(self.a_dust): pdb.set_trace()
        if len(uboundi) > 0:
            goodi = np.where(self.a_dust > 0)[0]
            self.a_dust = self.a_dust[goodi]
            self.e_dust = self.e_dust[goodi]
            self.I_dust = self.I_dust[goodi]
            self.Omega_dust = self.Omega_dust[goodi]
            self.omega_dust = self.omega_dust[goodi]
            self.beta_dust = self.beta_dust[goodi]

        self.SaveValues()
        print("Dust grains computed.")

    def ComputeBackgroundDustGrains(self, manual=False, beta=0.3, Nlaunchback=10):
        # Compute orbital parameters of launched dust grains
        # beta = Prad/Pgrav
        # Nback: number of parent background orbits
        # Nlaunchback = launch points per parent body orbit
        if "Nlaunchback" not in self.inputdata or "Nback" not in self.inputdata:
            return
        print("Computing Background Dust Grain Orbits...")
        if manual:
            self.beta = beta
            Nlaunchback = Nlaunchback
        elif self.inputdata["betadistrb"] == 0:
            self.beta = self.inputdata["beta"]
            Nlaunchback = int(self.inputdata["Nlaunchback"])
            self.beta_dust = np.ones((len(self.h), Nlaunchback)) * self.beta
        elif "beta_bounded" in self.inputdata and self.inputdata["beta_bounded"] == 1:
            Nlaunchback = int(self.inputdata["Nlaunchback"])
            betapow = self.inputdata["betadistrb"]
            betamin, betamax = self.inputdata["betamin"], self.inputdata["betamax"]
            self.beta_dust = np.zeros((len(self.h), Nlaunchback))
        else:
            Nlaunchback = int(self.inputdata["Nlaunchback"])
            betapow = self.inputdata["betadistrb"]
            betamin, betamax = self.inputdata["betamin"], self.inputdata["betamax"]
            self.beta_dust = np.zeros((len(self.h), Nlaunchback))

        self.a_dust = np.zeros((len(self.h), Nlaunchback))  # len(self.h) is Nback - this is Nback x Nlaunchback array
        self.e_dust = np.zeros((len(self.h), Nlaunchback))
        self.I_dust = np.zeros((len(self.h), Nlaunchback))
        self.Omega_dust = np.zeros((len(self.h), Nlaunchback))
        self.omega_dust = np.zeros((len(self.h), Nlaunchback))

        lps = np.zeros((len(self.h), Nlaunchback))

        for i in range(len(self.h)):  # for each parent body
            print("%i/%i background parent body" % (i + 1, len(self.h)))
            fp = nr.uniform(0, 2 * np.pi, Nlaunchback)  # get random true anomaly for each
            lps[i] = fp
            cosfp = np.cos(fp)
            sinfp = np.sin(fp)
            if self.inputdata["betadistrb"] != 0:
                for j in range(len(cosfp)):  # for each dust grain
                    if self.inputdata["sizedistrb"] == 1:
                        self.beta_dust[i, j] = bd.Donhanyi(self.e[i], cosfp[j], betapow=betapow, betamin=betamin,
                                                           betamax=betamax, Ndust=1, beta_bounded=self.beta_bounded,
                                                           a=self.a[i], Mstar=self.Mstar, Tage=self.age)
                    elif self.inputdata["sizedistrb"] == 2:
                        self.beta_dust[i, j] = bd.OrbTimeCorr(self.e[i], cosfp[j], betapow=betapow, betamin=betamin,
                                                              betamax=betamax, Ndust=1, beta_bounded=self.beta_bounded,
                                                              a=self.a[i], Mstar=self.Mstar, Tage=self.age)
                    elif self.inputdata["sizedistrb"] == 3:
                        self.beta_dust[i, j] = bd.OrbTimeCorrCirc(betapow=betapow, betamin=betamin, Ndust=1,
                                                                  beta_bounded=self.beta_bounded, a=self.a[i],
                                                                  Mstar=self.Mstar, Tage=self.age)
            self.a_dust[i, :] = (1 - self.beta_dust[i, :]) * self.a[i] * (1 - self.e[i] ** 2) / \
                                (1 - self.e[i] ** 2 - 2 * self.beta_dust[i, :] * (1 + self.e[i] * cosfp))
            self.e_dust[i, :] = np.sqrt(
                self.e[i] ** 2 + 2 * self.beta_dust[i, :] * self.e[i] * cosfp + self.beta_dust[i, :] ** 2) / (
                                            1 - self.beta_dust[i, :])
            self.omega_dust[i, :] = self.omega[i] + np.arctan2(self.beta_dust[i, :] * sinfp,
                                                               self.e[i] + self.beta_dust[i, :] * cosfp)

            self.I_dust[i, :] = self.I[i]
            self.Omega_dust[i, :] = self.Omega[i]

            uboundi = np.where(self.a_dust[i, :] < 0)[0]
            if len(uboundi) > 0 and self.inputdata["betadistrb"] != 0:
                pdb.set_trace()
        lps = lps.flatten()
        np.savetxt('backgroundlaunchpoints.txt', lps)

        self.a_dust = self.a_dust.flatten()
        self.e_dust = self.e_dust.flatten()
        self.I_dust = self.I_dust.flatten()
        self.Omega_dust = self.Omega_dust.flatten()
        self.omega_dust = self.omega_dust.flatten()
        self.beta_dust = self.beta_dust.flatten()

        uboundi = np.where(self.a_dust < 0)[0]
        if len(uboundi) == len(self.a_dust): pdb.set_trace()
        if len(uboundi) > 0:
            goodi = np.where(self.a_dust > 0)[0]
            self.a_dust = self.a_dust[goodi]
            self.e_dust = self.e_dust[goodi]
            self.I_dust = self.I_dust[goodi]
            self.Omega_dust = self.Omega_dust[goodi]
            self.omega_dust = self.omega_dust[goodi]
            self.beta_dust = self.beta_dust[goodi]
        self.SaveValues()

    # Save parameters. Useful for persistence between primary and background disk
    def SaveValues(self):
        self.a_all.extend(self.a_dust)
        #print(self.a_all)
        self.e_all.extend(self.e_dust)
        self.I_all.extend(self.I_dust)
        self.Omega_all.extend(self.Omega_dust)
        self.omega_all.extend(self.omega_dust)
        self.beta_all.extend(self.beta_dust)


    ##Output ASCII file
    def OutputParentOrbit(self, outfile):
        np.savetxt(outfile + "_pqhk.txt", list(zip(self.p, self.q, self.h, self.k, self.p0, \
                                                   self.q0, self.h0, self.k0, self.a)))
        np.savetxt(outfile + "_orbit.txt", list(zip(self.a, self.e, self.I, self.Omega, self.omega)))
        np.savetxt(outfile + "_freq.txt", list(zip(self.A, self.Aj, self.B, self.Bj)))

    def OutputDustOrbit(self, outfile):
        if self.inputdata["betadistrb"] == 0:
            np.savetxt(outfile + "_dustorbit.txt", list(zip(self.a_dust, self.e_dust, \
                                                            self.I_dust, self.Omega_dust, \
                                                            self.omega_dust)))
        else:
            np.savetxt(outfile + "_dustorbit.txt", list(zip(self.a_dust, self.e_dust, \
                                                            self.I_dust, self.Omega_dust, \
                                                            self.omega_dust, self.beta_dust)))

    def OutputDustAndBackOrbit(self, outfile):
        if self.inputdata["betadistrb"] == 0:
            np.savetxt(outfile + "_dustorbit.txt", list(zip(self.a_all, self.e_all, \
                                                            self.I_all, self.Omega_all, \
                                                            self.omega_all)))
        else:
            np.savetxt(outfile + "_dustorbit.txt", list(zip(self.a_all, self.e_all, \
                                                            self.I_all, self.Omega_all, \
                                                            self.omega_all, self.beta_all)))
    ##Plot
    def Plot_pqhk(self, pqxlim=(-2, 2), pqylim=(-2, 2)):
        pqxlim = (-np.max(np.sqrt(self.q ** 2 + self.p ** 2)), np.max(np.sqrt(self.q ** 2 + self.p ** 2)))
        pl.figure(figsize=(10, 5))
        pl.subplot(121, aspect='equal')
        pl.plot(self.q / np.pi, self.p / np.pi, 'ko')
        pl.plot(self.q0 / np.pi, self.p0 / np.pi, 'rx', markersize=14, markeredgewidth=2)
        pl.plot(self.Ip * np.cos(self.Omegap) / np.pi, self.Ip * np.sin(self.Omegap) / np.pi, 'bo', markersize=12, \
                markeredgewidth=4, markerfacecolor='none', markeredgecolor='blue')
        pl.plot((-7, 7), (0, 0), 'k--')
        pl.plot((0, 0), (-7, 7), 'k--')
        pl.xlim(pqxlim)
        pl.ylim(pqxlim)
        pl.xlabel(r'$I\cos\Omega/\pi$')
        pl.ylabel(r'$I\sin\Omega/\pi$')

        pl.subplot(122, aspect='equal')
        pl.plot(self.k, self.h, 'ko')
        pl.plot(self.k0, self.h0, 'rx', markersize=14, markeredgewidth=2)
        pl.plot(self.ep * np.cos(self.Omegap + self.omegap), \
                self.ep * np.sin(self.Omegap + self.omegap), 'bo', \
                markersize=12, markeredgewidth=4, markerfacecolor='none', markeredgecolor='blue')
        pl.plot(0.5 * np.cos(np.linspace(0, 2 * np.pi, 100)), 0.5 * np.sin(np.linspace(0, 2 * np.pi, 100)), 'r--')
        pl.plot((-1, 1), (0, 0), 'k--')
        pl.plot((0, 0), (-1, 1), 'k--')
        pl.xlabel(r'$e\cos\varpi$')
        pl.ylabel(r'$e\sin\varpi$')
        pl.ylim((-1, 1))
        pl.xlim((-1, 1))
        pl.tight_layout()
        pl.show()

    def Plot3DOrbit(self, elev=0, azim=0, dim=(-100, 100), single_planet=True):
        fig = pl.figure(2, figsize=(8, 8))  # initialize figure
        ax = Axes3D(fig)
        ax.plot(np.zeros(1), np.zeros(1), np.zeros(1), 'y*', markersize=20)  # draw the central star
        for i in range(len(self.h)):
            fd.draw_wire(ax, self.Omega[i], self.I[i], argu_peri=self.omega[i], \
                         e=self.e[i], r0=self.a[i], format='b-', lw=0.1)  # draw parent body orbits
        # for i in range(len(self.a_dust.flatten())):
        #    fd.draw_particles(ax, self.Omega_dust.flatten()[i], self.I_dust.flatten()[i], \
        #                   argu_peri=self.omega_dust.flatten()[i], e=self.e_dust.flatten()[i], \
        #                      r0=self.a_dust.flatten()[i], format='ro', msz=0.5) #draw dust particles
        fd.draw_wire(ax, self.Omegap[0], self.Ip[0], argu_peri=self.omegap[0], \
                     e=self.ep[0], r0=self.ap[0])  # draw planet orbit
        if not (single_planet):
            fd.draw_wire(ax, self.Omegap[1], self.Ip[1], argu_peri=self.omegap[1], \
                         e=self.ep[1], r0=self.ap[1])  # draw planet orbit
        ax.view_init(elev=elev, azim=azim)
        ax.set_ylim3d(dim)
        ax.set_zlim3d(dim)
        ax.set_xlim3d(dim)
        ax.set_position([0., 0., 1., 1.])
        fig.show()


    def ComputeBackgroundParentSingle(self, manual=False, Nback=500, amin=4., amax=6., I0=5. * (np.pi / 180.), e0=0.2,
                            Omega0=0., omega0=0., random=False):
        # Compute p, q, h, k of parent bodies for background disk
        # Single planet only
        # Nback: number of parent bodies
        # amin, amax: min, max semi-major axes of parent bodies
        # I0, e0: initial mutual inclination and eccentricities
        # Omega0, omega0: initial nodal angle and argument of periapse
        if "Nlaunchback" not in self.inputdata or "Nback" not in self.inputdata:
            return

        print("Computing Background Parent Orbits (single planet)...")
        if self.freeelem:
            self.a, self.I0, self.e0, self.Omega0, self.omega0 = np.loadtxt(self.freeelemtxt, unpack=True)
        elif manual:
            self.a = np.linspace(amin, amax, int(self.inputdata["Nback"]))
            self.I0 = I0
            self.e0 = e0
            self.Omega0 = Omega0
            self.omega0 = omega0
        else:
            if self.inputdata["InnerPlanet"] == 1:
                amin = self.ap * (1 + 2.0 * (self.Mp * consts.Mearth / self.Mstar) ** (2. / 7.))
                if self.inputdata["OuterPlanet"] == 0:
                    # amax = amin/(1-self.inputdata["awidth"])
                    amax = amin * (1. + self.inputdata["awidth"])
            if self.inputdata["OuterPlanet"] == 1:
                amax = self.ap * (1 - 2.0 * (self.Mp * consts.Mearth / self.Mstar) ** (2. / 7.))
                if self.inputdata["InnerPlanet"] == 0:
                    amin = amax * (1 - self.inputdata["awidth"])

            self.a = np.linspace(amin, amax, int(self.inputdata["Nback"]))

            if not (self.inputdata["Random"]):
                self.I0 = self.inputdata["I0"] * np.pi / 180.
                self.e0 = self.inputdata["e0"]
                self.Omega0 = self.inputdata["Omega0"] * np.pi / 180.
                self.omega0 = self.inputdata["omega0"] * np.pi / 180.
            elif not (self.inputdata["EndState"]):
                if not (self.inputdata["SingleCollision"]):
                    self.I0 = nr.uniform(0, self.inputdata["I0"] * np.pi / 180., int(self.inputdata["Nback"]))
                    self.e0 = nr.uniform(0, self.inputdata["e0"], int(self.inputdata["Nback"]))
                    self.Omega0 = nr.uniform(0, 2 * np.pi, int(self.inputdata["Nback"]))
                    self.omega0 = nr.uniform(0, 2 * np.pi, int(self.inputdata["Nback"]))
                else:
                    self.I0 = nr.uniform(
                        np.max((0, self.inputdata["Icent"] * np.pi / 180. - self.inputdata["I0"] * np.pi / 180.)), \
                        self.inputdata["Icent"] * np.pi / 180. + self.inputdata["I0"] * np.pi / 180.,
                        int(self.inputdata["Nback"]))
                    self.e0 = nr.uniform(np.max((0, self.inputdata["ecent"] - self.inputdata["e0"])), \
                                         self.inputdata["ecent"] + self.inputdata["e0"],
                                         int(self.inputdata["Nback"]))
                    self.Omega0 = nr.uniform(self.inputdata["Omcent"] - self.inputdata["Om0"], \
                                             self.inputdata["Omcent"] + self.inputdata["Om0"],
                                             int(self.inputdata["Nback"]))
                    self.omega0 = nr.uniform(self.inputdata["omcent"] - self.inputdata["om0"], \
                                             self.inputdata["omcent"] + self.inputdata["om0"],
                                             int(self.inputdata["Nback"]))
        self.A = ff.A(self.a, self.Mp[0], self.ap[0], Mstar=self.Mstar)
        self.Aj = ff.Aj(self.a, self.Mp[0], self.ap[0], Mstar=self.Mstar)
        self.B = ff.B(self.a, self.Mp[0], self.ap[0], Mstar=self.Mstar)
        self.Bj = ff.Bj(self.a, self.Mp[0], self.ap[0], Mstar=self.Mstar)

        # forced component
        self.h0 = -(self.Aj / self.A) * self.ep[0] * np.sin(self.Omegap[0] + self.omegap[0])
        self.k0 = -(self.Aj / self.A) * self.ep[0] * np.cos(self.Omegap[0] + self.omegap[0])
        self.p0 = -(self.Bj / self.B) * self.Ip[0] * np.sin(self.Omegap[0])
        self.q0 = -(self.Bj / self.B) * self.Ip[0] * np.cos(self.Omegap[0])

        if not (self.inputdata["EndState"]):
            Ifreesin = self.I0 * np.sin(self.Omega0)
            Ifreecos = self.I0 * np.cos(self.Omega0)
            efreesin = self.e0 * np.sin(self.Omega0 + self.omega0)
            efreecos = self.e0 * np.cos(self.Omega0 + self.omega0)

            Ifree = np.sqrt(Ifreesin ** 2 + Ifreecos ** 2)
            gamma = np.arctan2(Ifreesin, Ifreecos)
            efree = np.sqrt(efreesin ** 2 + efreecos ** 2)
            beta = np.arctan2(efreesin, efreecos)

            self.h = efree * np.sin(self.A * self.age + beta) + self.h0
            self.k = efree * np.cos(self.A * self.age + beta) + self.k0
            self.p = Ifree * np.sin(self.B * self.age + gamma) + self.p0
            self.q = Ifree * np.cos(self.B * self.age + gamma) + self.q0
        else:
            efree = nr.uniform(0, self.inputdata["e0"], int(self.inputdata["Nback"]))
            Ifree = nr.uniform(
                np.max((0, self.inputdata["Icent"] * np.pi / 180. - self.inputdata["I0"] * np.pi / 180.)), \
                self.inputdata["I0"] * np.pi / 180. + self.inputdata["Icent"] * np.pi / 180.,
                int(self.inputdata["Nback"]))
            omega = nr.uniform(0, 2 * np.pi, int(self.inputdata["Nback"]))
            Omega = nr.uniform(0, 2 * np.pi, int(self.inputdata["Nback"]))
            self.h = efree * np.sin(omega + Omega) + self.h0
            self.k = efree * np.cos(omega + Omega) + self.k0
            self.p = Ifree * np.sin(Omega) + self.p0
            self.q = Ifree * np.cos(Omega) + self.q0


    def ComputeBackgroundParentOrbital(self):
        # Compute orbital parameters of parent bodies
        self.Omega = np.arctan2(self.p, self.q)
        pomega = np.arctan2(self.h, self.k)
        self.omega = pomega - self.Omega
        self.I = np.sqrt(self.p ** 2 + self.q ** 2)
        self.e = np.sqrt(self.h ** 2 + self.k ** 2)


    def get_orbital_elements(self, r0, v, dv, mu, eps):
        # returns the orbital elements (a, e, i, O, w, f) of an orbit
        # calculated based on  r0, v0 (= v + dv), mu (= GM)
        # eps controls for numerical precision
        def get_norm(vector):
            return np.linalg.norm(vector)

        i_dir = np.array([1, 0, 0])
        j_dir = np.array([0, 1, 0])
        k_dir = np.array([0, 0, 1])
        v0 = v + dv
        speed = get_norm(v0)
        radius = get_norm(r0)

        h_vec = np.cross(r0, v0)
        h = get_norm(h_vec)
        h_z = h_vec @ k_dir
        node_vec = np.cross(k_dir, h_vec)
        n = get_norm(node_vec)
        n_x = node_vec @ i_dir

        scaled_r_vec = (speed ** 2 - mu / radius) * r0
        scaled_v_vec = (r0 @ v0) * v0
        e_vec = (scaled_r_vec - scaled_v_vec) / mu
        ecc = get_norm(e_vec)
        Energy = (speed ** 2) / 2 - mu / radius

        if abs(ecc - 1) < eps:
            a = float('inf')
            peri = h ** 2 / mu
        else:
            a = - mu / (2 * Energy)
            peri = a * (1 - ecc ** 2)

        I = np.arccos(h_z / h)

        if abs(ecc) < eps:  # circular
            if abs(I) < eps or abs(I - np.pi) < eps:  # equatorial
                f = np.arccos(r0[0] / radius)
                if v0[0] > 0:
                    f = 2 * np.pi - f
                w = 0  # "periapsis" at launch
                Omega = 0
            else:
                f = np.arccos((node_vec @ r0) / (n * radius))
                if r0[2] < 0:
                    f = 2 * np.pi - f
                w = f
                Omega = np.arccos(n_x / n)
        else:
            if abs(I) < eps or abs(I - np.pi) < eps:
                Omega = 0  # convention
                w = np.arctan2(e_vec[1], e_vec[0])
                if np.cross(r0, v0)[2] < 0:
                    w = 2 * np.pi - w
            else:
                Omega = np.arccos(n_x / n)
                w = np.arccos((node_vec @ e_vec) / (n * ecc))
            f = np.arccos((e_vec @ r0) / (ecc * radius))

            if node_vec[1] < 0:
                Omega = 2 * np.pi - Omega

            if r0 @ v0 < 0:
                f = 2 * np.pi - f
        return (a, ecc, I, Omega, w, f)


    def get_orbital_elements_opt(self, r0, v, dv, mu, eps):
        # returns the orbital elements (a, e, i, O, w, f) of an orbit
        # calculated based on  r0, v0 (= v + dv), mu (= GM)
        # eps controls for numerical precision
        def get_norm(vector):
            return np.linalg.norm(vector)

        k_dir = np.array([0, 0, 1])
        v0 = v + dv
        speed = get_norm(v0)
        radius = get_norm(r0)

        h_vec = np.cross(r0, v0)
        h = get_norm(h_vec)
        h_z = h_vec[2]
        node_vec = np.cross(k_dir, h_vec)
        n = get_norm(node_vec)
        n_x = node_vec[0]

        scaled_r_vec = (speed ** 2 - mu / radius) * r0
        scaled_v_vec = (r0 @ v0) * v0
        e_vec = (scaled_r_vec - scaled_v_vec) / mu
        ecc = get_norm(e_vec)
        Energy = (speed ** 2) / 2 - mu / radius

        if abs(ecc - 1) < eps:
            a = float('inf')
        else:
            a = - mu / (2 * Energy)

        I = np.arccos(h_z / h)

        if abs(ecc) < eps:  # circular
            if abs(I) < eps or abs(I - np.pi) < eps:  # equatorial
                f = np.arccos(r0[0] / radius)
                if v0[0] > 0:
                    f = 2 * np.pi - f
                w = 0  # "periapsis" at launch
                Omega = 0
            else:
                f = np.arccos((node_vec @ r0) / (n * radius))
                if r0[2] < 0:
                    f = 2 * np.pi - f
                w = f
                Omega = np.arccos(n_x / n)
        else:
            if abs(I) < eps or abs(I - np.pi) < eps:
                Omega = 0  # convention
                w = np.arctan2(e_vec[1], e_vec[0])
                if np.cross(r0, v0)[2] < 0:
                    w = 2 * np.pi - w
            else:
                Omega = np.arccos(n_x / n)
                w = np.arccos((node_vec @ e_vec) / (n * ecc))
            f = np.arccos((e_vec @ r0) / (ecc * radius))

            if node_vec[1] < 0:
                Omega = 2 * np.pi - Omega

            if r0 @ v0 < 0:
                f = 2 * np.pi - f
        return (a, ecc, I, Omega, w, f)

    # Returns orbital elements of an orbit with random dv of a certain RATIO
    def get_orbital_elements_rand_dv(self, r0, v, ratio, mu, eps, optimized=True):
        i_dir = np.array([1, 0, 0])
        j_dir = np.array([0, 1, 0])
        k_dir = np.array([0, 0, 1])
        dv_mag = ratio * np.linalg.norm(v)
        dv_x = np.random.normal()
        dv_y = np.random.normal()
        dv_z = np.random.normal()

        norm = math.pow(dv_x ** 2 + dv_y ** 2 + dv_z ** 2, 0.5)
        dv = dv_mag * (dv_x * i_dir + dv_y * j_dir + dv_z * k_dir) / norm
        if optimized:
            return self.get_orbital_elements_opt(r0, v, dv, mu, eps)
        else:
            return self.get_orbital_elements(r0, v, dv, mu, eps)

    # returns the P1 rotation matrix. See Murray and Dermott p.50 for details
    def get_p1_matrix(self, w):
        P1 = np.zeros((3, 3))
        P1[0][0] = np.cos(w)
        P1[0][1] = -1 * np.sin(w)
        P1[1][0] = np.sin(w)
        P1[1][1] = np.cos(w)
        P1[2][2] = 1
        return P1

    def get_p2_matrix(self, I):
        P2 = np.zeros((3, 3))
        P2[0][0] = 1
        P2[1][1] = np.cos(I)
        P2[1][2] = -1 * np.sin(I)
        P2[2][1] = np.sin(I)
        P2[2][2] = np.cos(I)
        return P2

    def get_p3_matrix(self, om):
        P3 = np.zeros((3, 3))
        P3[0][0] = np.cos(om)
        P3[0][1] = -1 * np.sin(om)
        P3[1][0] = np.sin(om)
        P3[1][1] = np.cos(om)
        P3[2][2] = 1
        return P3