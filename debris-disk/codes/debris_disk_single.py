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
import datetime
import time
import random
import scipy.integrate as sint
import scipy.interpolate as si


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
        self.beta_per_launch = 1
        self.beta_per_launch_back = 1
        self.stabfac = 0.997
        self.stabfac_back = 0.997
        self.num_e = 10
        self.Nlaunchback = 0
        if "verbose" in self.inputdata:
            self.verbose = self.inputdata["verbose"]
            if "print_every_x_dust" in self.inputdata:
                self.print_every_x_dust = self.inputdata["print_every_x_dust"]
        if "beta_per_launch" in self.inputdata:
            self.beta_per_launch = int(self.inputdata["beta_per_launch"])
        if "beta_per_launch_back" in self.inputdata:
            self.beta_per_launch_back = int(self.inputdata["beta_per_launch_back"])
        if "stabfac" in self.inputdata:
            self.stabfac = self.inputdata["stabfac"]
        if "stabfac_back" in self.inputdata:
            self.stabfac_back = self.inputdata["stabfac_back"]
        if "num_e" in self.inputdata:
            self.num_e = self.inputdata["num_e"]
        if "Nlaunchback" in self.inputdata:
            self.Nlaunchback = self.inputdata["Nlaunchback"]


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


    def ComputeParentSingle(self, manual=False, amin=4., amax=6., Icoll=5. * (np.pi / 180.), ecoll=0.2,
                            Omega0=0., omega0=0., fix_a=False, fixed_a=50, coll_in_middle=False):
        # Compute p, q, h, k of parent body
        # Single planet only
        # amin, amax: min, max semi-major axes of parent bodies
        # Icoll, ecoll: initial mutual inclination (in radians) and eccentricities
        # Omega0, omega0: initial nodal angle and argument of periapse
        # fix_a: flag for fixing all of the a-values of each particle (should only be one) to fixed_a
        # coll_in_middle: flag that determines whether the a value of the parent is constrained to be in the middle of amin, amax
        print("Computing Parent Orbits (single planet)...")
        if self.freeelem:
            self.a, self.Icoll, self.ecoll, self.Omega0, self.omega0 = np.loadtxt(self.freeelemtxt, unpack=True)
        elif manual:
            self.a = np.linspace(amin, amax, int(self.inputdata["Nparticles"]))
            self.Icoll = Icoll
            self.ecoll = ecoll
            self.Omega0 = Omega0
            self.omega0 = omega0
        else:
            if fix_a:
                self.a = fixed_a * np.ones(int(self.inputdata["Nparticles"]))
            else:
                if self.inputdata["InnerPlanet"] == 1:
                    amin = self.ap * (1 + 2.0 * (self.Mp * consts.Mearth / self.Mstar) ** (2. / 7.))
                    if self.inputdata["OuterPlanet"] == 0:
                        amax = amin * (1. + self.inputdata["awidth"])
                if self.inputdata["OuterPlanet"] == 1:
                    amax = self.ap * (1 - 2.0 * (self.Mp * consts.Mearth / self.Mstar) ** (2. / 7.))
                    if self.inputdata["InnerPlanet"] == 0:
                        amin = amax * (1 - self.inputdata["awidth"])
                if coll_in_middle:
                    self.a = (amin + amax) / 2 * np.ones(int(self.inputdata["Nparticles"]))
                else:
                    self.a = np.linspace(amin, amax, int(self.inputdata["Nparticles"]))

            if not (self.inputdata["Random"]):
                self.Icoll = self.inputdata["Icoll"] * np.pi / 180.
                self.ecoll = self.inputdata["ecoll"]
                self.Omega0 = self.inputdata["Omega0"] * np.pi / 180.
                self.omega0 = self.inputdata["omega0"] * np.pi / 180.
            elif not (self.inputdata["EndState"]):
                if not (self.inputdata["SingleCollision"]):
                    self.Icoll = nr.uniform(0, self.inputdata["Icoll"] * np.pi / 180., int(self.inputdata["Nparticles"]))
                    self.ecoll = nr.uniform(0, self.inputdata["ecoll"], int(self.inputdata["Nparticles"]))
                    self.Omega0 = nr.uniform(0, 2 * np.pi, int(self.inputdata["Nparticles"]))
                    self.omega0 = nr.uniform(0, 2 * np.pi, int(self.inputdata["Nparticles"]))
                else:
                    self.Icoll = nr.uniform(
                        np.max((0, self.inputdata["Icent"] * np.pi / 180. - self.inputdata["Icoll"] * np.pi / 180.)), \
                        self.inputdata["Icent"] * np.pi / 180. + self.inputdata["Icoll"] * np.pi / 180.,
                        int(self.inputdata["Nparticles"]))
                    self.ecoll = nr.uniform(np.max((0, self.inputdata["ecent"] - self.inputdata["ecoll"])), \
                                         self.inputdata["ecent"] + self.inputdata["ecoll"],
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
            Ifreesin = self.Icoll * np.sin(self.Omega0)
            Ifreecos = self.Icoll * np.cos(self.Omega0)
            efreesin = self.ecoll * np.sin(self.Omega0 + self.omega0)
            efreecos = self.ecoll * np.cos(self.Omega0 + self.omega0)

            Ifree = np.sqrt(Ifreesin ** 2 + Ifreecos ** 2)
            gamma = np.arctan2(Ifreesin, Ifreecos)
            efree = np.sqrt(efreesin ** 2 + efreecos ** 2)
            beta = np.arctan2(efreesin, efreecos)

            self.h = efree * np.sin(self.A * self.age + beta) + self.h0
            self.k = efree * np.cos(self.A * self.age + beta) + self.k0
            self.p = Ifree * np.sin(self.B * self.age + gamma) + self.p0
            self.q = Ifree * np.cos(self.B * self.age + gamma) + self.q0
        else:
            efree = nr.uniform(0, self.inputdata["ecoll"], int(self.inputdata["Nparticles"]))
            if "fixed_inc" in self.inputdata and self.inputdata["fixed_inc"] == 1:
                Ifree = self.inputdata["Icoll"] * np.pi / 180
            else:
                Ifree = nr.uniform(
                    np.max((0, self.inputdata["Icent"] * np.pi / 180. - self.inputdata["Icoll"] * np.pi / 180.)), \
                    self.inputdata["Icoll"] * np.pi / 180. + self.inputdata["Icent"] * np.pi / 180.,
                    int(self.inputdata["Nparticles"]))
            if self.inputdata["launchstyle"] == 4:
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
        # Compute orbital parameters of parent body
        self.e = np.sqrt(self.h ** 2 + self.k ** 2)
        self.I = np.sqrt(self.p ** 2 + self.q ** 2)
        if "acoll" in self.inputdata:
            self.a = [float(self.inputdata["acoll"])]
        if self.inputdata["launchstyle"] == 7: # custom input
            R = self.inputdata["radius"]
            e = np.array([self.inputdata["ecoll"]])
            self.e = e
            f = self.inputdata["launch_angle"] * np.pi / 180
            self.a = R*(1.+e*np.cos(f))/(1.-e**2)
            self.Omega = np.array([0.0])
            if "Omcoll" in self.inputdata:
                self.Omega = np.array([self.inputdata["Omcoll"]])
            self.omega = np.array([0.0]) - self.Omega - f

        elif self.inputdata["launchstyle"] == 8: # quadrature 1
            self.a = np.array([float(self.inputdata["acoll"])])
            self.e = np.array([self.inputdata["ecoll"]])
            self.Omega = np.array([0.0])
            f = 0 # could be pi/2 or 0
            self.Omega = np.array([0.0])
            self.omega = np.array([0.0]) - self.Omega - f
        else:
            self.Omega = np.arctan2(self.p, self.q)
            pomega = np.arctan2(self.h, self.k)
            self.omega = pomega - self.Omega

        if self.verbose:
            print("e = ", self.e)
            print("I = ", self.I)
            print("a = ", self.a)
            print("Omega = ", self.Omega)
            print("omega = ", self.omega)
            print("------------------------")


    def ComputeDustGrains_Optimized(self, manual=False, beta=0.3, Nlaunch=10):
        # Compute orbital parameters of launched dust grains
        # beta = Prad/Pgrav
        # Nlaunch = launch points per parent body orbit
        # print("Computing Dust Grain Orbits...")
        # ct = datetime.datetime.now()
        # print("current time:-", ct)
        start_time = time.time()
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
        if ("beta_bounded" in self.inputdata and self.inputdata["beta_bounded"] == 1) or ("beta_change" in self.inputdata and self.inputdata["beta_change"] == 1):
            self.beta_bounded = True
        self.a_dust = np.zeros((len(self.h), Nlaunch))  # len(self.h) is NParticles - this is Nparticles x Nlaunch array
        self.e_dust = np.zeros((len(self.h), Nlaunch))
        self.I_dust = np.zeros((len(self.h), Nlaunch))
        self.Omega_dust = np.zeros((len(self.h), Nlaunch))
        self.omega_dust = np.zeros((len(self.h), Nlaunch))

        lps = np.zeros((len(self.h), Nlaunch))
        mu = consts.G * self.Mstar
        matrix_time = 0
        ejecta_time = 0
        beta_calcs_time = 0
        beta_application_time = 0
        time_in_func = 0
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
                sinfp = -1 * np.ones(Nlaunch)
            elif self.inputdata["launchstyle"] == 6:
                # all at apo
                cosfp = -1 * np.ones(Nlaunch)
                sinfp = np.zeros(Nlaunch)
            elif self.inputdata["launchstyle"] == 7:
                # variable launch angle
                theta = self.inputdata["launch_angle"] * np.pi/180
                cosfp = np.cos(theta) * np.ones(Nlaunch)
                sinfp = np.sin(theta) * np.ones(Nlaunch)
            elif self.inputdata["launchstyle"] == 8:
                # all from quadrature
                cosfp = np.zeros(Nlaunch)
                sinfp = np.ones(Nlaunch)
            elif self.inputdata["launchstyle"] == 9:
                # all from other quadrature
                cosfp = np.zeros(Nlaunch)
                sinfp = -1 * np.ones(Nlaunch)

            #initial: orbital elements of dust grains before applying radiation pressure
            self.a_initial = np.array(Nlaunch * [self.a[i]])
            self.e_initial = np.array(Nlaunch * [self.e[i]])
            self.I_initial = np.array(Nlaunch * [self.I[i]])
            self.Omega_initial = np.array(Nlaunch * [self.Omega[i]])
            self.omega_initial = np.array(Nlaunch * [self.omega[i]])
            self.cosf_initial = cosfp
            self.sinf_initial = sinfp
            self.e_max = 0
            self.cosf_max = 0
            start_ejecta = time.time()

            #If ejected: calculate initial orbits after ejection, before radiation
            if "ejected" in self.inputdata and self.inputdata["ejected"] == 1: #including ejecta velocity
                start_matrices = time.time()
                dv_ratio = self.inputdata["dv_ratio"]
                P1 = self.get_p1_matrix(w=self.omega[i])
                P2 = self.get_p2_matrix(I=self.I[i])
                P3 = self.get_p3_matrix(om=self.Omega[i])
                M = P3 @ P2 @ P1 # M matrix translates coordinates from orbital plane to equatorial system

                end_matrices = time.time()
                matrix_time += end_matrices - start_matrices
                for j in range(len(self.a_dust[i])): #for each dust particle launched from this parent body, l
                    if self.verbose and j % self.print_every_x_dust == 0:
                        print(f"Computing {j}th dust grain...")

                    radius = self.a[i] * (1 - self.e[i] ** 2) / (1 + self.e[i] * cosfp[j])
                    #print(radius)


                    velocity = np.sqrt(consts.G * self.Mstar * (2 / radius - 1 / self.a[i]))


                    # drdt = self.a[i] * (1 - self.e[i] ** 2) * self.e[i] * sinfp[j] / ((1 + self.e[i] * cosfp[j]) ** 2)
                    # if abs(-1 * radius * sinfp[j] + drdt * cosfp[j]) < 1e-40:
                    #     if cosfp[j] > 0:
                    #         velocity_orbplane = [0, velocity, 0]
                    #     else:
                    #         velocity_orbplane = [0, -velocity, 0]
                    # else:
                    #     dydx = (radius * cosfp[j] - drdt * sinfp[j]) / (-1 * radius * sinfp[j] + drdt * cosfp[j])
                    #     velocity_orbplane = [1, dydx, 0]
                    #     ratio_ = velocity / np.linalg.norm(velocity_orbplane)
                    #     velocity_orbplane = [v * ratio_ for v in velocity_orbplane]
                    #
                    # old_orbplane = np.array(list(velocity_orbplane))
                    if abs(sinfp[j]) < 1e-40:
                        velocity_orbplane = [0, cosfp[j] * velocity, 0]
                    else:
                        dydx = - (cosfp[j] + self.e[i])/(sinfp[j])
                        velocity_orbplane = [1, dydx, 0]
                        ratio_ = velocity / np.linalg.norm(velocity_orbplane)
                        velocity_orbplane = [v * ratio_ for v in velocity_orbplane]
                        if cosfp[j] < 0: #cosfp or sinfp?
                            velocity_orbplane = [-1 * v for v in velocity_orbplane]
                        # dydx = (cosfp[j] + self.e[i]) / (sinfp[j])
                        # velocity_orbplane = [-1, dydx, 0]
                        # ratio_ = velocity / np.linalg.norm(velocity_orbplane)
                        # velocity_orbplane = [v * ratio_ for v in velocity_orbplane]
                        # if sinfp[j] < 0:  # cosfp or sinfp?
                        #     velocity_orbplane = [-1 * v for v in velocity_orbplane]

                    # print(old_orbplane)
                    # print(velocity_orbplane)
                    # print(velocity_orbplane/np.linalg.norm(velocity_orbplane))
                    # pdb.set_trace()







                    velocity_eq = M @ velocity_orbplane
                    coords_orbplane = [radius * cosfp[j], radius * sinfp[j], 0]
                    coords_eq = M @ coords_orbplane

                    # print(velocity_eq)
                    # print(coords_eq)
                    # print(self.Omega[i])

                    #hard-coding for double collision
                    if "hardcode" in self.inputdata:
                        if self.inputdata["hardcode"] == 1:
                            coords_eq = np.array([-radius * np.cos(self.I[i]), 0, -radius * np.sin(self.I[i])])
                            velocity_eq = np.array([0, -velocity, 0])
                        elif self.inputdata["hardcode"] == 2:
                            coords_eq = np.array([radius * np.cos(self.I[i]), 0, -radius * np.sin(self.I[i])])
                            velocity_eq = np.array([0, velocity, 0])

                    start_time_func = time.time()
                    a, e, I, O, w, f = self.get_orbital_elements_rand_dv(coords_eq, velocity_eq, dv_ratio, mu, 1e-40)
                    end_time_func = time.time()
                    time_in_func += end_time_func - start_time_func
                    self.a_initial[j] = a
                    self.e_initial[j] = e
                    self.I_initial[j] = I
                    self.Omega_initial[j] = O
                    self.omega_initial[j] = w
                    self.sinf_initial[j] = np.sin(f)
                    self.cosf_initial[j] = np.cos(f)


            end_ejecta = time.time()
            ejecta_time += end_ejecta - start_ejecta
            start_beta_calcs = time.time()

            self.a_dust, self.e_dust, self.I_dust, self.Omega_dust, self.omega_dust, self.beta_dust = \
                bd.OrbTimeCorr_MidOptimized(a_launch=self.a_initial, e_launch=self.e_initial, I_launch=self.I_initial,
                                            Omega_launch=self.Omega_initial, omega_launch=self.omega_initial, cosf_launch=self.cosf_initial,
                                            sinf_launch=self.sinf_initial, beta_per_launch=self.beta_per_launch, stabfac=self.stabfac, beta_limit=betamax)
            end_beta_calcs = time.time()
            beta_calcs_time += end_beta_calcs - start_beta_calcs


            uboundi = np.where(self.a_dust[:] < 0)[0]
            if len(uboundi) > 0 and self.inputdata["betadistrb"] != 0:
                pdb.set_trace()

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
        print("Time spent computing matrices:  ", matrix_time)
        print("Time spent computing orbital elements and ejecta:  ", ejecta_time)
        print("Time spent specifically in orb_elements function:    ", time_in_func)
        print("Time spent computing betas:  ", beta_calcs_time)
        end_time = time.time()
        print("Total time: ", end_time - start_time)
        print(len(self.a_dust))

    def ComputeForkDust_Optimized3(self, fileName):
        code_start = time.time()
        # launch sites --- NO radiation pressure
        Nlaunch = int(self.inputdata["Nfork"])
        I_launch = np.ones(Nlaunch) * 0.05
        # a_launch = np.ones(Nlaunch) * 200.
        # e_launch = np.ones(Nlaunch) * 0.7
        a_launch = np.ones(Nlaunch) * self.inputdata["afork"]
        e_launch = np.ones(Nlaunch) * self.inputdata["efork"]
        # omega_launch = np.ones(Nlaunch)*-1.*np.pi/4.
        # omega_launch = np.random.uniform(0., 2.*np.pi, Nlaunch)
        Omega_launch = np.random.uniform(0., 2. * np.pi, Nlaunch)
        if "exp_C" in self.inputdata:
            # NEW CODE BLOCK
            # Laplace-Lagrange collisional family modification
            ecc_forced = 0.3
            pomega_forced = 0.
            h0 = ecc_forced * np.cos(pomega_forced)
            k0 = ecc_forced * np.sin(pomega_forced)
            ecc_free = 0.05
            pomega_free = np.random.uniform(0., 2. * np.pi, Nlaunch)
            h = h0 + ecc_free * np.cos(pomega_free)
            k = k0 + ecc_free * np.sin(pomega_free)
            e_launch = np.sqrt(h ** 2 + k ** 2)
            pomega_launch = np.arctan2(k, h)
            print(pomega_launch * 180. / np.pi)
        elif "exp_B" in self.inputdata:
            pomega_launch = np.random.uniform(-np.pi / 8., np.pi / 8., Nlaunch)
        else:
            pomega_launch = 0.0 * np.pi / 2.

        omega_launch = pomega_launch - Omega_launch
        if "exp_A" in self.inputdata or "exp_B" in self.inputdata or "exp_C" in self.inputdata:
            f_launch = np.random.uniform(-np.pi / 2., np.pi / 2., Nlaunch)
        else:
            f_launch = np.random.uniform(0., 2. * np.pi, Nlaunch)
        cosf_launch = np.cos(f_launch)
        sinf_launch = np.sin(f_launch)

        ## beta distribution parameters
        if "stabfac_fork" in self.inputdata:
            stabfac = self.inputdata[
                "stabfac_fork"]  # needs to be less than 1 since integral diverges. if 0.997 then max Q_dust ~ 20000 au
        else:
            stabfac = 0.997
        betamin = 0.001  # betamin is a scalar
        betamax = (1. - e_launch ** 2) / (2. * (1. + e_launch * cosf_launch))  # betamax is an array!
        betamax = betamax * stabfac  # otherwise integral diverges

        betapow = 1.5
        beta_per_launch = int(self.inputdata["beta_per_fork"])  # number of dust particles = beta_per_launch * Nlaunch

        a_dust = np.zeros((Nlaunch, beta_per_launch))
        e_dust = np.zeros((Nlaunch, beta_per_launch))
        I_dust = np.zeros((Nlaunch, beta_per_launch))
        Omega_dust = np.zeros((Nlaunch, beta_per_launch))
        omega_dust = np.zeros((Nlaunch, beta_per_launch))
        beta_dust = np.zeros((Nlaunch, beta_per_launch))

        for j in range(len(betamax)):  # loop over Nlaunch directions (betamax.size = Nlaunch)
            print(f"Computing {j}/{len(betamax)} fork parent bodies...")

            e_scalar = e_launch[
                j]  # creating these scalars is not strictly necessary but there seems to be slight speed advantage for dNdbeta
            cosf_scalar = cosf_launch[j]
            betamax_scalar = betamax[j]

            # dNdbeta should account for non-zero e and non-unity cosf_scalar
            dNdbeta = lambda beta: beta ** betapow * (1 - beta) ** 1.5 * (
                    1 - e_scalar ** 2 - 2 * beta * (1 + e_scalar * cosf_scalar)) ** -1.5

            # create a smooth beta-grid that concentrates precision near betamax_scalar
            how_close = 1.e-3  # how fractionally close the SECOND-TO-LAST beta grid point comes to betamax_scalar
            n_beta_grid = 50  # number of points in the beta-grid = number of points in the cumulative Nbeta
            epsmax = np.log10(betamax_scalar - betamin)
            epsmin = np.log10(how_close * betamax_scalar)
            eps = np.linspace(epsmax, epsmin, n_beta_grid - 1)
            beta = betamax_scalar - 10. ** eps
            beta = np.append(beta,
                             betamax_scalar)  # need to include betamax as part of the array for later interpolation to work
            # beta = np.linspace(betamin,betamax_scalar,1000)

            norm = sint.quad(dNdbeta, betamin, betamax_scalar)[0]
            Nbeta = lambda beta: sint.quad(dNdbeta, betamin, beta)[0] / norm
            Nbeta = np.vectorize(Nbeta)

            x = nr.uniform(0, 1, beta_per_launch)
            invNbeta = si.interp1d(Nbeta(beta), beta)

            beta_segment = invNbeta(x)  # there are beta_per_launch values in beta_segment
            # plt.hist(beta_segment,bins=100,log=True)

            # calculate orbital elements corresponding to the jth direction, varying beta
            # beg = j * beta_per_launch
            # end = beg + beta_per_launch
            a_dust[j] = (1 - beta_segment) * a_launch[j] * (1 - e_launch[j] ** 2) / (
                    1 - e_launch[j] ** 2 - 2 * beta_segment * (1 + e_launch[j] * cosf_launch[j]))
            omega_dust[j] = omega_launch[j] + np.arctan2(beta_segment * sinf_launch[j],
                                                         e_launch[j] + beta_segment * cosf_launch[j])
            e_dust[j] = np.sqrt(
                e_launch[j] ** 2 + 2 * beta_segment * e_launch[j] * cosf_launch[j] + beta_segment ** 2) / (
                                1 - beta_segment)
            I_dust[j] = I_launch[j] * np.ones(beta_per_launch)
            Omega_dust[j] = Omega_launch[j] * np.ones(beta_per_launch)
            beta_dust[j] = beta_segment
            # end for loop over j

        code_stop = time.time()
        print('Number of dust particles = Nlaunch * beta_per_launch = ', Nlaunch * beta_per_launch)
        print('Time to run code:', code_stop - code_start)

        # omega_dust must be calculated before e_dust is updated!

        # plt.hist(Q_dust[np.where(Q_dust < 2000)],bins=100)
        # plt.hist(Q_dust[np.where(beta_dust > 0.47)],bins=1000)

        # print(a_dust,e_dust,I_dust,Omega_dust,omega_dust)

        # true_anom = np.random.uniform(0, 2. * np.pi,
        #                               Nlaunch * beta_per_launch)  # JOSH: THIS LINE IS FOR CONVENIENCE AND SHOULD BE REPLACED BY A PROPER SOLUTION OF KEPLER'S EQUATION THAT SAMPLES MEAN ANOMALIES RANDOMLY AND THEN CONVERTS TO TRUE ANOMALIES.
        a_dust = a_dust.flatten()
        omega_dust = omega_dust.flatten()
        e_dust = e_dust.flatten()
        I_dust = I_dust.flatten()
        Omega_dust = Omega_dust.flatten()
        beta_dust = beta_dust.flatten()

        # R = a_dust * (1 - e_dust ** 2) / (1. + e_dust * np.cos(true_anom))
        # X = R * (np.cos(Omega_dust) * np.cos(omega_dust + true_anom) - np.sin(Omega_dust) * np.sin(
        #     omega_dust + true_anom) * np.cos(I_dust))
        # Y = R * (np.sin(Omega_dust) * np.cos(omega_dust + true_anom) + np.cos(Omega_dust) * np.sin(
        #     omega_dust + true_anom) * np.cos(I_dust))
        # Z = R * np.sin(omega_dust + true_anom) * np.sin(I_dust)
        # np.savetxt("/Users/sjosh/PycharmProjects/Research/img_2/debris-disk/parentorbit/" + f"{fileName}_X.txt", X)
        # np.savetxt("/Users/sjosh/PycharmProjects/Research/img_2/debris-disk/parentorbit/" + f"{fileName}_Y.txt", Y)
        # np.savetxt("/Users/sjosh/PycharmProjects/Research/img_2/debris-disk/parentorbit/" + f"{fileName}_Z.txt", Z)

        # mean_anom = np.random.uniform(0, 2.*np.pi, Nlaunch * beta_per_launch)
        # true_anom = np.zeros(Nlaunch * beta_per_launch)
        # for idx in range(len(mean_anom)):
        #     if idx % 100 == 0:
        #         print("IDX:  ", idx)
        #     true_anom[idx], _ = ppo.OutputPosition(e=e_dust[idx], Npts=1)
        #
        # R = a_dust * (1 - e_dust ** 2) / (1. + e_dust * np.cos(true_anom))
        # X = R * (np.cos(Omega_dust) * np.cos(omega_dust + true_anom) - np.sin(Omega_dust) * np.sin(
        #     omega_dust + true_anom) * np.cos(I_dust))
        # Y = R * (np.sin(Omega_dust) * np.cos(omega_dust + true_anom) + np.cos(Omega_dust) * np.sin(
        #     omega_dust + true_anom) * np.cos(I_dust))
        # Z = R * np.sin(omega_dust + true_anom) * np.sin(I_dust)
        # np.savetxt("/Users/sjosh/PycharmProjects/Research/img_2/debris-disk/parentorbit/" + f"{fileName}_Xtrue.txt", X)
        # np.savetxt("/Users/sjosh/PycharmProjects/Research/img_2/debris-disk/parentorbit/" + f"{fileName}_Ytrue.txt", Y)
        # np.savetxt("/Users/sjosh/PycharmProjects/Research/img_2/debris-disk/parentorbit/" + f"{fileName}_Ztrue.txt", Z)
        #
        # plt.figure(figsize=(10, 4))
        # plt.xlim([-2000, 1000])
        # plt.ylim([-200, 200])
        # plt.scatter(X, Z, s=0.01, color='black')
        self.a_dust = a_dust
        self.e_dust = e_dust
        self.I_dust = I_dust
        self.Omega_dust = Omega_dust
        self.omega_dust = omega_dust
        self.beta_dust = beta_dust
        self.SaveValues()


    def ComputeBackgroundDustGrains_Optimized(self, manual=False, beta=0.3, Nlaunch=10):
        # Compute orbital parameters of launched dust grains
        # beta = Prad/Pgrav
        # Nlaunch = launch points per parent body orbit
        # print("Computing Dust Grain Orbits...")
        # ct = datetime.datetime.now()
        # print("current time:-", ct)
        if "Nlaunchback" not in self.inputdata or "Nback" not in self.inputdata:
            return
        print("Computing Background Dust Grain Orbits...")
        start_time = time.time()
        if manual:
            self.beta = beta
            Nlaunchback = Nlaunchback
        elif self.inputdata["betadistrb"] == 0:
            self.beta = self.inputdata["beta"]
            Nlaunchback = int(self.inputdata["Nlaunchback"])
            self.beta_dust = np.ones((len(self.h), Nlaunchback)) * self.beta
        else:
            Nlaunchback = int(self.inputdata["Nlaunchback"])
            betapow = self.inputdata["betadistrb"]
            betamin, betamax = self.inputdata["betamin"], self.inputdata["betamax"]
            self.beta_dust = np.zeros((len(self.h), Nlaunchback, self.beta_per_launch_back))

        self.a_dust = np.zeros((len(self.h), Nlaunchback, self.beta_per_launch_back))  # len(self.h) is Nback - this is Nback x Nlaunchback array
        self.e_dust = np.zeros((len(self.h), Nlaunchback, self.beta_per_launch_back))
        self.I_dust = np.zeros((len(self.h), Nlaunchback, self.beta_per_launch_back))
        self.Omega_dust = np.zeros((len(self.h), Nlaunchback, self.beta_per_launch_back))
        self.omega_dust = np.zeros((len(self.h), Nlaunchback, self.beta_per_launch_back))

        lps = np.zeros((len(self.h), Nlaunchback))
        mu = consts.G * self.Mstar
        matrix_time = 0
        ejecta_time = 0
        beta_calcs_time = 0
        beta_application_time = 0
        time_in_func = 0




        for i in range(len(self.h)):  # for each parent body (Nback)
            print("%i/%i background parent body" % (i + 1, len(self.h)))
            fp = list(self.f_vals)
            #fp = nr.uniform(0, 2 * np.pi, Nlaunchback)  # get random true anomaly for each
            lps[i] = fp
            cosfp = np.cos(fp)
            sinfp = np.sin(fp)

            # initial: orbital elements of dust grains before applying radiation pressure
            for j in range(len(fp)): #for each launch point
                # launch beta_per_launch_back points
                for k in range(self.beta_per_launch_back):
                    beta = self.inv_map[round(self.e[i],5)][round(cosfp[j], 5)](nr.uniform(0, 1))
                    self.a_dust[i][j][k] = (1 - beta) * self.a[i] * (1 - self.e[i] ** 2) / (
                            1 - self.e[i] ** 2 - 2 * beta * (1 + self.e[i] * cosfp[j]))
                    self.omega_dust[i][j][k] = self.omega[i] + np.arctan2(beta * sinfp[j],
                                                                       self.e[i] + beta * cosfp[j])
                    self.e_dust[i][j][k] = np.sqrt(
                        self.e[i] ** 2 + 2 * beta * self.e[i] * cosfp[j] + beta ** 2) / (
                                                  1 - beta)
                    self.I_dust[i][j][k] = self.I[j]
                    self.Omega_dust[i][j][k] = self.Omega[j]
                    self.beta_dust[i][j][k] = beta
        self.a_dust = self.a_dust.flatten()
        self.e_dust = self.e_dust.flatten()
        self.I_dust = self.I_dust.flatten()
        self.Omega_dust = self.Omega_dust.flatten()
        self.omega_dust = self.omega_dust.flatten()
        self.beta_dust = self.beta_dust.flatten()







        uboundi = np.where(self.a_dust[:] < 0)[0]
        if len(uboundi) > 0 and self.inputdata["betadistrb"] != 0:
            pdb.set_trace()

        lps = lps.flatten()
        np.savetxt('launchpoints.txt', lps)

        # self.a_dust = self.a_dust.flatten()
        # self.e_dust = self.e_dust.flatten()
        # self.I_dust = self.I_dust.flatten()
        # self.Omega_dust = self.Omega_dust.flatten()
        # self.omega_dust = self.omega_dust.flatten()
        # self.beta_dust = self.beta_dust.flatten()

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
        # print("Time spent computing betas:  ", beta_calcs_time)
        end_time = time.time()
        print("Total time: ", end_time - start_time)
        print(len(self.a_dust))

    def ComputeBackgroundDustGrains_BetaOptimized2(self, manual=False, beta=0.3, Nlaunch=10):
        # Compute orbital parameters of launched dust grains
        # beta = Prad/Pgrav
        # Nback: number of parent background orbits
        # Nlaunchback = launch points per parent body orbit
        if "Nlaunchback" not in self.inputdata or "Nback" not in self.inputdata:
            return
        print("Computing Background Dust Grain Orbits...")
        start_time = time.time()
        if manual:
            self.beta = beta
            Nlaunchback = Nlaunchback
        elif self.inputdata["betadistrb"] == 0:
            self.beta = self.inputdata["beta"]
            Nlaunchback = int(self.inputdata["Nlaunchback"])
            self.beta_dust = np.ones((len(self.h), Nlaunchback)) * self.beta
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
        mu = consts.G * self.Mstar
        matrix_time = 0
        ejecta_time = 0
        beta_calcs_time = 0
        beta_application_time = 0
        time_in_func = 0
        self.e_max=  0.02
        self.cosf_max = 1
        start_beta_calcs = time.time()
        end_beta_calcs = time.time()
        beta_calcs_time += end_beta_calcs - start_beta_calcs
        for i in range(len(self.h)):  # for each parent body
            print("%i/%i background parent body" % (i + 1, len(self.h)))
            fp = nr.uniform(0, 2 * np.pi, Nlaunchback)  # get random true anomaly for each
            lps[i] = fp
            cosfp = np.cos(fp)
            sinfp = np.sin(fp)


            self.e_max = self.e[i]
            self.cosf_max = np.max(cosfp)

            start_beta_calcs = time.time()
            inverseCDF = bd.get_inverse_CDF_stab(e=self.e_max, cosf=self.cosf_max, betapow=betapow,
                                                 betamin=betamin, betamax=betamax, stabfac=self.stabfac_back, precision=500)

            self.beta_dust[i] = bd.OrbTimeCorr_Vectorized(inverseCDF, Nlaunchback)
            end_beta_calcs = time.time()
            beta_calcs_time += end_beta_calcs - start_beta_calcs

            start_beta_application_time = time.time()
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
            end_beta_application_time = time.time()
            beta_application_time += end_beta_application_time - start_beta_application_time

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
        print("Time spent computing matrices:  ", matrix_time)
        print("Time spent computing betas:  ", beta_calcs_time)
        print("Time applying radiation pressure formulae:    ", beta_application_time)
        end_time = time.time()
        print("Total time: ", end_time - start_time)
        #print("Betamax:     ", approx_betamax)


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
    def OutputQbeta(self, outfile):
        Qmax = self.a_dust * (1 + self.e_dust)
        Qbeta = np.array([self.beta_dust, Qmax])
        np.savetxt(outfile + "_qbeta.txt", list(Qbeta))

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
        for i in range(len(self.a_dust.flatten())):
           fd.draw_particles(ax, self.Omega_dust.flatten()[i], self.I_dust.flatten()[i], \
                          argu_peri=self.omega_dust.flatten()[i], e=self.e_dust.flatten()[i], \
                             r0=self.a_dust.flatten()[i], format='ro', msz=0.5) #draw dust particles
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

    def ComputeForkParentSingle(self, manual=False, Nback=500, amin=4., amax=6., I0=5. * (np.pi / 180.),
                                                e0=0.2,
                                                Omega0=0., omega0=0., random=False):
        # Compute p, q, h, k of parent bodies for background disk
        # Single planet only
        # Nback: number of parent bodies
        # amin, amax: min, max semi-major axes of parent bodies
        # I0, e0: initial mutual inclination and eccentricities
        # Omega0, omega0: initial nodal angle and argument of periapse
        if "Nforkback" not in self.inputdata or "Nfork" not in self.inputdata:
            return

        Nfork = int(self.inputdata["Nfork"])
        print("Computing Fork Parent Orbits (single planet)...")
        self.a = self.inputdata["afork"] * np.ones(Nfork)
        #self.a = np.array([self.inputdata["afork"] for _ in range(Nfork)])
        self.Ifork = self.inputdata["Ifork"]
        self.efork = self.inputdata["efork"]
        self.Omega = np.linspace(0, 2 * np.pi, Nfork)
        self.omega = -self.Omega


        self.A = ff.A(self.a, self.Mp[0], self.ap[0], Mstar=self.Mstar)
        self.Aj = ff.Aj(self.a, self.Mp[0], self.ap[0], Mstar=self.Mstar)
        self.B = ff.B(self.a, self.Mp[0], self.ap[0], Mstar=self.Mstar)
        self.Bj = ff.Bj(self.a, self.Mp[0], self.ap[0], Mstar=self.Mstar)

        # forced component
        self.h0 = -(self.Aj / self.A) * self.ep[0] * np.sin(self.Omegap[0] + self.omegap[0])
        self.k0 = -(self.Aj / self.A) * self.ep[0] * np.cos(self.Omegap[0] + self.omegap[0])
        self.p0 = -(self.Bj / self.B) * self.Ip[0] * np.sin(self.Omegap[0])
        self.q0 = -(self.Bj / self.B) * self.Ip[0] * np.cos(self.Omegap[0])


        # self.h = efree * np.sin(omega + Omega) + self.h0
        # self.k = efree * np.cos(omega + Omega) + self.k0
        # self.p = Ifree * np.sin(Omega) + self.p0
        # self.q = Ifree * np.cos(Omega) + self.q0
        self.h = self.h0
        self.k = self.k0
        self.p = self.p0
        self.q = self.q0


    def ComputeForkParentOrbital(self):
        # Compute orbital parameters of parent bodies
        # compare with what we expect
        N = len(self.p)
        self.Omega = np.random.uniform(0, 2*np.pi, N)
        self.omega = np.zeros(N) - self.Omega
        self.I = 0.05 * np.ones(N) #np.array([0.05])
        self.e = 0.7 * np.ones(N) #np.array([0.7])
        # self.Omega = np.arctan2(self.p, self.q)
        # pomega = np.arctan2(self.h, self.k)
        # self.omega = pomega - self.Omega
        # self.I = np.sqrt(self.p ** 2 + self.q ** 2)
        # self.e = np.sqrt(self.h ** 2 + self.k ** 2)

    def ComputeBackgroundParentSingle_Optimized(self, manual=False, Nback=500, amin=4., amax=6., I0=5. * (np.pi / 180.), e0=0.2,
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
            print("IF")
            self.a, self.I0, self.e0, self.Omega0, self.omega0 = np.loadtxt(self.freeelemtxt, unpack=True)
        elif manual:
            print('ELIF MANUAL')
            self.a = np.linspace(amin, amax, int(self.inputdata["Nback"]))
            self.I0 = I0
            self.e0 = e0
            self.Omega0 = Omega0
            self.omega0 = omega0
        else:
            print("ELSE")
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
            #num_f should be Nlaunchback? change parameters



            # e_vals, f_vals, inv_map = bd.OrbTimeCorr_MidOptimized_Background(num_e=self.num_e,
            #                                                                  max_e=self.inputdata["e0"],
            #                                                                  num_f=self.Nlaunchback,
            #                                                                  n_beta_grid=50, betapow=1.5, betamin=0.001,
            #                                                                  stabfac=self.stabfac_back)
            #
            # self.e_vals = e_vals
            # self.f_vals = f_vals
            # self.inv_map = inv_map
            e_vals = np.linspace(0, self.inputdata["e0"], self.num_e)
            efree = np.zeros(int(self.inputdata["Nback"]))
            for i in range(len(efree)):
                e = e_vals[np.random.randint(0, len(e_vals))]
                efree[i] = e
            #efree = nr.uniform(0, self.inputdata["e0"], int(self.inputdata["Nback"]))
            Ifree = nr.uniform(
                np.max((0, self.inputdata["Icent"] * np.pi / 180. - self.inputdata["I0"] * np.pi / 180.)), \
                self.inputdata["I0"] * np.pi / 180. + self.inputdata["Icent"] * np.pi / 180.,
                int(self.inputdata["Nback"]))
            #print(Ifree)
            omega = nr.uniform(0, 2 * np.pi, int(self.inputdata["Nback"]))
            Omega = nr.uniform(0, 2 * np.pi, int(self.inputdata["Nback"]))
            self.h = efree * np.sin(omega + Omega) + self.h0
            self.k = efree * np.cos(omega + Omega) + self.k0
            self.p = Ifree * np.sin(Omega) + self.p0
            self.q = Ifree * np.cos(Omega) + self.q0
        print("MIN, MAX", amin, amax)
        # 0.02


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
            #print(Ifree)
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

    def ComputeBackgroundParentOrbital_Optimized(self):
        # Compute orbital parameters of parent bodies
        self.Omega = np.arctan2(self.p, self.q)
        pomega = np.arctan2(self.h, self.k)
        self.omega = pomega - self.Omega
        self.I = np.sqrt(self.p ** 2 + self.q ** 2)
        self.e = np.sqrt(self.h ** 2 + self.k ** 2)
        e_vals, f_vals, inv_map = bd.OrbTimeCorr_MidOptimized_Background_2(e_set=set(self.e),
                                                                         num_f=self.Nlaunchback,
                                                                         n_beta_grid=50, betapow=1.5, betamin=0.001,
                                                                         stabfac=self.stabfac_back)

        self.e_vals = e_vals
        self.f_vals = f_vals
        self.inv_map = inv_map
        # print("BACKGROUND OMEGA: ", self.Omega)
        # for _ in range(10):
        #     print('--------------------------------------')
        # print(self.I)
        # for _ in range(10):
        #     print('--------------------------------------')
        # print(self.e)



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
                print("Circular, equatorial case")
                pdb.set_trace()
                f = np.arccos(r0[0] / radius)
                if v0[0] > 0:
                    f = 2 * np.pi - f
                w = 0  # "periapsis" at launch
                Omega = 0
            else:
                print("Circular, inclined case")
                pdb.set_trace()
                f = np.arccos((node_vec @ r0) / (n * radius))
                if r0[2] < 0:
                    f = 2 * np.pi - f
                w = f # SUSPECT - w = 0? -f
                Omega = np.arccos(n_x / n)
        else: # eccentric
            if abs(I) < eps or abs(I - np.pi) < eps: # equatorial
                # print("Eccentric, equatorial case")
                # pdb.set_trace()
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

            if e_vec[2] < 0:
                w = 2 * np.pi - w


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