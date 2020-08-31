import numpy as np
import os
import errno
import pickle as pc
import hibpplotlib as hbplot
import copy
from matplotlib import path
from scipy.interpolate import RegularGridInterpolator
from matplotlib.patches import Rectangle

# %% define class for trajectories


class Traj():

    def __init__(self, q, m, Ebeam, r0, alpha, beta, U, dt=1e-7):
        ''' class for trajectories
        q - particle charge [Co]
        m - particle mass [kg]
        Ebeam - beam energy [keV]
        r0 - initial position [m]
        alpha - injection angle in XY plane [rad]
        beta - injection angle in XZ plane [rad]
        U - list of voltages in [kV] [A2 B2 A3 B3]
        dt - timestep for RK algorithm [s]
        '''
        self.q = q
        self.m = m
        self.Ebeam = Ebeam
        drad = np.pi/180.  # converts degrees to radians
        Vabs = np.sqrt(2 * Ebeam * 1.602176487E-16 / m)
        V0 = np.array([-Vabs * np.cos(alpha*drad)*np.cos(beta*drad),
                       -Vabs * np.sin(alpha*drad),
                       Vabs * np.cos(alpha*drad)*np.sin(beta*drad)])
        self.alpha = alpha
        self.beta = beta
        self.U = U
        self.RV0 = np.array([np.hstack((r0, V0))])  # initial condition

        self.RV_prim = self.RV0  # array with r,V for the whole trajectory
        self.tag_prim = [1]

        self.RV_sec = np.array([[]])
        self.tag_sec = [2]

        # list to contain RV of the whole fan
        self.Fan = []
        # time step for primary orbit
        self.dt1 = dt
        # time step for secondary orbit
        self.dt2 = dt
        self.IsAimXY = False
        self.IsAimZ = False
        self.IntersectGeometry = False
        self.IntersectGeometrySec = False
        self.B_out_of_bounds = False
        # multislit
        self.fan_to_slits = []
        self.RV_sec_toslits = []
        self.ion_zones = []

    def pass_prim(self, E_interp, B_interp, geom, tmax=1.0):
        ''' passing primary trajectory from initial point self.RV0
            geom - Geometry object
        '''
        print('\n Passing primary trajectory')
        self.IntersectGeometry = False
        t = 0.
        dt = self.dt1
        RV_old = self.RV0  # initial position
        RV = self.RV0  # array to collect all r,V
        tag_column = [10]

        while t <= tmax:
            r = RV_old[0, :3]
            # Electric field
            E_local = return_E(r, E_interp, self.U)
            # Magnetic field
            try:
                B_local = return_B(r, B_interp)
            except ValueError:
                print('Btor Out of bounds for primaries, r = ', r)
                break
            # runge-kutta step
            RV_new = runge_kutt(self.q, self.m, RV_old, dt, E_local, B_local)
            RV = np.vstack((RV, RV_new))

            tag_column = np.hstack((tag_column, 10))

            if geom.check_chamb_ent_intersect(RV_old[0, 0:3], RV_new[0, 0:3]):
                print('Primary intersected chamber entrance')
                self.IntersectGeometry = True
                break

            if geom.check_plates_intersect(RV_old[0, 0:3], RV_new[0, 0:3]):
                print('Primary intersected plates')
                self.IntersectGeometry = True
                break

            RV_old = RV_new
            t = t + dt
            # print('t = ', t)

        self.RV_prim = RV
        self.tag_prim = tag_column

    def pass_sec(self, RV0, r_aim, E_interp, B_interp, geom,
                 stop_plane_n=np.array([1, 0, 0]), tmax=7e-5,
                 eps_xy=1e-3, eps_z=1e-3):
        ''' passing secondary trajectory from initial point RV0 to point r_aim
            with accuracy eps
            geom - Geometry object
        '''
        # print('Passing secondary trajectory')
        self.IsAimXY = False
        self.IsAimZ = False
        self.IntersectGeometrySec = False
        self.B_out_of_bounds = False
        t = 0
        dt = self.dt2
        RV_old = RV0  # initial position
        RV = RV0  # array to collect all r,V
        tag_column = [20]

        while t <= tmax:  # to witness curls initiate tmax as 1 sec
            r = RV_old[0, :3]
            # Electric field
            E_local = return_E(r, E_interp, self.U)
            # Magnetic field
            try:
                B_local = return_B(r, B_interp)
            except ValueError:
                print('Btor Out of bounds for secondaries, r = ',
                      np.round(r, 3))
                self.B_out_of_bounds = True
                break
            # runge-kutta step
            RV_new = runge_kutt(2*self.q, self.m, RV_old, dt, E_local, B_local)

            if geom.check_chamb_ext_intersect(RV_old[0, 0:3], RV_new[0, 0:3]):
                # print('Secondary intersected chamber exit')
                self.IntersectGeometrySec = True

            if geom.check_plates_intersect(RV_old[0, 0:3], RV_new[0, 0:3]):
                print('Secondary intersected plates')
                self.IntersectGeometrySec = True

            if (RV_new[0, 0] > 0.5) and (RV_new[0, 1] < 0.4):
                # find last point of the secondary trajectory
                # intersection with the stop plane
                planeNormal = stop_plane_n
                planePoint = r_aim
                rayDirection = RV_new[0, :3] - RV_old[0, :3]
                rayPoint = RV_new[0, :3]
                r_intersect = line_plane_intersect(planeNormal, planePoint,
                                                   rayDirection, rayPoint)
                # check if r_intersect is between RV_old and RV_new
                if is_between(RV_old[0, :3], RV_new[0, :3], r_intersect):
                    RV_new[0, :3] = r_intersect
                    RV = np.vstack((RV, RV_new))
                    # check XY plane
                    if (np.linalg.norm(RV_new[0, :2] - r_aim[:2]) <= eps_xy):
                        # print('aim XY!')
                        self.IsAimXY = True
                    # check XZ plane
                    if (np.linalg.norm(RV_new[0, [0, 2]] - r_aim[[0, 2]]) <=
                            eps_z):
                        # print('aim Z!')
                        self.IsAimZ = True
                    break

            # continue trajectory calculation
            RV_old = RV_new
            t = t + dt
            RV = np.vstack((RV, RV_new))
            tag_column = np.hstack((tag_column, 20))
            # print('t secondary = ', t)

        self.RV_sec = RV
        self.tag_sec = tag_column

    def pass_fan(self, r_aim, E_interp, B_interp, geom,
                 stop_plane_n=np.array([1, 0, 0]),
                 r_plasma=0.3, R=0.36, elon=1.8, eps_xy=1e-3, eps_z=1e-3,
                 no_intersect=False, no_out_of_bounds=False):
        ''' passing fan from initial point self.RV0
        geom - geometry object
        '''
        print('\n Passing fan of trajectories')
        self.pass_prim(E_interp, B_interp, geom)
        # create a list fro secondary trajectories
        list_sec = []
        # check intersection of primary trajectory
        if self.IntersectGeometry:
            print('Fan list is empty')
            self.Fan = list_sec
            return

        # check eliptical radius of particle
        # R - major radius of a torus, elon - size along Y
        mask = np.sqrt((self.RV_prim[:, 0] - R)**2 +
                       (self.RV_prim[:, 1] / elon)**2) <= r_plasma
        self.tag_prim[mask] = 11

        # list of initial points of secondary trajectories
        RV0_sec = self.RV_prim[(self.tag_prim == 11)]

        for RV02 in RV0_sec:
            RV02 = np.array([RV02])
            self.pass_sec(RV02, r_aim, E_interp, B_interp, geom,
                          stop_plane_n=stop_plane_n,
                          eps_xy=eps_xy, eps_z=eps_z)
            if (no_intersect and self.IntersectGeometrySec) or \
               (no_out_of_bounds and self.B_out_of_bounds):
                continue
            list_sec.append(self.RV_sec)

        self.Fan = list_sec

    def add_slits(self, n_slits):
        ''' empty list for secondary trajectories, which go to different slits
        '''
        if len(self.RV_sec_toslits) == n_slits:
            pass
        else:
            self.RV_sec_toslits = [None]*n_slits
            self.ion_zones = [None]*n_slits

    def plot_prim(self, ax, axes='XY', color='k', full_primary=False):
        ''' plot primary trajectory
        '''
        axes_dict = {'XY': (0, 1), 'XZ': (0, 2), 'ZY': (2, 1)}
        index_X, index_Y = axes_dict[axes]
        if full_primary:
            index = -1
        else:
            # find where secondary trajectory starts
            for i in range(self.RV_prim.shape[0]):
                if np.linalg.norm(self.RV_prim[i, :3]
                                  - self.RV_sec[0, :3]) < 1e-4:
                    index = i+1
        ax.plot(self.RV_prim[:index, index_X],
                self.RV_prim[:index, index_Y],
                color=color, linewidth=2)

    def plot_sec(self, ax, axes='XY', color='r'):
        ''' plot secondary trajectory
        '''
        axes_dict = {'XY': (0, 1), 'XZ': (0, 2), 'ZY': (2, 1)}
        index_X, index_Y = axes_dict[axes]
        ax.plot(self.RV_sec[:, index_X], self.RV_sec[:, index_Y],
                color=color, linewidth=2)

    def plot_fan(self, ax, axes='XY', color='r'):
        ''' plot fan of secondary trajectories
        '''
        axes_dict = {'XY': (0, 1), 'XZ': (0, 2), 'ZY': (2, 1)}
        index_X, index_Y = axes_dict[axes]
        for i in self.Fan:
            ax.plot(i[:, index_X], i[:, index_Y], color=color)


# %% define class for geometry
class Geometry():
    '''
    object containing geometry points
    '''

    def __init__(self):
        # lists with coordinates of 4 points,
        # determining chamber entrance and exit
        self.chamb_ent = []
        self.chamb_ext = []
        # dictionary for arrays of plates coordinates
        self.plates_edges = dict()
        # dictionary for positions of all plates
        self.r_dict = dict()
        # dictionary for poloidal field coils
        self.pf_coils = dict()
        # Tor Field coil contour
        self.coil = np.array([])
        # vacuum vessel contour
        self.camera = np.array([])
        # separatrix contour
        self.sep = np.array([])
        # inner and outer first wall contours
        self.in_fw = np.array([])
        self.out_fw = np.array([])
        # array for slits coordinates
        self.slits_edges = np.array([])
        # array for slit normal
        self.slit_plane_n = np.array([])
        # array for slit polygon
        self.slits_spot = np.array([])
        # arrays for primary and secondary beamline angles
        self.prim_angles = np.array([])
        self.sec_angles = np.array([])

    def check_chamb_ent_intersect(self, point1, point2):
        ''' check the intersection between segment 1->2 and the chamber
        '''
        if (segm_intersect(point1[0:2], point2[0:2],
                           self.chamb_ent[0], self.chamb_ent[1])) or \
            (segm_intersect(point1[0:2], point2[0:2],
                            self.chamb_ent[2], self.chamb_ent[3])):
            return True
        else:
            return False

    def check_chamb_ext_intersect(self, point1, point2):
        ''' check the intersection between segment 1->2 and the chamber
        '''
        if (segm_intersect(point1[0:2], point2[0:2],
                           self.chamb_ext[0], self.chamb_ext[1])) or \
            (segm_intersect(point1[0:2], point2[0:2],
                            self.chamb_ext[2], self.chamb_ext[3])):
            return True
        else:
            return False

    def check_plates_intersect(self, point1, point2):
        segment_coords = np.array([point1, point2])
        for key in self.plates_edges:
            # maxpoint = np.max(np.max(self.plates_edges[key], axis=1), axis=0)
            # minpoint = np.min(np.min(self.plates_edges[key], axis=1), axis=0)

            # if np.all(point1 < maxpoint) and np.all(point1 > minpoint):
            if segm_poly_intersect(self.plates_edges[key][0],
                                   segment_coords) or \
                segm_poly_intersect(self.plates_edges[key][1],
                                    segment_coords):
                return True
            else:
                continue
        return False

    def add_slits(self, n_slits, slit_dist, slit_w, slit_l, slit_gamma):
        ''' add slits to Geometry
        n_slits - number of slits
        slit_dist - distance between slits [m]
        slit_w - slit width (along Y) [m]
        slit_l - slit length (along Z)
        slit_gamma - angle of ratation around X [deg]
        '''
        # angle of slit line with X-axis in XZ plane
        slit_alpha = self.sec_angles[0]
        # angle of slit line with XZ plane
        slit_beta = self.sec_angles[1]
        # convert degrees to radians

        # calculate slits coordinates
        r_slits = np.zeros([n_slits, 5, 3])
        rs = self.r_dict['slit']
        for i_slit in range(n_slits):
            # calculate coords of slit center
            y0 = (n_slits//2 - i_slit)*slit_dist
            r_slits[i_slit, 0, :] = [0., y0, 0.]
            # calculate slit edges
            r_slits[i_slit, 1, :] = [0., y0 + slit_w/2, slit_l/2]
            r_slits[i_slit, 2, :] = [0., y0 - slit_w/2, slit_l/2]
            r_slits[i_slit, 3, :] = [0., y0 - slit_w/2, -slit_l/2]
            r_slits[i_slit, 4, :] = [0., y0 + slit_w/2, -slit_l/2]
            # rotate and shift to slit position
            for j in range(5):
                r_slits[i_slit, j, :] = rotate(r_slits[i_slit, j, :],
                                               axis=(1, 0, 0), deg=slit_gamma)
                r_slits[i_slit, j, :] = rotate(r_slits[i_slit, j, :],
                                               axis=(0, 0, 1), deg=slit_alpha)
                r_slits[i_slit, j, :] = rotate(r_slits[i_slit, j, :],
                                               axis=(0, 1, 0), deg=slit_beta)
                r_slits[i_slit, j, :] += rs

        # calculate normal to slit plane
        slit_plane_n = np.cross(r_slits[0, 0, :] - r_slits[0, 1, :],
                                r_slits[0, 0, :] - r_slits[0, 2, :])
        slit_plane_n = slit_plane_n/np.linalg.norm(slit_plane_n)

        # create polygon, which contains all slits (slits spot)
        slits_spot = 1.5*np.vstack([r_slits[0, [1, 4], :] - rs,
                                    r_slits[-1, [3, 2], :] - rs]) + rs

        self.slits_edges = r_slits
        self.slit_plane_n = slit_plane_n
        self.slits_spot = slits_spot

    def plot_geom(self, ax, major_radius=1.5, plot_sep=True):
        '''
        plot toroidal and poloidal field coils, camera and
        separatrix in XY plane
        '''
        # plot toroidal coil
        ax.plot(self.coil[:, 0], self.coil[:, 1], '--', color='k')
        ax.plot(self.coil[:, 2], self.coil[:, 3], '--', color='k')

        # get T-15 camera and plasma contours
        for i in range(self.camera.shape[0]):
            ax.plot(self.camera[i, [3, 0]],
                    self.camera[i, [4, 1]], color='tab:blue')

        # plot first wall
        # ax.plot(self.in_fw[:, 0], self.in_fw[:, 1], color='k')
        # ax.plot(self.out_fw[:, 0], self.out_fw[:, 1], color='k')

        if plot_sep:
            ax.plot(self.sep[:, 0], self.sep[:, 1], 'o', markersize=2,
                    color='b')  # 'tab:orange')

        for coil in self.pf_coils.keys():
            xc = self.pf_coils[coil][0]
            yc = self.pf_coils[coil][1]
            dx = self.pf_coils[coil][2]
            dy = self.pf_coils[coil][3]
            ax.add_patch(Rectangle((xc-dx/2, yc-dy/2), dx, dy,
                                   linewidth=1, edgecolor='tab:gray',
                                   facecolor='tab:gray'))

    def plot_plates(self, ax, axes='XY'):
        '''plot plates
        axes='XY', 'XZ', 'ZY'
        '''
        axes_dict = {'XY': (0, 1), 'XZ': (0, 2), 'ZY': (2, 1)}
        index_X, index_Y = axes_dict[axes]
        for name in self.plates_edges.keys():
            ax.fill(self.plates_edges[name][0][:, index_X],
                    self.plates_edges[name][0][:, index_Y], fill=False,
                    hatch='\\', linewidth=2)
            ax.fill(self.plates_edges[name][1][:, index_X],
                    self.plates_edges[name][1][:, index_Y], fill=False,
                    hatch='/', linewidth=2)

    def plot_aim(self, ax, axes='XY', color='r'):
        ''' plot aim dot
        '''
        axes_dict = {'XY': (0, 1), 'XZ': (0, 2), 'ZY': (2, 1)}
        index_X, index_Y = axes_dict[axes]
        ax.plot(self.r_dict['aim'][index_X], self.r_dict['aim'][index_Y],
                '*', color=color)

    def plot_slit(self, ax, axes='XY', color='g'):
        ''' plot slit dot
        '''
        axes_dict = {'XY': (0, 1), 'XZ': (0, 2), 'ZY': (2, 1)}
        index_X, index_Y = axes_dict[axes]
        ax.plot(self.r_dict['slit'][index_X], self.r_dict['slit'][index_Y],
                '*', color=color)

    def plot_slits(self, ax, axes='XY'):
        ''' plot slits contours
        '''
        axes_dict = {'XY': (0, 1), 'XZ': (0, 2), 'ZY': (2, 1)}
        index_X, index_Y = axes_dict[axes]

        r_slits = self.slits_edges
        n_slits = r_slits.shape[0]
        for i in range(n_slits):
            # plot center
            ax.plot(r_slits[i, 0, index_X], r_slits[i, 0, index_Y], '*')
            # plot edge
            ax.fill(r_slits[i, 1:, index_X], r_slits[i, 1:, index_Y],
                    fill=False)


# %% Runge-Kutta
def runge_kutt(q, m, RV, dt, E, B):
    '''
    Calculate one step using Runge-Kutta algorithm
    :param q: particle charge [Co]
    :param m: particle mass [kg]
    :param RV: 7 dimensial vector
           array[[x,y,z,vx,vy,vz,Flag]]
           Flag = 10 primary out of plasma
           Flag = 11 primary in plasma
           Flag = 20 secondary
    :param dt: time discretisation step
    :param E: values of electric field in current point
              np.array([Ex, Ey, Ez])
    :param B: values of magnetic field in current point
              np.array([Bx, By, Bz])
    :return: new RV

     V' = k(E + [VxB]) == K(E + np.cross(V,B)) == f
     r' = V == g

    V[n+1] = V[n] + (h/6)(m1 + 2m2 + 2m3 + m4)
    r[n+1] = r[n] + (h/6)(k1 + 2k2 + 2k3 + k4)
    m[1] = f(t[n], V[n], r[n])
    k[1] = g(t[n], V[n], r[n])
    m[2] = f(t[n] + (h/2), V[n] + (h/2)m[1], r[n] + (h/2)k[1])
    k[2] = g(t[n] + (h/2), V[n] + (h/2)m[1], r[n] + (h/2)k[1])
    m[3] = f(t[n] + (h/2), V[n] + (h/2)m[2], r[n] + (h/2)k[2])
    k[3] = g(t[n] + (h/2), V[n] + (h/2)m[2], r[n] + (h/2)k[2])
    m[4] = f(t[n] + h, V[n] + h*m[3], r[n] + h*k[3])
    k[4] = g(t[n] + h, V[n] + h*m[3], r[n] + h*k[3])

    E - np.array([Ex, Ey, Ez])
    B - np.array([Bx, By, Bz])
    '''
    k = q/m
    r = RV[0, :3]
    V = RV[0, 3:]

    # define equations of movement
    def f(E, V, B): return k*(E + np.cross(V, B))
    def g(V): return V
    ''' m1,k1 '''
    m1 = f(E, V, B)
    k1 = g(V)
    ''' m2,k2 '''
    fV2 = V + (dt / 2.) * m1
    gV2 = V + (dt / 2.) * m1
    m2 = f(E, fV2, B)
    k2 = g(gV2)
    ''' m3,k3 '''
    fV3 = V + (dt / 2.) * m2
    gV3 = V + (dt / 2.) * m2
    m3 = f(E, fV3, B)
    k3 = g(gV3)
    ''' m4,k4 '''
    fV4 = V + dt * m3
    gV4 = V + dt * m3
    m4 = f(E, fV4, B)
    k4 = g(gV4)
    ''' all together! '''
    V = V + (dt / 6.) * (m1 + (2. * m2) + (2. * m3) + m4)
    r = r + (dt / 6.) * (k1 + (2. * k2) + (2. * k3) + k4)

    RV = np.hstack((r, V))

    return RV


# %%
def optimize_B2(tr, r_aim, geom, UB2, dUB2, E, B, dt,
                stop_plane_n, r_plasma, R, elon, eps_xy=1e-3, eps_z=1e-3):
    ''' get voltages on B2 plates and choose secondary trajectory
    which goes into r_aim
    '''
    attempts_high = 0
    attempts_low = 0
    while True:
        tr.U[1], tr.dt1, tr.dt2 = UB2, dt, dt
        # pass fan of trajectories
        tr.pass_fan(r_aim, E, B, geom, stop_plane_n=stop_plane_n,
                    r_plasma=r_plasma, R=R, elon=elon, eps_xy=eps_xy,
                    eps_z=eps_z, no_intersect=True, no_out_of_bounds=True)
        if tr.IntersectGeometry:
            break
        if len(tr.Fan) == 0:
            print('NO secondary trajectories')
            break
        # reset flags in order to let the algorithm work properly
        tr.IsAimXY = False
        # tr.IsAimZ = False
        tr.IntersectGeometrySec = False

        # find which secondaries are higher/lower than r_aim
        # sign = -1 means higher, 1 means lower
        signs = np.array([np.sign(np.cross(RV[-1, :3], r_aim)[-1])
                          for RV in tr.Fan])
        are_higher = np.argwhere(signs == -1)
        are_lower = np.argwhere(signs == 1)
        twisted_fan = False

        if are_higher.shape[0] == 0:
            attempts_high += 1
            n = int(are_lower[are_lower.shape[0]//2])
            if attempts_high > 5:
                attempts_high = 0
                print('Aim is too HIGH along Y!')
                # return tr
                break
        elif are_lower.shape[0] == 0:
            attempts_low += 1
            n = int(are_higher[are_higher.shape[0]//2])
            if attempts_low > 5:
                attempts_low = 0
                print('Aim is too LOW along Y!')
                # return tr
                break
        else:
            attempts_high = 0
            attempts_low = 0
            if are_higher[-1] > are_lower[0]:
                print('Fan is twisted!')
                twisted_fan = True
                n = int(are_lower[-1])
            else:
                n = int(are_higher[-1])  # find one which is higher
        RV_old = np.array([tr.Fan[n][0]])

        # find secondary, which goes directly into r_aim
        tr.dt1 = tr.dt1/2.

        while True:
            # make a small step along primary trajectory
            r = RV_old[0, :3]
            B_local = return_B(r, B)
            E_local = np.array([0., 0., 0.])
            RV_new = runge_kutt(tr.q, tr.m, RV_old, tr.dt1, E_local, B_local)
            # pass new secondary trajectory
            tr.pass_sec(RV_new, r_aim, E, B, geom,
                        stop_plane_n=stop_plane_n, eps_xy=eps_xy, eps_z=eps_z)

            # check XY flag
            if tr.IsAimXY:
                # insert RV_new into primary traj
                index = np.argmin(np.linalg.norm(tr.RV_prim[:, :3] -
                                                 RV_new[0, :3], axis=1))
                tr.RV_prim = np.insert(tr.RV_prim, index+1, RV_new, axis=0)
                tr.tag_prim = np.insert(tr.tag_prim, index+1, 11, axis=0)
                break
            # check if the new secondary traj is lower than r_aim
            if (not twisted_fan and
                    np.sign(np.cross(tr.RV_sec[-1, :3], r_aim)[-1]) > 0):
                # if lower, halve the timestep and try once more
                tr.dt1 = tr.dt1/2.
                print('dt1={}'.format(tr.dt1))
                if tr.dt1 < 1e-10:
                    print('dt too small')
                    break
            else:
                # if higher, continue steps along the primary
                RV_old = RV_new

        if not tr.IsAimZ:
            dz = r_aim[2]-tr.RV_sec[-1, 2]
            print('UB2 OLD = {:.2f}, z_aim - z_curr = {:.4f} m'
                  .format(UB2, dz))
            UB2 = UB2 - dUB2*dz
            print('UB2 NEW = {:.2f}'.format(UB2))
        else:
            break

    return tr


# %%
def optimize_A3B3(tr, rs, geom, UA3, UB3, dUA3, dUB3,
                  E, B, dt, eps_xy=1e-3, eps_z=1e-3):
    ''' get voltages on A3 and B3 plates
    rs - aim point
    '''
    print('\nEb = {}, UA2 = {}'.format(tr.Ebeam, tr.U[0]))
    tr.dt1 = dt
    tr.dt2 = dt
    tmax = 9e-5
    tr.IsAimXY = False
    tr.IsAimZ = False
    RV0 = np.array([tr.RV_sec[0]])

    n_stepsA3 = 0
    while not (tr.IsAimXY and tr.IsAimZ):
        tr.U[2:] = [UA3, UB3]
        tr.pass_sec(RV0, rs, E, B, geom, tmax=tmax,
                    eps_xy=eps_xy, eps_z=eps_z)

        drXY = np.linalg.norm(rs[:2]-tr.RV_sec[-1, :2]) * \
            np.sign(np.cross(tr.RV_sec[-1, :2], rs[:2]))
        print('\n UA3 OLD = {:.2f} kV, dr XY = {:.4f} m'.format(UA3, drXY))
        print('IsAimXY = ', tr.IsAimXY)
        # if drXY < 1e-2:
        #     dUA3 = 10.0

        UA3 = UA3 + dUA3*drXY
        print('UA3 NEW = {:.2f} kV'.format(UA3))
        n_stepsA3 += 1

        if abs(UA3) > 200.:
            print('ALPHA3 failed, voltage too high')
            return tr
        if n_stepsA3 > 100:
            print('ALPHA3 failed, too many steps')
            return tr

        # dz = rs[2] - tr.RV_sec[-1, 2]
        # print('\n UB3 OLD = {:.2f} kV, dZ = {:.4f} m'.format(UB3, dz))
        if abs(drXY) < 0.01:
            n_stepsZ = 0
            while not tr.IsAimZ:
                print('pushing Z direction')
                tr.U[2:] = [UA3, UB3]
                tr.pass_sec(RV0, rs, E, B, geom,
                            eps_xy=eps_xy, eps_z=eps_z, tmax=tmax)
                # tr.IsAimZ = True  # if you want to skip UB3 calculation
                dz = rs[2] - tr.RV_sec[-1, 2]
                print(' UB3 OLD = {:.2f} kV, dZ = {:.4f} m'
                      .format(UB3, dz))
                print('IsAimXY = ', tr.IsAimXY)
                print('IsAimZ = ', tr.IsAimZ)
                # if abs(dz) < 1e-2 and n_stepsZ < 50:
                #     dUB3 = 20.0
                #     print('dUoct = ', dUB3)
                # if n_stepsZ > 100:
                #     dUB3 = 7.0
                #     print('dUoct = ', dUB3)
                #     print('n_stepsZ = ', n_stepsZ)
                # if n_stepsZ > 200:
                #     print('too long!')
                #     break

                UB3 = UB3 - dUB3*dz
                n_stepsZ += 1
                if abs(UB3) > 190.:
                    print('BETA3 failed, voltage too high')
                    return tr
                if n_stepsZ > 100:
                    print('BETA3 failed, too many steps')
                    return tr
                # print('UB3 NEW = {:.2f} kV'.format(UB3))
            n_stepsA3 = 0
            print('n_stepsZ = ', n_stepsZ)
            dz = rs[2] - tr.RV_sec[-1, 2]
            print('UB3 NEW = {:.2f} kV, dZ = {:.4f} m'.format(UB3, dz))

    return tr


# %%
def pass_to_slits(tr, dt, E, B, geom, timestep_divider=10):
    ''' pass trajectories to slits and save secondaries which get into slits
    '''
    tr.dt1 = dt
    tr.dt2 = dt
    # find the number of slits
    n_slits = geom.slits_edges.shape[0]
    tr.add_slits(n_slits)
    # find slits position
    r_slits = geom.slits_edges
    rs = geom.r_dict['slit']
    slit_plane_n = geom.slit_plane_n
    slits_spot = geom.slits_spot

    # pass fan of trajectories
    tr.pass_fan(rs, E, B, geom, stop_plane_n=slit_plane_n,
                r_plasma=0.3, R=0.36, elon=1.8,
                no_intersect=True, no_out_of_bounds=True)
    # create slits polygon
    ax_index = np.argmax(slit_plane_n)
    slits_spot_flat = np.delete(slits_spot, ax_index, 1)
    slits_spot_poly = path.Path(slits_spot_flat)

    # find which secondaries get into slits spot
    # list of sec trajectories indexes which get into slits spot
    sec_ind = []
    for i in range(len(tr.Fan)):
        fan_tr = tr.Fan[i]
        intersect_coords_flat = np.delete(fan_tr[-1, :3], ax_index, 0)
        if slits_spot_poly.contains_point(intersect_coords_flat):
            sec_ind.append(i)

    print('\nStarting precise fan calculation')
    # divide the timestep
    tr.dt1 = dt/timestep_divider
    tr.dt2 = dt
    # number of steps during new fan calculation
    n_steps = timestep_divider * (len(sec_ind) + 1)
    # list for new trajectories
    fan_list = []
    # take the point to start fan calculation
    RV_old = tr.Fan[sec_ind[0]-1][0]
    RV_old = np.array([RV_old])
    RV_new = RV_old

    i_steps = 0
    while i_steps <= n_steps:
        # pass new secondary trajectory
        tr.pass_sec(RV_new, rs, E, B, geom,
                    stop_plane_n=slit_plane_n, tmax=9e-5,
                    eps_xy=1e-3, eps_z=1)
        # make a step on primary trajectory
        r = RV_old[0, :3]
        B_local = return_B(r, B)
        E_local = np.array([0., 0., 0.])
        RV_new = runge_kutt(tr.q, tr.m, RV_old, tr.dt1, E_local, B_local)
        RV_old = RV_new
        i_steps += 1
        if not (tr.IntersectGeometrySec or tr.B_out_of_bounds):
            fan_list.append(tr.RV_sec)
    print('\nPrecise fan calculated')

    # choose secondaries which get into slits
    # start slit cycle
    for i_slit in range(n_slits):
        print('\nslit = {}'.format(i_slit+1))
        print('center of the slit = ', r_slits[i_slit, 0, :], '\n')
        # create slit polygon
        slit_flat = np.delete(r_slits[i_slit, 1:, :], ax_index, 1)
        slit_poly = path.Path(slit_flat)
        zones_list = []  # list for ion zones coordinates
        rv_list = []  # list for RV arrays of secondaries
        for fan_tr in fan_list:
            # get last coordinates of the secondary trajectory
            intersect_coords_flat = np.delete(fan_tr[-1, :3], ax_index, 0)
            if slit_poly.contains_point(intersect_coords_flat):
                print('slit {} ok!\n'.format(i_slit+1))
                rv_list.append(fan_tr)
                zones_list.append(fan_tr[0, :3])

        tr.RV_sec_toslits[i_slit] = rv_list
        tr.ion_zones[i_slit] = np.array(zones_list)
    tr.fan_to_slits = fan_list

    return tr


# %%
def translate(input_array, xyz):
    '''
    move the vector in space
    :param xyz: 3 component vector
    :return: translated input_array
    '''
    if input_array is not None:
        input_array += np.array(xyz)

    return input_array


def rot_mx(axis=(1, 0, 0), deg=0):
    '''
    function calculates rotation matrix
    :return: rotation matrix
    '''
    n = axis
    ca = np.cos(np.radians(deg))
    sa = np.sin(np.radians(deg))
    R = np.array([[n[0]**2*(1-ca)+ca, n[0]*n[1]*(1-ca)-n[2]*sa,
                   n[0]*n[2]*(1-ca)+n[1]*sa],

                  [n[1]*n[0]*(1-ca)+n[2]*sa, n[1]**2*(1-ca)+ca,
                   n[1]*n[2]*(1-ca)-n[0]*sa],

                  [n[2]*n[0]*(1-ca)-n[1]*sa, n[2]*n[1]*(1-ca)+n[0]*sa,
                   n[2]**2*(1-ca)+ca]])
    return R


def rotate(input_array, axis=(1, 0, 0), deg=0):
    '''
    rotate vector around given axis by deg degrees
    :param axis: axis of rotation
    :param deg: angle in degrees
    :return: rotated input_array
    '''
    if input_array is not None:
        R = rot_mx(axis, deg)
        input_array = np.dot(input_array, R.T)
    return input_array


# %% Intersection check functions
def line_plane_intersect(planeNormal, planePoint, rayDirection,
                         rayPoint, eps=1e-6):
    ''' function returns intersection point between plane and ray
    '''
    ndotu = np.dot(planeNormal, rayDirection)
    if abs(ndotu) < eps:
        # print('no intersection or line is within plane')
        return np.full_like(planeNormal, np.nan)
    else:
        w = rayPoint - planePoint
        si = -np.dot(planeNormal, w) / ndotu
        Psi = w + si * rayDirection + planePoint
        return Psi


def is_between(A, B, C, eps=1e-6):
    ''' function returns True if point C is on the segment AB
    (between A and B)'''
    if np.isnan(C).any():
        return False
    # check if the points are on the same line
    crossprod = np.cross(B-A, C-A)
    if np.linalg.norm(crossprod) > eps:
        return False
    # check if the point is between
    dotprod = np.dot(B-A, C-A)
    if dotprod < 0 or dotprod > np.linalg.norm(B-A)**2:
        return False
    return True


def segm_intersect(A, B, C, D):  # doesn't work with collinear case
    ''' function returns true if line segments AB and CD intersect
    '''
    def order(A, B, C):
        ''' If counterclockwise return True
            If clockwise return False '''
        return (C[1]-A[1]) * (B[0]-A[0]) > (B[1]-A[1]) * (C[0]-A[0])

    # Return true if line segments AB and CD intersect
    return order(A, C, D) != order(B, C, D) and \
        order(A, B, C) != order(A, B, D)


def segm_poly_intersect(polygon_coords, segment_coords):
    ''' check segment and polygon intersection'''
    polygon_normal = np.cross(polygon_coords[2, 0:3]-polygon_coords[0, 0:3],
                              polygon_coords[1, 0:3]-polygon_coords[0, 0:3])
    polygon_normal = polygon_normal/np.linalg.norm(polygon_normal)

    intersect_coords = line_plane_intersect(polygon_normal,
                                            polygon_coords[2, 0:3],
                                            segment_coords[1, 0:3] -
                                            segment_coords[0, 0:3],
                                            segment_coords[0, 0:3])
    if np.isnan(intersect_coords).any():
        return False
    else:
        i = np.argmax(polygon_normal)
        polygon_coords_flat = np.delete(polygon_coords, i, 1)
        intersect_coords_flat = np.delete(intersect_coords, i, 0)
        p = path.Path(polygon_coords_flat)
        return p.contains_point(intersect_coords_flat) and \
            is_between(segment_coords[0, 0:3], segment_coords[1, 0:3],
                       intersect_coords)

    # check projections on XY and XZ planes
    # pXY = path.Path(polygon_coords[:, [0, 1]])  # XY plane
    # pXZ = path.Path(polygon_coords[:, [0, 2]])  # XZ plane
    # return pXY.contains_point(intersect_coords[[0, 1]]) and \
    #     pXZ.contains_point(intersect_coords[[0, 2]]) and \
    #         is_between(segment_coords[0, 0:3], segment_coords[1, 0:3],
    #                   intersect_coords)


# %%
def plate_flags(range_x, range_y, range_z, U, plts_geom, plts_angles):

    alpha, beta, gamma = plts_angles
    length, width, thick, gap = plts_geom

    # Geometry rotated in system based on central point between plates
    # upper plate
    UP1 = np.array([-length/2., gap/2. + thick, width/2.])
    UP2 = np.array([-length/2., gap/2. + thick, -width/2.])
    UP3 = np.array([length/2., gap/2. + thick, -width/2.])
    UP4 = np.array([length/2., gap/2. + thick, width/2.])
    UP5 = np.array([-length/2., gap/2., width/2.])
    UP6 = np.array([-length/2., gap/2., -width/2.])
    UP7 = np.array([length/2., gap/2., -width/2.])
    UP8 = np.array([length/2., gap/2., width/2.])
    UP_full = np.array([UP1, UP2, UP3, UP4, UP5, UP6, UP7, UP8])
    UP_rotated = UP_full.copy()
    for i in range(UP_full.shape[0]):
        UP_rotated[i, :] = rotate(UP_rotated[i, :],
                                  axis=(1, 0, 0), deg=gamma)
        UP_rotated[i, :] = rotate(UP_rotated[i, :],
                                  axis=(0, 0, 1), deg=alpha)
        UP_rotated[i, :] = rotate(UP_rotated[i, :],
                                  axis=(0, 1, 0), deg=beta)

    # lower plate
    LP1 = np.array([-length/2., -gap/2. - thick, width/2.])
    LP2 = np.array([-length/2., -gap/2. - thick, -width/2.])
    LP3 = np.array([length/2., -gap/2. - thick, -width/2.])
    LP4 = np.array([length/2., -gap/2. - thick, width/2.])
    LP5 = np.array([-length/2., -gap/2., width/2.])
    LP6 = np.array([-length/2., -gap/2., -width/2.])
    LP7 = np.array([length/2., -gap/2., -width/2.])
    LP8 = np.array([length/2., -gap/2., width/2.])
    LP_full = np.array([LP1, LP2, LP3, LP4, LP5, LP6, LP7, LP8])
    LP_rotated = LP_full.copy()
    for i in range(LP_full.shape[0]):
        LP_rotated[i, :] = rotate(LP_rotated[i, :],
                                  axis=(1, 0, 0), deg=gamma)
        LP_rotated[i, :] = rotate(LP_rotated[i, :],
                                  axis=(0, 0, 1), deg=alpha)
        LP_rotated[i, :] = rotate(LP_rotated[i, :],
                                  axis=(0, 1, 0), deg=beta)

    # Find coords of 'cubes' containing each plate
    upper_cube = np.array([np.min(UP_rotated, axis=0),
                           np.max(UP_rotated, axis=0)])
    lower_cube = np.array([np.min(LP_rotated, axis=0),
                           np.max(LP_rotated, axis=0)])

    # create mask for plates
    upper_plate_flag = np.full_like(U, False, dtype=bool)
    lower_plate_flag = np.full_like(U, False, dtype=bool)
    for i in range(range_x.shape[0]):
        for j in range(range_y.shape[0]):
            for k in range(range_z.shape[0]):
                x = range_x[i]
                y = range_y[j]
                z = range_z[k]
                # check upper cube
                if (x >= upper_cube[0, 0]) and (x <= upper_cube[1, 0]) and \
                   (y >= upper_cube[0, 1]) and (y <= upper_cube[1, 1]) and \
                   (z >= upper_cube[0, 2]) and (z <= upper_cube[1, 2]):
                    r_rot = rotate(rotate(rotate(
                        np.array([x, y, z]), axis=(0, 1, 0), deg=-beta),
                        axis=(0, 0, 1), deg=-alpha),
                        axis=(1, 0, 0), deg=-gamma)
                    # define masks for upper and lower plates
                    upper_plate_flag[i, j, k] = (r_rot[0] >= -length/2.) and \
                        (r_rot[0] <= length/2.) and (r_rot[2] >= -width/2.) and \
                        (r_rot[2] <= width/2.) and (r_rot[1] >= gap/2.) and \
                        (r_rot[1] <= gap/2. + thick)
                # check lower cube
                if (x >= lower_cube[0, 0]) and (x <= lower_cube[1, 0]) and \
                   (y >= lower_cube[0, 1]) and (y <= lower_cube[1, 1]) and \
                   (z >= lower_cube[0, 2]) and (z <= lower_cube[1, 2]):
                    r_rot = rotate(rotate(rotate(
                        np.array([x, y, z]), axis=(0, 1, 0), deg=-beta),
                        axis=(0, 0, 1), deg=-alpha),
                        axis=(1, 0, 0), deg=-gamma)
                    # define masks for upper and lower plates
                    lower_plate_flag[i, j, k] = (r_rot[0] >= -length/2.) and \
                        (r_rot[0] <= length/2.) and (r_rot[2] >= -width/2.) and \
                        (r_rot[2] <= width/2.) and \
                        (r_rot[1] >= -gap/2. - thick) and \
                        (r_rot[1] <= -gap/2.)

    return UP_rotated, LP_rotated, upper_plate_flag, lower_plate_flag


def return_E(r, Ein, U):
    '''
    take dot and try to interpolate electiric fields in it
    return: interpolated Electric field
    :param Ein: list of interpolants for Ex, Ey, Ez
    '''
    Eout = np.zeros(3)
    for i in range(len(Ein)):
        try:
            Eout[0] += Ein[i][0](r)*U[i]
            Eout[1] += Ein[i][1](r)*U[i]
            Eout[2] += Ein[i][2](r)*U[i]
            # print('U = ', U[i])
            # print('E = ', Eout)
        except:
            continue
    return Eout


def return_B(r, Bin):
    Bx_interp, By_interp, Bz_interp = Bin[0], Bin[1], Bin[2]
    Bout = np.c_[Bx_interp(r), By_interp(r), Bz_interp(r)]
    return Bout


def save_E(beamline, plts_name, Ex, Ey, Ez, angles, geom,
           domain, plate1, plate2, dirname='elecfield'):
    '''
    save Ex, Ey, Ez arrays to file
    '''
    dirname = dirname + '/' + beamline + \
        '_alpha_{}_beta_{}'.format(int(angles[0]), int(angles[1]))

    if not os.path.exists(dirname):
        try:
            os.makedirs(dirname, 0o700)
            print("Directory ", dirname, " created ")
        except OSError as e:
            if e.errno != errno.EEXIST:
                raise

    fname = plts_name + '_geometry.dat'
    # erases data from file before writing
    open(dirname + '/' + fname, 'w').close()
    with open(dirname + '/' + fname, 'w') as myfile:
        myfile.write('{} {} {} {} # plate\'s length, thic, width and gap'
                     '\n'.format(geom[0], geom[1], geom[2], geom[3]))
        myfile.write('{} {} {} # plate\'s alpha, beta and gamma angle'
                     '\n'.format(angles[0], angles[1], angles[2]))
        myfile.write('{} {} {} # number of dots (x,y,z)'
                     '\n'.format(Ex.shape[0], Ex.shape[1], Ex.shape[2]))
        myfile.write('{} {} {} {} # border x, border y, border z, delta'
                     '\n'.format(domain[0], domain[1], domain[2], domain[3]))
        for i in range(plate1.shape[0]):
            myfile.write(np.array2string(plate1[i], precision=4)[1:-1] +
                         ' # 1st plate rotated \n')
        for i in range(plate2.shape[0]):
            myfile.write(np.array2string(plate2[i], precision=4)[1:-1] +
                         ' # 2nd plate rotated \n')

    np.save(dirname + '/' + plts_name + '_Ex', Ex)
    np.save(dirname + '/' + plts_name + '_Ey', Ey)
    np.save(dirname + '/' + plts_name + '_Ez', Ez)

    print('Electric field saved, ' + plts_name + '\n')


def read_E(beamline, geom, dirname='elecfield'):
    '''
    read plate's shape and angle parametres along with electric field values
    from provided file (should lbe in the same directory)
    return: intrepolator function for electric field
    :param fname: filename
    :param r_dict: dict (key: name of plate,
                         value: coordinate array of plates centre)
    '''
    r_dict = geom.r_dict
    if beamline == 'prim':
        plts_angles = geom.prim_angles
    elif beamline == 'sec':
        plts_angles = geom.sec_angles

    E = []
    edges_dict = {}
    dirname = dirname + '/' + beamline + \
        '_alpha_{}_beta_{}'.format(int(plts_angles[0]),
                                   int(plts_angles[1]))

    for filename in os.listdir(dirname):
        if filename.endswith('.dat'):
            plts_name = filename[0:2]
            r_new = r_dict[plts_name]

            edges_list = []
            print('\n Reading geometry {} ...'.format(plts_name))

            with open(dirname + '/' + filename, 'r') as f:
                geom = [float(i) for i in f.readline().split()[0:4]]
                angles = [float(i) for i in f.readline().split()[0:3]]
                dim = [int(i) for i in f.readline().split()[0:3]]
                domain = [float(i) for i in f.readline().split()[0:4]]
                for line in f:
                    edges_list.append([float(i) for i in line.split()[0:3]])

            edges_list = np.array(edges_list)
            edges_dict[plts_name] = np.array([edges_list[0:4, :] + r_new,
                                              edges_list[4:, :] + r_new])
            print('position', r_new)
            Ex = np.load(dirname + '/' + plts_name + '_Ex.npy')
            Ey = np.load(dirname + '/' + plts_name + '_Ey.npy')
            Ez = np.load(dirname + '/' + plts_name + '_Ez.npy')

            x = np.arange(-domain[0]/2., domain[0]/2., domain[3])
            y = np.arange(-domain[1]/2., domain[1]/2., domain[3])
            z = np.arange(-domain[2]/2., domain[2]/2., domain[3])

            # make interpolation for Ex, Ey, Ez
            Ex_interp = RegularGridInterpolator((x + r_new[0], y + r_new[1],
                                                 z + r_new[2]), Ex)
            Ey_interp = RegularGridInterpolator((x + r_new[0], y + r_new[1],
                                                 z + r_new[2]), Ey)
            Ez_interp = RegularGridInterpolator((x + r_new[0], y + r_new[1],
                                                 z + r_new[2]), Ez)
            E_read = [Ex_interp, Ey_interp, Ez_interp]

            E.append(E_read)

    return E, edges_dict


def read_B(Btor, Ipl, PF_dict, dirname='magfield'):
    '''
    read magnetic field values from provided file (should be in
                                                   the same directory)
    return: list of intrepolator functions for Bx, By, Bz
    :param dirname: name of directory with magfield dats
    '''
    print('\n Reading Magnetic field')
    B_dict = {}
    for filename in os.listdir(dirname):
        if filename.endswith('.dat'):  # load geometry .dat file
            with open(dirname + '/' + filename, 'r') as f:
                volume_corner1 = [float(i) for i in f.readline().split()[0:3]]
                volume_corner2 = [float(i) for i in f.readline().split()[0:3]]
                resolution = float(f.readline().split()[0])
            continue
        elif 'Tor' in filename:
            print('Reading toroidal magnetic field...')
            B_read = np.load(dirname + '/' + filename) * Btor
            name = 'Tor'

        elif 'Plasm' in filename:
            print('Reading plasma field...')
            B_read = np.load(dirname + '/' + filename)  #* Ipl
            name = 'Plasm'

        else:
            name = filename.replace('magfield', '').replace('.npy', '')
            print('Reading {} magnetic field...'.format(name))
            # Icir = PF_dict[name]
            Icir = PF_dict[name][4] * PF_dict[name][5]
            print('Current = ', Icir)
            B_read = np.load(dirname + '/' + filename) * Icir

        B_dict[name] = B_read

    # create grid of points
    grid = np.mgrid[volume_corner1[0]:volume_corner2[0]:resolution,
                    volume_corner1[1]:volume_corner2[1]:resolution,
                    volume_corner1[2]:volume_corner2[2]:resolution]

    B = np.zeros_like(B_read)
    for key in B_dict.keys():
        B += B_dict[key]

#    cutoff = 10.0
#    Babs = np.linalg.norm(B, axis=1)
#    B[Babs > cutoff] = [np.nan, np.nan, np.nan]
    # make an interpolation of B
    x = np.arange(volume_corner1[0], volume_corner2[0], resolution)
    y = np.arange(volume_corner1[1], volume_corner2[1], resolution)
    z = np.arange(volume_corner1[2], volume_corner2[2], resolution)
    Bx = B[:, 0].reshape(grid.shape[1:])
    By = B[:, 1].reshape(grid.shape[1:])
    Bz = B[:, 2].reshape(grid.shape[1:])
    Bx_interp = RegularGridInterpolator((x, y, z), Bx)
    By_interp = RegularGridInterpolator((x, y, z), By)
    Bz_interp = RegularGridInterpolator((x, y, z), Bz)
    print('Interpolants for magnetic field created')

    hbplot.plot_B_stream(B, volume_corner1, volume_corner2, resolution, grid,
                         plot_sep=True)

    B_list = [Bx_interp, By_interp, Bz_interp]

    return B_list


# %% poloidal field coils
def import_PFcoils(filename):
    ''' import a dictionary with poloidal field coils parameters
    {'NAME': (x center, y center, width along x, width along y [m],
               current [MA-turn], N turns)}
    Andreev, VANT 2014, No.3
    '''
    d = {}  # defaultdict(list)
    with open(filename) as f:
        for line in f:
            if line[0] == '#':
                continue
            lineList = line.split(', ')
            key, val = lineList[0], tuple([float(i) for i in lineList[1:]])
            d[key] = val
    return d


def import_PFcur(filename, pf_coils):
    '''
    Creates dictionary with coils names and currents from TOKAMEQ file
    :param filename: Tokameqs filename
    :param coils: coil dict (we only take keys)
    :return: PF dictianary with currents
    '''
    with open(filename, 'r') as f:
        data = f.readlines()  # read tokameq file
    PF_dict = {}  # Here we will store coils names and currents
    pf_names = list(pf_coils)  # get coils names
    n_coil = 0  # will be used for getting correct coil name
    for i in range(len(data)):
        if data[i].strip() == 'External currents:':
            n_line = i + 2  # skip 2 lines and read from the third
            break
    while float(data[n_line].strip().split()[3]) != 0:
        key = pf_names[n_coil]
        val = data[n_line].strip().split()[3]
        PF_dict[key] = float(val)
        n_line += 1
        n_coil += 1

    return PF_dict


# %%
def import_Bflux(filename):
    ''' import magnetic flux from Tokameq file'''
    with open(filename, 'r') as f:
        data = f.readlines()

    # R coordinate corresponds to X, Z coordinate corresponds to Y
    NrNz = []
    for i in data[2].strip().split():
        if i.isdigit():
            NrNz.append(i)
    Nx = int(NrNz[0]) + 1
    Ny = int(NrNz[1]) + 1

    for i in range(len(data)):
        if ' '.join(data[i].strip().split()[:4]) == 'Flux at the boundary':
            bound_flux = float(data[i].strip().split()[-1])
        if data[i].strip() == 'Poloidal flux F(r,z)':
            index = i

    x_vals = [float(r) for r in data[index+1].strip().split()[1:]]
    x_vals = np.array(x_vals)

    Psi_data = [i.strip().split() for i in data[index+2:index+2+Ny]]
    Psi_vals = []
    y_vals = []
    for line in Psi_data:
        y_vals.append(float(line[0]))
        Psi_vals.append([float(j) for j in line[1:]])

    y_vals = np.array(y_vals)
    Psi_vals = np.array(Psi_vals)
    return Psi_vals, x_vals, y_vals, bound_flux


# %%
def save_traj_list(traj_list, Btor, Ipl, r_aim, dirname='output'):
    '''
    Save list of Traj objects to pickle file
    :param traj_list: list of trajectories
    '''

    if len(traj_list) == 0:
        print('traj_list empty!')
        return

    Ebeam_list = []
    UA2_list = []

    for traj in traj_list:
        Ebeam_list.append(traj.Ebeam)
        UA2_list.append(traj.U[0])

    dirname = dirname + '/' + 'B{}_I{}'.format(int(Btor), int(Ipl))

    if not os.path.exists(dirname):
        try:
            os.makedirs(dirname, 0o700)
            print("Directory ", dirname, " created ")
        except OSError as e:
            if e.errno != errno.EEXIST:
                raise

    fname = dirname + '/' + 'E{}-{}'.format(int(min(Ebeam_list)),
                                            int(max(Ebeam_list))) + \
        '_UA2{}-{}'.format(int(min(UA2_list)), int(max(UA2_list))) + \
        '_alpha{}_beta{}'.format(int(round(traj.alpha)),
                                 int(round(traj.beta))) +\
        '_x{}y{}z{}.pkl'.format(int(r_aim[0]*100), int(r_aim[1]*100),
                                int(r_aim[2]*100))

    with open(fname, 'wb') as f:
        pc.dump(traj_list, f, -1)

    print('\nSAVED LIST: \n' + fname)


# %%
def read_traj_list(fname, dirname='output'):
    '''
    import list of Traj objects from .pkl file
    '''
    with open(dirname + '/'+fname, 'rb') as f:
        traj_list = pc.load(f)
    return traj_list


# %%
def save_png(fig, name, save_dir='Results/Grids'):
    """
    Saves picture as name.png
    Args:
    :fig - array of figures to save
    :name - array of picture names
    :save_dir - directory used to store results
    """

    # check wether directory exist and if not - create one
    if not os.path.isdir(save_dir):
        os.makedirs(save_dir)
        print('LOG: {} directory created'.format(save_dir))
    print('LOG: Saving pictures to {}'.format(save_dir+'/'))
    for fig, name in zip(fig, name):
        # save fig with tight layout
        fig_savename = str(name + '.png')
        fig.savefig(save_dir + '/' + fig_savename, bbox_inches='tight')
        print('LOG: Figure ' + fig_savename + ' saved')
