import numpy as np
import hibplib as hb
import hibpplotlib as hbplot
import copy

# %%
''' MAIN '''

if __name__ == '__main__':
    # timestep [sec]
    dt = 0.6e-7  # 0.4e-7

    # toroidal field on the axis
    Btor = 0.7  # [T]
    Ipl = 0.5  # Plasma current [MA]
    q = 1.60217662e-19  # electron charge [Co]
    # m_ion = 204.3833 * 1.6605e-27  # Tl ion mass [kg]
    m_ion = 133.0 * 1.6605e-27  # Cs ion mass [kg]

    # initial beam energy range
    dEbeam = 5.
    Ebeam_range = np.arange(40., 40. + dEbeam, dEbeam)  # [keV]

    # A2 plates voltage
    dUA2 = 1.0  # [kV]
    UA2_range = np.arange(-7., 7. + dUA2, dUA2)  # [kV]

    # B2 plates voltage
    UB2 = 0.0  # [kV]
    dUB2 = 5.0  # [kV/m]

    # B3 voltages
    UB3 = -1.0  # [kV]
    dUB3 = 40.0  # [kV/m]

    # A3 voltages
    UA3 = 0.0  # [kV]
    dUA3 = 40.0  # [kV/m]

# %% PRIMARY beamline geometry
    geomGlob = hb.Geometry()

    # plasma parameters
    geomGlob.R = 0.36  # tokamak major radius [m]
    geomGlob.r_plasma = 0.3  # plasma minor radius [m]
    geomGlob.elon = 1.8  # plasma elongation

    # alpha and beta angles of the PRIMARY beamline
    alpha_prim = 60.  # deg
    beta_prim = -5  # deg
    gamma_prim = 0.  # deg
    geomGlob.prim_angles = np.array([alpha_prim, beta_prim, gamma_prim])

    # coordinates of the injection port [m]
    # coordinates of injection port and chamber entrance coordinates [m]
#   port 75 deg
    xpatr = 0.4672459
    ypatr = 0.5753405
    zpatr = 0.0
    geomGlob.chamb_ent = [(0.425, 0.583), (0.4, 0.497),
                          (0.501, 0.491), (0.486, 0.433)]
    geomGlob.chamb_ext = [(0.63, 0.226), (0.63, 0.085),
                          (0.63, -0.085), (0.63, -0.226)]

#   port 25 deg
#    xpatr = 0.4659326
#    ypatr = 0.5678362
#    zpatr = 0.0
#    geomGlob.chamb_ent = [(0.428, 0.588), (0.31, 0.533),
#                          (0.458, 0.525), (0.399, 0.498)]
#    geomGlob.chamb_ext = [(0.63, 0.226), (0.63, 0.145),
#                          (0.63, -0.145), (0.63, -0.226)]

#   port 90 deg
#    xpatr = 0.451501 + 0.02
#    ypatr = 0.5720209
#    zpatr = 0.0
#    geomGlob.chamb_ent = [(0.397, 0.55), (0.397, 0.499),
#                          (0.50, 0.55), (0.5, 0.415)]
#    geomGlob.chamb_ext = [(0.63, 0.226), (0.63, 0.085),
#                          (0.63, -0.085), (0.63, -0.226)]

    geomGlob.r_dict['patr'] = np.array([xpatr, ypatr, zpatr])

    # distance from the injection port to the Alpha2 plates
    dist_A2 = 0.15/2 + 0.05  # [m]
    # distance from the injection port to the Beta2 plates
    dist_B2 = dist_A2 + 0.15 + 0.05  # [m]
    # distance from the injection port to the initial piont of the traj [m]
    dist_0 = dist_B2 + 0.15/2 + 0.1

    # convert degrees to radians
    drad = np.pi/180
    # coordinates of the center of the ALPHA2 plates
    xA2 = xpatr + dist_A2*np.cos(alpha_prim*drad) * \
          np.cos(beta_prim*drad)
    yA2 = ypatr + dist_A2*np.sin(alpha_prim*drad)
    zA2 = zpatr - dist_A2*np.cos(alpha_prim*drad) * \
          np.sin(beta_prim*drad)
    rA2 = np.array([xA2, yA2, zA2])
    geomGlob.r_dict['A2'] = rA2

    # coordinates of the center of the BETA2 plates
    xB2 = xpatr + dist_B2*np.cos(alpha_prim*drad) * \
          np.cos(beta_prim*drad)
    yB2 = ypatr + dist_B2*np.sin(alpha_prim*drad)
    zB2 = zpatr - dist_B2*np.cos(alpha_prim*drad) * \
          np.sin(beta_prim*drad)
    rB2 = np.array([xB2, yB2, zB2])
    geomGlob.r_dict['B2'] = rB2

    # coordinates of the initial point of the trajectory [m]
    x0 = xpatr + dist_0*np.cos(alpha_prim*drad) * \
         np.cos(beta_prim*drad)
    y0 = ypatr + dist_0*np.sin(alpha_prim*drad)
    z0 = zpatr - dist_0*np.cos(alpha_prim*drad) * \
         np.sin(beta_prim*drad)
    r0 = np.array([x0, y0, z0])
    geomGlob.r_dict['r0'] = r0

# %% AIM position (BEFORE the Secondary beamline)
    xaim = 0.75
    yaim = 0.0
    zaim = 0.0
    r_aim = np.array([xaim, yaim, zaim])
    geomGlob.r_dict['aim'] = r_aim

    # angles of aim plane normal
    alpha_aim = 0.
    beta_aim = 0.
    stop_plane_n = np.array([np.cos(alpha_aim*drad)*np.cos(beta_aim*drad),
                             np.sin(alpha_aim*drad),
                             np.cos(alpha_aim*drad)*np.sin(beta_aim*drad)])
    stop_plane_n = stop_plane_n/np.linalg.norm(stop_plane_n)

# %% SECONDARY beamline geometry
    # alpha and beta angles of the SECONDARY beamline
    alpha_sec = -15.  # deg
    beta_sec = 20.  # deg
    gamma_sec = 0.  # deg
    geomGlob.sec_angles = np.array([alpha_sec, beta_sec, gamma_sec])

    # distance from r_aim to the ALPHA3 center
    dist_A3 = 0.05  # 1/2 of plates length
    # distance from r_aim to the BETA3 center
    dist_B3 = dist_A3 + 0.15 +0.05
    # distance from r_aim the entrance slit of the analyzer
    dist_s = dist_B3 + 0.15/2 + 0.1

    # coordinates of the center of the ALPHA3 plates
    xA3 = xaim + dist_A3*np.cos(alpha_sec*drad) * \
        np.cos(beta_sec*drad)
    yA3 = yaim + dist_A3*np.sin(alpha_sec*drad)
    zA3 = zaim - dist_A3*np.cos(alpha_sec*drad) * \
        np.sin(beta_sec*drad)
    rA3 = np.array([xA3, yA3, zA3])
    geomGlob.r_dict['A3'] = rA3

    # coordinates of the center of the BETA3 plates
    xB3 = xaim + dist_B3*np.cos(alpha_sec*drad) * \
        np.cos(beta_sec*drad)
    yB3 = yaim + dist_B3*np.sin(alpha_sec*drad)
    zB3 = zaim - dist_B3*np.cos(alpha_sec*drad) * \
        np.sin(beta_sec*drad)
    rB3 = np.array([xB3, yB3, zB3])
    geomGlob.r_dict['B3'] = rB3

    # Coordinates of the CENTRAL slit
    xs = xaim + dist_s*np.cos(alpha_sec*drad) * \
        np.cos(beta_sec*drad)
    ys = yaim + dist_s*np.sin(alpha_sec*drad)
    zs = zaim - dist_s*np.cos(alpha_sec*drad) * \
        np.sin(beta_sec*drad)
    rs = np.array([xs, ys, zs])
    geomGlob.r_dict['slit'] = rs

# %% print info
    print('\nShot parameters: Btor = {} T, Ipl = {} MA'. format(Btor, Ipl))
    print('Primary beamline angles: ', geomGlob.prim_angles[0:2])
    print('r0 = ', np.round(r0, 3))
    print('r_aim = ', r_aim)
    print('r_slit = ', np.round(rs, 3))

# %% GEOMETRY
    # Toroidal Field coil
    geomGlob.coil = np.loadtxt('TFCoil.dat')  # [m]
    # Poloidal Field coils
    geomGlob.pf_coils = hb.import_PFcoils('PFCoils.dat')
    # Camera contour
    geomGlob.camera = np.loadtxt('port75.txt', skiprows=1) / 1000
    # Separatrix contour
    geomGlob.sep = np.loadtxt('Globus_sep.txt')
    # First wall innner and outer contours
    # geomGlob.in_fw = np.loadtxt('infw.txt') / 1000  # [m]
    # geomGlob.out_fw = np.loadtxt('outfw.txt') / 1000  # [m]

# %% ELECTRIC field
    ''' Electric field part '''
    # load E for primary beamline
    E_prim, edges_prim = hb.read_E('prim', geomGlob)
    geomGlob.plates_edges.update(edges_prim)
    print('Primary Beamline loaded')

    # load E for secondary beamline
    try:
        E_sec, edges_sec = hb.read_E('sec', geomGlob)
        geomGlob.plates_edges.update(edges_sec)
    except FileNotFoundError:
        print('Secondary Beamline NOT FOUND')
        E_sec = []

    E = E_prim + E_sec

# %% MAGNETIC field
    ''' Magnetic field part '''
    pf_coils = hb.import_PFcoils('PFCoils.dat')

    # PF_dict = hb.import_PFcur('{}MA_sn.txt'.format(int(abs(Ipl))), pf_coils)
    PF_dict = pf_coils
    if 'B' not in locals():
        dirname = 'magfield'
        B = hb.read_B(Btor, Ipl, PF_dict, dirname)

# %% PRIMARY beamline optimization
    print('\n Primary beamline optimization')
    # define list of trajectores that hit r_aim
    traj_list = []

    for Ebeam in Ebeam_range:
        for UA2 in UA2_range:
            print('\n\nE = {} keV; UA2 = {} kV\n'.format(Ebeam, UA2))
            # list of starting voltages
            U_list = [UA2, UB2, UA3, UB3]

            # create new trajectory object
            tr = hb.Traj(q, m_ion, Ebeam, r0, alpha_prim, beta_prim,
                         U_list, dt)

            tr = hb.optimize_B2(tr, r_aim, geomGlob, UB2, dUB2, E, B, dt,
                                stop_plane_n, eps_xy=1e-3, eps_z=1e-3)

            if tr.IntersectGeometry:
                print('NOT saved, primary intersected geometry')
                continue
            if tr.IsAimXY and tr.IsAimZ:
                traj_list.append(tr)
                print('\n Trajectory saved, UB2={:.2f} kV'.format(tr.U[1]))
                UB2 = tr.U[1]
            else:
                print('NOT saved, sth wrong')

# %%
    traj_list_passed = copy.deepcopy(traj_list)

# %% Additonal plots

    hbplot.plot_grid(traj_list_passed, geomGlob, Btor, Ipl, marker_E='')
    hbplot.plot_fan(traj_list_passed, geomGlob, 50., UA2, Btor, Ipl,
                    plot_slits=True, plot_traj=True, plot_all=True)

    hbplot.plot_scan(traj_list_passed, geomGlob, 30., Btor, Ipl)
    # hbplot.plot_scan(traj_list_passed, geomGlob, 120., Btor, Ipl)
    # hbplot.plot_sec_angles(traj_list_passed, Btor, Ipl, Ebeam='all')
    # hbplot.plot_fan(traj_list_passed, geomGlob, 240., 40., Btor, Ipl)

# %% SECONDARY beamline optimization
    print('\n Secondary beamline optimization')
    traj_list_oct = []
    for tr in copy.deepcopy(traj_list_passed):
        tr = hb.optimize_A3B3(tr, rs, geomGlob,
                              UA3, UB3, dUA3, dUB3, E, B, dt,
                              eps_xy=1e-3, eps_z=1e-3)
        if not tr.IntersectGeometrySec:
            traj_list_oct.append(tr)
            print('\n Trajectory saved')
            UA3 = tr.U[2]
            UB3 = tr.U[3]

# %%
    hbplot.plot_traj(traj_list_oct, geomGlob, 40., 0., Btor, Ipl)
    hbplot.plot_scan(traj_list_oct, geomGlob, 40., Btor, Ipl)

# %% Save list of trajectories

    # hb.save_traj_list(traj_list_passed, Btor, Ipl, r_aim)
