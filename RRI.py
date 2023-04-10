import numpy as np
from RRI_Sub import *
from RRI_Dam import *
from RRI_GW import *

def rri():
    """
    RRI model as a function

    :Author: T.Sayama
    :Translator: Y. Chemin
    :Version: 1.4.2.7
    """
    ##########################################################
    ### STEP 0: FILE NAME AND PARAMETER SETTING
    ##########################################################
    # This is now "from RRI_Input import *" at the beginning of the file
    version = "RRI_Input_Format_Ver1_4_2"
    datadir = "/home/yann/dev/RRI_PY/solo30s/"

    rainfile = "rain/rain_.dat"
    demfile = "topo/adem.txt"
    accfile = "topo/acc.txt"
    dirfile = "topo/adir.txt"

    utm = 0				# utm(1) or latlon(0)
    eight_dir = 1			# 4-direction (0), 8-direction(1)
    lasth = 360			        # lasth(hour)
    dt = 600			        # dt(second)
    dt_riv = 60			        # dt_riv
    outnum = 60			        # outnum [-]"
    xllcorner_rain = 110.400000	        # xllcorner_rain
    yllcorner_rain = -8.158333	        # yllcorner_rain
    cellsize_rain_x = 0.0083333333	# cellsize_rain"
    cellsize_rain_y = 0.0083333333	# cellsize_rain"

    ns_river            = 3.000e-2      # ns_river
    num_of_landuse 	= 1		# num_of_landuse
    dif 		= [1]		# diffusion(1) or kinematic(0), one per LU
    ns_slope 	        = [4.000e-1]    # ns_slope, one per LU
    soildepth 	        = [1.0000] 	# soildepth, one per LU
    gammaa 		= [4.750e-1]    # gammaa, one per LU

    ksv 		= [0.0000]     	# kv (m/s), one per LU
    faif 		= [3.163e-1]   	# Sf (m), one per LU

    ka 		        = [0.0000]     	# ka, one per LU
    gammam 		= [0.0000]	# gammam, one per LU
    beta 		= [8.0000]     	# beta, one per LU

    ksg		        = [0.0]         # ksg (m/s), 1/LU -- set zero for no bedrock gw
    gammag		= [0.0370]      # gammag (-), one per LU
    kg0		        = [5.7e-5]      # kg0 (m/s), one per LU
    fpg		        = [0.10]        # fg (-), one per LU
    rgl		        = [0.0]         # rgl, one per LU

    riv_thresh = 20			# riv_thresh
    width_param_c = 5.000		# width_param_c (2.5)
    width_param_s = 3.50e-1		# width_param_s (0.4)
    depth_param_c = 9.50e-1		# depth_param_c (0.1)
    depth_param_s = 2.00e-1		# depth_param_s (0.4)
    height_param = 0.000		# height_param
    height_limit_param = 20		# height_limit_param

    rivfile_switch = 1
    widthfile = "riv/width.txt"
    depthfile = "riv/depth.txt"
    heightfile = "riv/height.txt"

    init_slo_switch = 0
    init_riv_switch = 0
    init_gw_switch = 0
    init_gampt_ff_switch = 0
    initfile_slo = "init/hs_init.out"
    initfile_riv = "init/hr_init.out"
    initfile_gw = "init/hg_init.out"
    initfile_gampt_ff = "init/gamptff_init.out"

    bound_slo_wlev_switch = 0
    bound_riv_wlev_switch = 0
    boundfile_slo_wlev = "bound/hs_wlev_bound.txt"
    boundfile_riv_wlev = "bound/hr_wlev_bound.txt"

    bound_slo_disc_switch = 0
    bound_riv_disc_switch = 0
    boundfile_slo_disc = "bound/qs_bound.txt"
    boundfile_riv_disc = "bound/qr_bound.txt"

    land_switch = 0
    landfile = "topo/landuse.txt"

    dam_switch = 0
    damfile = "damcnt.txt"

    div_switch = 0
    div_switch = "div.txt"

    evp_switch = 0
    evpfile = "rain/Evp.dat"
    xllcorner_evp = 100.000000		# xllcorner_evp
    yllcorner_evp = 10.000000		# yllcorner_evp
    cellsize_evp_x = 0.0083333300	# cellsize x
    cellsize_evp_y = 0.0083333300       # cellsize y

    sec_length_switch = 0
    sec_length_file = "riv/length.txt"

    sec_switch = 0
    sec_map_file = "riv/sec_map.txt"
    sec_file = "riv/section/sec_"

    outswitch_hs = 1
    outswitch_hr = 1
    outswitch_hg = 0
    outswitch_qr = 1
    outswitch_qu = 0
    outswitch_qv = 0
    outswitch_gu = 0
    outswitch_gv = 0
    outswitch_gampt_ff = 1
    outswitch_storage = 1
    outfile_hs = "out/hs_"
    outfile_hr = "out/hr_"
    outfile_hg = "out/hg_"
    outfile_qr = "out/qr_"
    outfile_qu = "out/qu_"
    outfile_qv = "out/qv_"
    outfile_gu = "out/gu_"
    outfile_gv = "out/gv_"
    outfile_gampt_ff = "out/gampt_ff_"
    outfile_storage = "out/storage.dat"

    hydro_switch = 1
    hydro_file = "hydro.txt"
    hydro_hr_file = "hydro_hr.txt"

    location_file = "location.txt"

    # Parameter Check
    for i in range(0, num_of_landuse):
        if( ksv[i] > 0.0 and ka[i] > 0.0 ):
            raise Exception ("Error: both ksv and ka are non-zero.")
        if( gammam[i] > gammaa[i] ):
            raise Exception ("Error: gammam must be smaller than gammaa.")
    #enddo

    # Set da, dm and infilt_limit
    da = np.zeros(num_of_landuse)
    dm = np.zeros(num_of_landuse)
    infilt_limit = np.zeros(num_of_landuse)
    for i in range(0, num_of_landuse):
        if( soildepth[i] > 0.0 and ksv[i] > 0.0 ):
            infilt_limit[i] = soildepth[i] * gammaa[i]
        if( soildepth[i] > 0.0 and ka[i] > 0.0 ):
            da[i] = soildepth[i] * gammaa[i]
        if( soildepth[i] > 0.0 and ka[i] > 0.0 and gammam[i] > 0.0 ):
            dm[i] = soildepth[i] * gammam[i]
    #enddo

    # if ksg(i) = 0.d0 -> no gw calculation
    gw_switch = 0
    for i in range(0, num_of_landuse):
        if( ksg[i] > 0.0 ):
            gw_switch = 1
        else:
            gammag[i] = 0.0
            kg0[i] = 0.0
            fpg[i] = 0.0
            rgl[i] = 0.0
        #endif
    #enddo

    ##########################################################
    ### STEP 1: FILE READING 
    ##########################################################

    # max timestep
    maxt = lasth * 3600 / dt

    # dem file metadata (6 lines, kvp, space separated)
    f10 = open(datadir+demfile,"r")
    lines_list = f10.readlines()
    nx = int(lines_list[0].split(" ")[1])
    ny = int(lines_list[1].split(" ")[1])
    xllcorner = float(lines_list[2].split(" ")[1])
    yllcorner = float(lines_list[3].split(" ")[1])
    cellsize = float(lines_list[4].split(" ")[1])
    nodata = int([s for s in lines_list[5].split(' ') if s][1])
    print(nx, ny, xllcorner, yllcorner, cellsize, nodata)
    # Load DEM data
    zs = np.zeros((ny, nx))
    for l in range(6,len(lines_list)):
        line = lines_list[l]
        #line_list = line.split(" ")
        line_list = [s for s in line.split(' ') if s]
        for m in range(nx):
            zs[l-6][m] = float(line_list[m])

    f10.close()

    # Create remaining arrays
    zb = np.zeros((ny, nx))
    zb_riv = np.zeros((ny, nx))
    domain = np.zeros((ny, nx))

    # flow accumulation file
    riv = np.zeros((ny, nx))
    acc = np.zeros((ny, nx))
    f10 = open(datadir+accfile,"r")
    lines_list = f10.readlines()
    for l in range(6,len(lines_list)):
        line = lines_list[l]
        #line_list = line.split(" ")
        line_list = [s for s in line.split(' ') if s]
        for m in range(nx):
            acc[l-6][m] = line_list[m]

    f10.close()

    # flow direction file
    direc = np.zeros((ny, nx))
    f10 = open(datadir+dirfile,"r")
    lines_list = f10.readlines()
    for l in range(6,len(lines_list)):
        line = lines_list[l]
        #line_list = line.split(" ")
        line_list = [s for s in line.split(' ') if s]
        for m in range(nx):
            direc[l-6][m] = line_list[m]

    f10.close()
    #call read_gis_int(dirfile, dir)

    # landuse file
    land = np.ones((ny, nx))
    if( land_switch == 1 ):
        f10 = open(datadir+landfile,"r")
        lines_list = f10.readlines()
        for l in range(6,len(lines_list)):
            line = lines_list[l]
            #line_list = line.split(" ")
            line_list = [s for s in line.split(' ') if s]
            for m in range(nx):
                direc[l-6][m] = line_list[m]

        f10.close()
        #call read_gis_int(landfile, land)
    #endif

    # land : 1 ... num_of_landuse
    print( "num_of_landuse : %d", num_of_landuse )
    land = np.where( land <= 0, num_of_landuse, land)
    land = np.where( land > num_of_landuse, num_of_landuse, land)

    ##########################################################
    ### STEP 2: CALC PREPARATION 
    ##########################################################

    # dx, dy calc
    # d1: south side length
    x1 = xllcorner
    y1 = yllcorner
    x2 = xllcorner + nx * cellsize
    y2 = yllcorner
    if( utm == 0 ):
        d1 = hubeny_sub( x1, y1, x2, y2)

    # d2: north side length
    x1 = xllcorner
    y1 = yllcorner + ny * cellsize
    x2 = xllcorner + nx * cellsize
    y2 = yllcorner + ny * cellsize
    if( utm == 0 ):
        d2 = hubeny_sub( x1, y1, x2, y2 )

    # d3: west side length
    x1 = xllcorner
    y1 = yllcorner
    x2 = xllcorner
    y2 = yllcorner + ny * cellsize
    if( utm == 0 ):
        d3 = hubeny_sub( x1, y1, x2, y2 )

    # d4: east side length
    x1 = xllcorner + nx * cellsize
    y1 = yllcorner
    x2 = xllcorner + nx * cellsize
    y2 = yllcorner + ny * cellsize
    if( utm == 0 ):
        d4 = hubeny_sub( x1, y1, x2, y2 )

    if( utm == 1 ):
        dx = cellsize
        dy = cellsize
    else:
        dx = (d1 + d2) / 2.0 / float(nx)
        dy = (d3 + d4) / 2.0 / float(ny)
    #endif
    print( "dx [m] : %f" % (dx))
    print( "dy [m] : %f" % (dy))

    # length and area of each cell
    length = np.sqrt(dx * dy)
    area = dx * dy

    # river widhth, depth, leavy height, river length, river area ratio
    width = np.zeros((ny, nx))
    depth = np.zeros((ny, nx))
    height = np.zeros((ny, nx))
    len_riv = np.zeros((ny, nx))
    area_ratio = np.zeros((ny, nx))

    # slope cell
    if( riv_thresh > 0 ):
        acc = np.where(acc > riv_thresh, 1, acc) # river cell
    #endif

    width = np.where(riv == 1, width_param_c * ( acc * dx * dy * 1.e-6 ) ** width_param_s, width)
    depth = np.where(riv == 1, depth_param_c * ( acc * dx * dy * 1.e-6 ) ** depth_param_s, depth)
    height = np.where((riv == 1) & (acc > height_limit_param), height_param, height)

    # river data is replaced by the information in files
    if( rivfile_switch >= 1 ):
        riv_thresh = 1
        f10 = open(datadir+widthfile,"r")
        lines_list = f10.readlines()
        for l in range(6,len(lines_list)):
            line = lines_list[l]
            #line_list = line.split(" ")
            line_list = [s for s in line.split(' ') if s]
            for m in range(nx):
                width[l-6][m] = line_list[m]

        f10.close()
        #call read_gis_real(widthfile, width)
        riv = np.where(width > 0, 1, riv) # river cells (if width >= 0.)

        f10 = open(datadir+depthfile,"r")
        lines_list = f10.readlines()
        for l in range(6,len(lines_list)):
            line = lines_list[l]
            #line_list = line.split(" ")
            line_list = [s for s in line.split(' ') if s]
            for m in range(nx):
                depth[l-6][m] = line_list[m]

        f10.close()
        #call read_gis_real(depthfile, depth)
        f10 = open(datadir+heightfile,"r")
        lines_list = f10.readlines()
        for l in range(6,len(lines_list)):
            line = lines_list[l]
            #line_list = line.split(" ")
            line_list = [s for s in line.split(' ') if s]
            for m in range(nx):
                height[l-6][m] = line_list[m]

        f10.close()
        #call read_gis_real(heightfile, height)
        height = np.where( height < 0.0, 0.0, height)
    #endif 
    len_riv = np.where(riv == 1, length, len_riv)

    # river cross section is set by section files
    sec_map = np.zeros((ny, nx))
    if( sec_switch == 1 ):
        f10 = open(datadir+sec_map_file,"r")
        lines_list = f10.readlines()
        for l in range(6,len(lines_list)):
            line = lines_list[l]
            #line_list = line.split(" ")
            line_list = [s for s in line.split(' ') if s]
            for m in range(nx):
                sec_map[l-6][m] = line_list[m]

        f10.close()
        #call read_gis_int(sec_map_file, sec_map)
        sec_id_max = np.max( sec_map )
        # call set_section
        width, sec_width, depth, sec_depth, height, sec_height, sec_hr, sec_area, sec_peri, sec_b, sec_ns_river, riv = set_section(sec_id_max, sec_map, sec_div, secfile, depth, width, height, riv)
    #endif
    len_riv = np.where(riv == 1, length, len_riv) # added on Dec. 27, 2021

    # river length is set by input file
    sec_length = np.zeros((ny, nx))
    if( sec_length_switch == 1 ):
        f10 = open(datadir+sec_length_file,"r")
        lines_list = f10.readlines()
        for l in range(6,len(lines_list)):
            line = lines_list[l]
            #line_list = line.split(" ")
            line_list = [s for s in line.split(' ') if s]
            for m in range(nx):
                sec_length[l-6][m] = line_list[m]

        f10.close()
        #call read_gis_real(sec_length_file, sec_length)
        len_riv = np.where(sec_length > 0.0, sec_length, len_riv)
    #endif

    if( rivfile_switch == 2 ):
        # levee on both river and slope grid cells : zs is increased with height
        zs = np.where( height > 0.0, zs + height, zs)
    else:
        # levee on only slope grid cells : zs is increased with height
        zs = np.where( (height > 0.0) & (riv == 0), zs + height, zs)
    #endif

    #where(riv == 1) area_ratio = width / length
    #where(riv == 1) area_ratio = width / len_riv
    area_ratio = np.where(riv == 1, width * len_riv / area, area_ratio) # modified by T.Sayama on Nov. 27, 2021

    zb_riv = zs
    for i in range(ny):
        for j in range(nx):
            zb[i, j] = zs[i, j] - soildepth[int(land[i,j])-1]
            if(riv[i, j] == 1):
                zb_riv[i, j] = zs[i, j] - depth[i, j]
        #enddo
    #enddo

    # domain setting
    # domain = 0 : outside the domain
    # domain = 1 : inside the domain
    # domain = 2 : outlet point (where dir[i,j] = 0 or dir[i,j] = -1),
    #              and cells located at edges
    domain.fill( 0 )
    num_of_cell = 0 
    for i in range(ny):
        for j in range(nx):
            if( zs[i, j] > -100.0 ):
                domain[i, j] = 1
                if( direc[i, j] == 0 ):
                    domain[i, j] = 2
                if( direc[i, j] == -1 ):
                    domain[i, j] = 2
            num_of_cell = num_of_cell + 1
            #endif
        #enddo
    #enddo
    print( "num_of_cell : ", num_of_cell)
    print( "total area [km2] : ", num_of_cell * area / (10.0 ** 6.00))

    # river index setting
    riv_idx2i, riv_idx2j, riv_ij2idx, down_riv_idx, domain_riv_idx, width_idx, depth_idx, height_idx, area_ratio_idx, zb_riv_idx, dis_riv_idx, dif_riv_idx, sec_map_idx, len_riv_idx, riv_count = riv_idx_setting(ny, nx, dx, dy, domain, riv, width, depth, height, area_ratio, zb_riv, dif, land, sec_map, len_riv, direc, sec_length_switch)

    # slope index setting
    slo_idx2i, slo_idx2j, slo_ij2idx, down_slo_idx, domain_slo_idx, zb_slo_idx, dis_slo_idx, len_slo_idx, acc_slo_idx, down_slo_1d_idx, dis_slo_1d_idx, len_slo_1d_idx, land_idx, dif_slo_idx, ns_slo_idx, soildepth_idx, gammaa_idx, ksv_idx, faif_idx, infilt_limit_idx, ka_idx, gammam_idx, beta_idx, da_idx, dm_idx, ksg_idx, gammag_idx, kg0_idx, fpg_idx, rgl_idx, slo_count = slo_idx_setting(ny, nx, domain, zb, acc, land, dif, ns_slope, soildepth, gammaa, ksv, faif, infilt_limit, ka, gammaa, beta, da, dm, ksg, gammag, kg0, fpg, rgl, eight_dir, dy, dx, direc )

    # TODO reading dam file
    damflg, dam_qin, dam_num, dam_name, dam_kind, dam_ix, dam_iy, dam_vol, dam_vol_temp, dam_volmax, dam_state, dam_qout, dam_loc, dam_floodq, dam_maxfloodq, dam_rate = dam_read(riv_count, dam_switch, damfile)

    # initial condition
    hs = np.zeros((ny, nx))
    hr = np.zeros((ny, nx))
    hg = np.zeros((ny, nx))
    gampt_ff = np.zeros((ny, nx))
    gampt_f = np.zeros((ny, nx))
    qrs = np.zeros((ny, nx))

    hr.fill(-0.10)
    hs.fill(-0.10)
    hg.fill(-0.10)

    #where(riv == 1) hr = init_cond_riv
    #where(domain == 1) hs = init_cond_slo
    hr = np.where(riv == 1, 0.0, hr)
    hs = np.where(domain == 1, 0.0, hs)
    hs = np.where(domain == 2, 0.0, hs)

    # for slope cells
    # if init_slo_switch = 1 => read from file
    if(init_slo_switch == 1):
        inith = np.zeros((ny, nx))
        f13 = open(datadir+initfile_slo)
        lines_list = f13.readlines()
        for i in range( ny ):
            for j in range( nx ):
                inith[i,j] = lines_list[i].split(" ")[j]
        #enddo
        inith = np.where(inith <= 0.0, 0.0, inith)
        hs = np.where(domain == 1 and inith >= 0.0, inith, hs)
        del inith 
        f13.close()
    #endif

    # for river cells
    # if init_riv_switch = 1 => read from file
    if(init_riv_switch == 1):
        inith = np.zeros((ny, nx))
        f13 = open(datadir+initfile_riv)
        lines_list = f13.readlines()
        for i in range( ny ):
            for j in range ( nx ):
                inith[i,j] = lines_list[i].split(" ")[j]
        #enddo
        inith = np.where(inith <= 0.0, 0.0, inith)
        hr = np.where(riv == 1 and inith >= 0.0, inith, hr)
        del inith 
        f13.close()
    #endif

    # for slope cells
    # if init_gw_switch = 1 => read from file
    if(init_gw_switch == 1):
        inith = np.zeros((ny, nx))
        f13 = open(datadir+initfile_gw)
        lines_list = f13.readlines()
        for i in range( ny ):
            for j in range ( nx ):
                inith[i,j] = lines_list[i].split(" ")[j]
        #enddo
        inith = np.where(inith <= 0.0, 0.0, inith)
        hg = np.where(domain == 1 and inith >= 0.0, inith, hg)
        del inith
        f13.close()
    #endif

    # if init_gampt_ff_switch = 1 => read from file

    if(init_gampt_ff_switch == 1):
        inith = np.zeros((ny, nx))
        f13 = open(datadir+initfile_gampt_ff)
        lines_list = f13.readlines()
        for i in range( ny ):
            for j in range ( nx ):
                inith[i,j] = lines_list[i].split(" ")[j]
        #enddo
        inith = np.where(inith <= 0.0, 0.0, inith)
        gampt_ff = np.where(domain == 1, inith, gampt_ff)
        del inith 
        f13.close()
    #endif

    # TODO boundary conditions
    #call read_bound

    # div file
    if( div_switch == 1 ):
        f20 = open(datadir+divfile)
        lines_list = f20.readlines()
        div_id_max = len(lines_list)
        print( "div_id_max : ", div_id_max)
        div_org_idx = np.zeros(div_id_max)
        div_dest_idx = np.zeros(div_id_max)
        div_rate = np.zeros(div_id_max)
        for k in range( div_id_max ):
            div_org_i, div_org_j, div_dest_i, div_dest_j, div_rate[k] = lines_list[i].split(" ")
            div_org_idx[k] = riv_ij2idx[ div_org_i, div_org_j ]
            div_dest_idx[k] = riv_ij2idx[ div_dest_i, div_dest_j ]
        #enddo
        print( "done: reading div file" )
        f20.close()
    #endif

    # hydro file
    if( hydro_switch == 1 ):
        f5 = open(datadir+location_file)
        lines_list = f5.readlines()
        f1012 = open(datadir+hydro_file)
        f1013 = open(datadir+hydro_hr_file)
        maxhydro = len(lines_list) - 1
        hydro_i = np.zeros(maxhydro)
        hydro_j = np.zeros(maxhydro)
        for i in range( maxhydro ):
            ctemp = lines_list[i].split(" ")[0]
            hydro_i[i] = int(lines_list[i].split(" ")[1])
            hydro_j[i] = int(lines_list[i].split(" ")[2])
        #enddo
        f5.close()
    #endif

    # dynamic allocation
    qs_ave = np.zeros((4, ny, nx))
    qr_ave = np.zeros((ny, nx))
    qg_ave = np.zeros((4, ny, nx))

    qr_idx = np.zeros(riv_count)
    qr_ave_idx = np.zeros(riv_count)
    qr_ave_temp_idx = np.zeros(riv_count)
    hr_idx = np.zeros(riv_count)

    fr = np.zeros(riv_count)
    vr_temp = np.zeros(riv_count)
    hr_err = np.zeros(riv_count)
    vr_err = np.zeros(riv_count)
    vr_idx = np.zeros(riv_count)
    kr2 = np.zeros(riv_count)
    kr3 = np.zeros(riv_count)
    kr4 = np.zeros(riv_count)
    kr5 = np.zeros(riv_count)
    kr6 = np.zeros(riv_count)

    qs_idx = np.zeros((4, slo_count))
    qs_ave_idx = np.zeros((4, slo_count))
    qs_ave_temp_idx = np.zeros((4, slo_count))
    hs_idx = np.zeros(slo_count)

    qp_t_idx = np.zeros(slo_count)
    fs = np.zeros(slo_count)
    hs_temp = np.zeros(slo_count)
    hs_err = np.zeros(slo_count)
    ks2 = np.zeros(slo_count)
    ks3 = np.zeros(slo_count)
    ks4 = np.zeros(slo_count)
    ks5 = np.zeros(slo_count)
    ks6 = np.zeros(slo_count)

    qg_idx = np.zeros((4, slo_count))
    qg_ave_idx = np.zeros((4, slo_count))
    qg_ave_temp_idx = np.zeros((4, slo_count))
    hg_idx = np.zeros(slo_count)
    fg = np.zeros(slo_count)
    hg_temp = np.zeros(slo_count)
    hg_err = np.zeros(slo_count)
    kg2 = np.zeros(slo_count)
    kg3 = np.zeros(slo_count)
    kg4 = np.zeros(slo_count)
    kg5 = np.zeros(slo_count)
    kg6 = np.zeros(slo_count)
    gampt_ff_idx = np.zeros(slo_count)
    gampt_f_idx = np.zeros(slo_count)

    rain_i = np.zeros(ny)
    rain_j = np.zeros(nx)
    qe_t_idx = np.zeros(slo_count)
    evp_i = np.zeros(ny)
    evp_j = np.zeros(nx)
    aevp = np.zeros((ny, nx))
    aevp_tsas = np.zeros(slo_count)
    exfilt_hs_tsas = np.zeros(slo_count)
    rech_hs_tsas = np.zeros(slo_count)

    # gw initial setting
    if(init_gw_switch != 1):
        hg_idx = hg_init(slo_count)
        hg = sub_slo_idx2ij( hg_idx )
    #endif

    # initial storage calculation
    rain_sum = 0.0
    aevp_sum = 0.0
    pevp_sum = 0.0
    sout = 0.0
    si = 0.0
    sg = 0.0
    ss, sr, si, sg = storage_calc(hs, hr, hg, domain, area, riv_thresh, riv, gampt_ff, gammag_idx, slo_ij2idx)
    sinit = ss + sr + si + sg
    # Write to file 1000
    f1000 = open(outfile_storage)
    f1000.write( rain_sum, pevp_sum, aevp_sum, sout, ss + sr + si + sg, (rain_sum - aevp_sum - sout - (ss + sr + si + sg) + sinit), ss, sr, si, sg )
    f1000.close()

    # reading rainfall data
    f11 = open(rainfile)
    lines_list = f11.readlines()
    t_rain[0], nx_rain, ny_rain = lines_list[0].split(" ")
    tt_max_rain = (len(lines_list) / (nx_rain + 1)) - 1
    print("tt_max_rain = %f" %(tt_max_rain))
    t_rain = np.zeros(tt_max_rain)
    qp = np.zeros((tt_max_rain, ny_rain, nx_rain))
    qp_t = np.zeros((ny, nx))
    print(tt_max_rain, nx_rain, ny_rain)
    qp = 0.0
    qp_t = 0.0 # added by T.Sayamaa on Dec 7, 2022 v1.4.2.7
    for tt in range(0, tt_max_rain, nx_rain + 1):
        if ( tt % ( nx_rain + 1 ) == 0 ):
            t_rain[tt], nx_rain, ny_rain = lines_list[tt].split(" ")
            for i in range( ny_rain ):
                for j in range( nx_rain ):
                    qp[tt, i, j] = lines_list[tt+i].split(" ")[j]
                #enddo
            #enddo
        #ENDIF
    #enddo
    # unit convert from (mm/h) to (m/s)
    qp = qp / 3600.0 / 1000.0
    for j in range( nx ):
        rain_j[j] = int( (xllcorner + (float[j] - 0.50) * cellsize - xllcorner_rain) / cellsize_rain_x ) + 1
    #enddo
    for i in range( ny ):
        rain_i[i] = ny_rain - int( (yllcorner + (float[ny] - float[i] + 0.50) * cellsize - yllcorner_rain) / cellsize_rain_y )
    #enddo
    f11.close()
    print( "done: reading rain file" )


    # reading evp data
    if( evp_switch != 0 ):
        f11 = open(evpfile)
        lines_list = f11.readlines()
        t_evp[0], nx_evp, ny_evp = lines_list[0].split(" ")
        tt_max_evp = (len(lines_list) / ( nx_evp + 1 )) - 1
        t_evp = np.zeros(tt_max_evp)
        qe = np.zeros((tt_max_evp, ny_evp, nx_evp))
        qe_t = np.zeros((ny, nx))
        for tt in range(0, tt_max_evp, nx_evp + 1 ):
            if ( tt % ( nx_rain + 1 ) == 0 ):
                t_evp[tt], nx_evp, ny_evp = lines_list[tt].split(" ")
                for i in range( ny_evp ):
                    for j in range( nx_evp ):
                        qe[tt, i, j] = lines_list[tt+i].split(" ")[j]
                    #enddo
                #enddo
            #endif
        #enddo
        # unit convert from (mm/h) to (m/s)
        qe = qe / 3600.0 / 1000.0
        for j in range( nx ):
            evp_j[j] = int( (xllcorner + (float[j] - 0.50) * cellsize - xllcorner_evp) / cellsize_evp_x ) + 1
        #enddo
        for i in range( ny ):
            evp_i[i] = ny_evp - int( (yllcorner + (float[ny] - float[i] + 0.50) * cellsize - yllcorner_evp) / cellsize_evp_y )
        #enddo
        f11.close()
        print( "done: reading evp file")
    #endif

    # For TSAS Output (Initial Condition)
    hs_idx = sub_slo_ij2idx( hs )
    hr_idx = sub_riv_ij2idx( hr )
    #call RRI_TSAS(0, hs_idx, hr_idx, hg_idx, qs_ave_idx, qr_ave_idx, qg_ave_idx, qp_t_idx)

    ##########################################################
    ### STEP 3: CALCULATION 
    ##########################################################

    rain_sum = 0.0
    aevp_sum = 0.0
    sout = 0.0

    # output timestep
    out_dt = float(maxt) / float(outnum)
    out_dt = max(1.0, out_dt)
    out_next = round(out_dt)
    tt = 0

    for t in range( maxt ):
        if(t % 1 == 0):
            print( t, "/", maxt )
        #******* RIVER CALCULATION ******************************
        if( riv_thresh >= 0 ):
            # if (riv_thresh < 0) go to 2 # TODO Deal with GOTO 2
            # from time = (t - 1) * dt to t * dt
            time = (t - 1) * dt # (current time)
            # time step is initially set to be "dt_riv"
            ddt = dt_riv
            ddt_chk_riv = dt_riv

            qr_ave = 0.0
            qr_ave_idx = 0.0
            if( dam_switch == 1 ):
                dam_vol_temp.fill(0.0)
            # hr -> hr_idx
            # Memo: riv_ij2idx must be here. 
            # hr_idx cannot be replaced within the following for loop.
            hr_idx = sub_riv_ij2idx( hr )
            # from hr_idx -> vr_idx
            for k in range( riv_count ):
                vr_idx[k] = hr2vr(hr_idx[k], k)
            #enddo

            #do #-----------------*******************
            while( time < t * dt):
                # "time + ddt" should be less than "t * dt"
                if(time + ddt > t * dt ):
                    ddt = t * dt - time
                    # boundary condition for river (water depth boundary)
                    if( bound_riv_wlev_switch >= 1 ):
                        itemp = -1
                    for jtemp in range( tt_max_bound_riv_wlev):
                        if( t_bound_riv_wlev[jtemp-1] < (time + ddt) and (time + ddt) <= t_bound_riv_wlev[jtemp] ):
                            itemp = jtemp
                    #enddo
                    for k in range( riv_count ):
                        if( bound_riv_wlev_idx[itemp, k] <= -100.0 ):
                            continue # not boundary
                        hr_idx[k] = bound_riv_wlev_idx[itemp, k]
                        vr_idx[k] = hr2vr(hr_idx[k], k)
                    #enddo
                #endif

                #1 continue
                while( errmax > 1.0 and ddt > ddt_min_riv): # modified on Jan 7, 2021
                    # try smaller ddt
                    ddt = np.max( safety * ddt * (errmax ** pshrnk), 0.50 * ddt )
                    ddt = np.max( ddt, ddt_min_riv ) # added on Jan 7, 2021
                    ddt_chk_riv = ddt
                    print( "shrink (riv): %f, %f, %f" % (ddt, errmax, maxloc( vr_err )))
                    if(ddt == 0):
                        raise Exception ('stepsize underflow')
                    if(dam_switch == 1 ):
                        dam_vol_temp[:] = 0.0
                    #go to 1
                    qr_ave_temp_idx[:] = 0.0

                    # Adaptive Runge-Kutta 
                    # (1)
                    fr, qr_idx = funcr( riv_count, vr_idx )
                    vr_temp = vr_idx + b21 * ddt * fr
                    vr_temp = np.where(vr_temp < 0, 0.0, vr_temp)
                    qr_ave_temp_idx = qr_ave_temp_idx + qr_idx * ddt

                    # (2)
                    kr2, qr_idx = funcr( riv_count, vr_temp )
                    vr_temp = vr_idx + ddt * (b31 * fr + b32 * kr2)
                    vr_temp = np.where(vr_temp < 0, 0.0, vr_temp)
                    qr_ave_temp_idx = qr_ave_temp_idx + qr_idx * ddt

                    # (3)
                    kr3, qr_idx = funcr( riv_count, vr_temp )
                    vr_temp = vr_idx + ddt * (b41 * fr + b42 * kr2 + b43 * kr3)
                    vr_temp = np.where(vr_temp < 0, 0.0, vr_temp)
                    qr_ave_temp_idx = qr_ave_temp_idx + qr_idx * ddt

                    # (4)
                    kr4, qr_idx = funcr( riv_count, vr_temp )
                    vr_temp = vr_idx + ddt * (b51 * fr + b52 * kr2 + b53 * kr3 + b54 * kr4)
                    vr_temp = np.where(vr_temp < 0, 0.0, vr_temp)
                    qr_ave_temp_idx = qr_ave_temp_idx + qr_idx * ddt

                    # (5)
                    kr5, qr_idx = funcr( riv_count, vr_temp )
                    vr_temp = vr_idx + ddt * (b61 * fr + b62 * kr2 + b63 * kr3 + b64 * kr4 + b65 * kr5)
                    vr_temp = np.where(vr_temp < 0, 0.0, vr_temp)
                    qr_ave_temp_idx = qr_ave_temp_idx + qr_idx * ddt

                    # (6)
                    kr6, qr_idx = funcr( riv_count, vr_temp )
                    vr_temp = vr_idx + ddt * (c1 * fr + c3 * kr3 + c4 * kr4 + c6 * kr6)
                    vr_temp = np.where(vr_temp < 0, 0.0, vr_temp)
                    qr_ave_temp_idx = qr_ave_temp_idx + qr_idx * ddt

                    # (e)
                    vr_err = ddt * (dc1 * fr + dc3 * kr3 + dc4 * kr4 + dc5 * kr5 + dc6 * kr6)

                    hr_err[:] = vr_err[:] / (area * area_ratio_idx[:])

                    # error evaluation
                    hr_err = np.where( domain_riv_idx == 0, 0, hr_err)
                    errmax = maxval( hr_err ) / eps
                #else
                # modified on Jan 7, 2021
                if(ddt == ddt_min_riv):
                    kr6, qr_idx = funcr( riv_count, vr_temp )
                    qr_ave_temp_idx = qr_idx * ddt * 6.0
                #endif
                if(time + ddt > t * dt ):
                    ddt = t * dt - time
                time = time + ddt
                vr_idx = vr_temp
                qr_ave_idx = qr_ave_idx + qr_ave_temp_idx
                #endif
            #enddo --------------************
            qr_ave_idx = qr_ave_idx / float(dt) / 6.0

            for k in range( riv_count ):
                hr_idx[k] = vr2hr(vr_idx[k], k)
            #enddo

            # hr_idx -> hr, qr_ave_idx -> qr_ave
            hr = sub_riv_idx2ij( hr_idx )
            qr_ave = sub_riv_idx2ij( qr_ave_idx )

            if( dam_switch == 1 ):
                qr_ave = dam_checkstate(dam_num, dam_vol_temp, dam_vol, dam_state, dam_vol_max)

        #******* SLOPE CALCULATION ******************************
        #2 continue
        else:
            # from time = (t - 1) * dt to t * dt
            time = (t - 1) * dt  # (current time)
            # time step is initially set to be "dt"
            ddt = dt
            ddt_chk_slo = dt

            qs_ave = 0.0
            qs_ave_idx = 0.0

            # hs -> hs_idx
            # Memo: slo_ij2idx must be here. 
            # hs_idx cannot be replaced within the following for loop.
            hs_idx = sub_slo_ij2idx( hs )
            gampt_ff_idx = sub_slo_ij2idx( gampt_ff ) # modified by T.Sayama on June 10, 2017

            while (time < t * dt): #do
                if(time + ddt > t * dt ):
                    ddt = t * dt - time
                # rainfall
                itemp = -1
                for jtemp in range( tt_max_rain ):
                    if( t_rain(jtemp-1) < (time + ddt) and (time + ddt) <= t_rain(jtemp) ):
                        itemp = jtemp
                #enddo
                for i in range( ny ):
                    if(rain_i[i] < 1 or rain_i[i] > ny_rain ):
                        continue
                    for j in range( nx ):
                        if(rain_j[j] < 1 or rain_j[j] > nx_rain ):
                            continue
                        qp_t[i,j] = qp(itemp, rain_i[i], rain_j[j])
                    #enddo
                #enddo
                qp_t_idx = sub_slo_ij2idx( qp_t )
                # boundary condition for slope (water depth boundary)
                if( bound_slo_wlev_switch >= 1 ):
                    itemp = -1
                    for jtemp in rnage( tt_max_bound_slo_wlev ):
                        if( t_bound_slo_wlev(jtemp-1) < (time + ddt) and (time + ddt) <= t_bound_slo_wlev(jtemp) ):
                            itemp = jtemp
                    #enddo
                    for k in range( slo_count ):
                        if( bound_slo_wlev_idx(itemp, k) <= -100.0 ):
                            continue # not boundary
                        hs_idx[k] = bound_slo_wlev_idx(itemp, k)
                    #enddo
                #endif
                #3 continue
                #if(errmax > 1.0 and ddt >= ddt_min_slo):
                while (errmax > 1.0 and ddt > ddt_min_slo): # modified on Jan 7, 2021
                    # try smaller ddt
                    ddt = max( safety * ddt * (errmax ** pshrnk), 0.50 * ddt )
                    ddt = max( ddt, ddt_min_slo ) # added on Jan 7, 2021
                    ddt_chk_slo = ddt
                    print( "shrink (slo): ", ddt, errmax, maxloc( hs_err ))
                    if(ddt == 0):
                        raise Exception ('stepsize underflow')
                    #go to 3
                    qs_ave_temp_idx.fill(0.0)
                    # Adaptive Runge-Kutta 
                    # (1)
                    #call funcs( hs_idx, qp_t_idx, fs, qs_idx )
                    fs, qs_idx = funcs(hs_idx, qp_t_idx, fs, qs_idx, slo_count, zb_slo_idx, ns_slo_idx, ka_idx, da_idx, dm_idx, beta_idx, dif_slo_idx, lmax, down_slo_idx, down_slo_1d_idx, dis_slo_idx, dis_slo_1d_idx, len_slo_idx, len_slo_1d_idx, lev_p, lev_n, area, bound_slo_disc_switch, tt_max_bound_slo_disc, t_bound_slo_disc, time, ddt, bound_slo_disc_idx, direc, dif_slo_idx)
                    hs_temp = hs_idx + b21 * ddt * fs
                    hs_temp = np.where(hs_temp < 0, 0.0, hs_temp)
                    qs_ave_temp_idx = qs_ave_temp_idx + qs_idx * ddt
                    # (2)
                    #call funcs( hs_temp, qp_t_idx, ks2, qs_idx )
                    ks2, qs_idx = funcs(hs_temp, qp_t_idx, ks2, qs_idx, slo_count, zb_slo_idx, ns_slo_idx, ka_idx, da_idx, dm_idx, beta_idx, dif_slo_idx, lmax, down_slo_idx, down_slo_1d_idx, dis_slo_idx, dis_slo_1d_idx, len_slo_idx, len_slo_1d_idx, lev_p, lev_n, area, bound_slo_disc_switch, tt_max_bound_slo_disc, t_bound_slo_disc, time, ddt, bound_slo_disc_idx, direc, dif_slo_idx)
                    hs_temp = hs_idx + ddt * (b31 * fs + b32 * ks2)
                    hs_temp = np.where(hs_temp < 0, 0.0, hs_temp)
                    qs_ave_temp_idx = qs_ave_temp_idx + qs_idx * ddt
                    # (3)
                    #call funcs( hs_temp, qp_t_idx, ks3, qs_idx )
                    ks3, qs_idx = funcs(hs_temp, qp_t_idx, ks3, qs_idx, slo_count, zb_slo_idx, ns_slo_idx, ka_idx, da_idx, dm_idx, beta_idx, dif_slo_idx, lmax, down_slo_idx, down_slo_1d_idx, dis_slo_idx, dis_slo_1d_idx, len_slo_idx, len_slo_1d_idx, lev_p, lev_n, area, bound_slo_disc_switch, tt_max_bound_slo_disc, t_bound_slo_disc, time, ddt, bound_slo_disc_idx, direc, dif_slo_idx)
                    hs_temp = hs_idx + ddt * (b41 * fs + b42 * ks2 + b43 * ks3)
                    hs_temp = np.where(hs_temp < 0, 0.0, hs_temp)
                    qs_ave_temp_idx = qs_ave_temp_idx + qs_idx * ddt
                    # (4)
                    #call funcs( hs_temp, qp_t_idx, ks4, qs_idx )
                    ks4, qs_idx = funcs(hs_temp, qp_t_idx, ks4, qs_idx, slo_count, zb_slo_idx, ns_slo_idx, ka_idx, da_idx, dm_idx, beta_idx, dif_slo_idx, lmax, down_slo_idx, down_slo_1d_idx, dis_slo_idx, dis_slo_1d_idx, len_slo_idx, len_slo_1d_idx, lev_p, lev_n, area, bound_slo_disc_switch, tt_max_bound_slo_disc, t_bound_slo_disc, time, ddt, bound_slo_disc_idx, direc, dif_slo_idx)
                    hs_temp = hs_idx + ddt * (b51 * fs + b52 * ks2 + b53 * ks3 + b54 * ks4)
                    hs_temp = np.where(hs_temp < 0, 0.0, hs_temp)
                    qs_ave_temp_idx = qs_ave_temp_idx + qs_idx * ddt
                    # (5)
                    #call funcs( hs_temp, qp_t_idx, ks5, qs_idx )
                    ks5, qs_idx = funcs(hs_temp, qp_t_idx, ks5, qs_idx, slo_count, zb_slo_idx, ns_slo_idx, ka_idx, da_idx, dm_idx, beta_idx, dif_slo_idx, lmax, down_slo_idx, down_slo_1d_idx, dis_slo_idx, dis_slo_1d_idx, len_slo_idx, len_slo_1d_idx, lev_p, lev_n, area, bound_slo_disc_switch, tt_max_bound_slo_disc, t_bound_slo_disc, time, ddt, bound_slo_disc_idx, direc, dif_slo_idx)
                    hs_temp = hs_idx + ddt * (b61 * fs + b62 * ks2 + b63 * ks3 + b64 * ks4 + b65 * ks5)
                    hs_temp = np.where(hs_temp < 0, 0.0, hs_temp)
                    qs_ave_temp_idx = qs_ave_temp_idx + qs_idx * ddt
                    # (6)
                    #call funcs( hs_temp, qp_t_idx, ks6, qs_idx )
                    ks6, qs_idx = funcs(hs_temp, qp_t_idx, ks6, qs_idx, slo_count, zb_slo_idx, ns_slo_idx, ka_idx, da_idx, dm_idx, beta_idx, dif_slo_idx, lmax, down_slo_idx, down_slo_1d_idx, dis_slo_idx, dis_slo_1d_idx, len_slo_idx, len_slo_1d_idx, lev_p, lev_n, area, bound_slo_disc_switch, tt_max_bound_slo_disc, t_bound_slo_disc, time, ddt, bound_slo_disc_idx, direc, dif_slo_idx)
                    hs_temp = hs_idx + ddt * (c1 * fs + c3 * ks3 + c4 * ks4 + c6 * ks6)
                    hs_temp = np.where(hs_temp < 0, 0.0, hs_temp)
                    qs_ave_temp_idx = qs_ave_temp_idx + qs_idx * ddt
                    # (e)
                    hs_err = ddt * (dc1 * fs + dc3 * ks3 + dc4 * ks4 + dc5 * ks5 + dc6 * ks6)
                    # error evaluation
                    hs_err = np.where( domain_slo_idx == 0, 0.0, hs_err)
                    errmax = maxval( hs_err ) / eps
                #else:
                # "time + ddt" should be less than "t * dt"
                if(time + ddt > t * dt ):
                    ddt = t * dt - time
                time = time + ddt
                hs_idx = hs_temp
                qs_ave_idx = qs_ave_idx + qs_ave_temp_idx
                #endif

                # cumulative rainfall
                for i in range( ny ):
                    for j in range( nx ):
                        if( domain[i,j] != 0 ):
                            rain_sum = rain_sum + float(qp_t[i,j] * area * ddt)
                    #enddo
                #enddo
            #enddo
        qs_ave_idx = qs_ave_idx / float(dt) / 6.0 # modified on ver 1.4.1

        #******* GW CALCULATION ******************************
        #if( gw_switch == 0 ) go to 6
        if( gw_switch != 0):
            # from time = (t - 1) * dt to t * dt
            time = (t - 1) * dt  # (current time)
            # time step is initially set to be "dt"
            ddt = dt
            ddt_chk_slo = dt

            qg_ave = 0.0
            qg_ave_idx = 0.0

            # hg -> hg_idx
            # Memo: slo_ij2idx must be here. 
            # hg_idx cannot be replaced within the following for loop.
            hg_idx = sub_slo_ij2idx( hg )
            #call sub_slo_ij2idx( gampt_ff, gampt_ff_idx ) # modified by T.Sayama on June 10, 2017

            # GW Recharge
            hg_idx = gw_recharge( hs_idx, gampt_ff_idx )

            # GW Lose
            hg_idx = gw_lose( hg_idx )

            while( time < t * dt ):
                if( time + ddt > t * dt ):
                    ddt = t * dt - time
                #5 continue
                #if(errmax > 1.0 and ddt >= ddt_min_slo):
                while ( errmax > 1.0 and ddt > ddt_min_slo ): # modified on Jan 7, 2021
                    qg_ave_temp_idx[:,:] = 0.0

                    # Adaptive Runge-Kutta 
                    # (1)
                    qg_idx = funcg( hg_idx, fg )
                    hg_temp = hg_idx + b21 * ddt * fg
                    qg_ave_temp_idx = qg_ave_temp_idx + qg_idx * ddt

                    # (2)
                    qg_idx = funcg( hg_temp, kg2 )
                    hg_temp = hg_idx + ddt * (b31 * fg + b32 * kg2)
                    qg_ave_temp_idx = qg_ave_temp_idx + qg_idx * ddt

                    # (3)
                    qg_idx = funcg( hg_temp, kg3 )
                    hg_temp = hg_idx + ddt * (b41 * fg + b42 * kg2 + b43 * kg3)
                    qg_ave_temp_idx = qg_ave_temp_idx + qg_idx * ddt

                    # (4)
                    qg_idx = funcg( hg_temp, kg4 )
                    hg_temp = hg_idx + ddt * (b51 * fg + b52 * kg2 + b53 * kg3 + b54 * kg4)
                    qg_ave_temp_idx = qg_ave_temp_idx + qg_idx * ddt

                    # (5)
                    qg_idx = funcg( hg_temp, kg5 )
                    hg_temp = hg_idx + ddt * (b61 * fg + b62 * kg2 + b63 * kg3 + b64 * kg4 + b65 * kg5)
                    qg_ave_temp_idx = qg_ave_temp_idx + qg_idx * ddt

                    # (6)
                    qg_idx = funcg( hg_temp, kg6 )
                    hg_temp = hg_idx + ddt * (c1 * fg + c3 * kg3 + c4 * kg4 + c6 * kg6)
                    qg_ave_temp_idx = qg_ave_temp_idx + qg_idx * ddt

                    # (e)
                    hg_err = ddt * (dc1 * fg + dc3 * kg3 + dc4 * kg4 + dc5 * kg5 + dc6 * kg6)

                    # error evaluation
                    hg_err = np.where( domain_slo_idx == 0, 0.0, hg_err)
                    errmax = np.max( hg_err ) / eps

                    ##if(errmax > 1.0 and ddt >= ddt_min_slo):
                    #if(errmax > 1.0 and ddt > ddt_min_slo): # modified on Jan 7, 2021
                    # try smaller ddt
                    ddt = np.max( safety * ddt * (errmax ** pshrnk), 0.50 * ddt )
                    ddt = np.max( ddt, ddt_min_slo ) # added on Jan 7, 2021
                    ddt_chk_slo = ddt
                    print( "shrink (gw): %f, %f, %f " % (ddt, errmax, maxloc( hg_err )))
                    if(ddt == 0):
                        raise Exception ('stepsize underflow')
                    #go to 5
                #else:
                # "time + ddt" should be less than "t * dt"
                if( time + ddt > t * dt ):
                    ddt = t * dt - time
                time = time + ddt
                hg_idx = hg_temp
                qg_ave_idx = qg_ave_idx + qg_ave_temp_idx
                #endif

            #end while do

            qg_ave_idx = qg_ave_idx / float(dt) / 6.0

            time = t * dt

            #******* GW Exfiltration ********************************
            hg_idx = gw_exfilt( hs_idx, gampt_ff_idx )
        #ENDIF if( gw_switch == 0 ) go to 6
        #6 continue

        #******* Evapotranspiration *****************************
        if( evp_switch != 0 ):
            #call evp(hs_idx, gampt_ff_idx)
            aevp, aevp_tsas, hs_idx, gampt_ff_idx, aevp_sum, pevp_sum = evp( evp_switch, itemp, jtemp, t_evp, evp_i, evp_j, qe_t, slo_count, hs_idx, gampt_ff_idx, aevp, aevp_tsas, aevp_sum, tevp_sum )

        # hs_idx -> hs
        hs = sub_slo_idx2ij( hs_idx )
        qs_ave = sub_slo_idx2ij4( qs_ave_idx )
        hg = sub_slo_idx2ij( hg_idx )
        qg_ave = sub_slo_idx2ij4( qg_ave_idx )
        gampt_ff = sub_slo_idx2ij( gampt_ff_idx )

        #******* LEVEE BREAK ************************************
        #call levee_break(t, hr, hs, xllcorner, yllcorner, cellsize)

        #******* RIVER-SLOPE INTERACTIONS ***********************
        if( riv_thresh >= 0 ):
            #call funcrs(hr, hs)
            qrs, hr, hs = funcrs( hr, hs, ny, nx, domain, riv, depth, riv_ij2idx, len_riv_idx, height, dt, area, k, qrs )
        hr_idx = sub_riv_ij2idx( hr )
        hs_idx = sub_slo_ij2idx( hs )

        #******* INFILTRATION (Green Ampt) **********************
        #call infilt(hs_idx, gampt_ff_idx, gampt_f_idx)
        hs_idx, gampt_ff_idx, gampt_f_idx = infilt(hs_idx, gampt_ff_idx, gampt_f_idx, slo_count, ksv_idx, faif_idx, gamaa_idx)
        hs = sub_slo_idx2ij( hs_idx )
        gampt_ff = sub_slo_idx2ij( gampt_ff_idx )
        gampt_f = sub_slo_idx2ij( gampt_f_idx )

        #******* SET WATER DEPTH 0 AT DOMAIN = 2 ****************
        for i in range( ny ):
            for j in range( nx ):
                if( domain[i,j] == 2 ):
                    sout = sout + hs[i,j] * area
                    hs[i,j] = 0.0
                    if( riv[i,j] == 1):
                        vr_out = hr2vr(hr[i, j], riv_ij2idx[i,j])
                        sout = sout + vr_out
                        hr[i,j] = 0.0
                    #endif
                #endif
            #enddo
        #enddo

        # hs -> hs_idx, hr -> hr_idx, hg -> hg_idx
        hr_idx = sub_riv_ij2idx( hr )
        hs_idx = sub_slo_ij2idx( hs )
        hg_idx = sub_slo_ij2idx( hg )

        print( "max hr: ", maxval(hr), "loc : ", maxloc(hr))
        print( "max hs: ", maxval(hs), "loc : ", maxloc(hs))
        if(gw_switch == 1):
            print( "max hg: ", maxval(hg), "loc : ", maxloc(hg))

        #******* OUTPUT *****************************************
        # For TSAS Output
        #call RRI_TSAS(t, hs_idx, hr_idx, hg_idx, qs_ave_idx, qr_ave_idx, qg_ave_idx, qp_t_idx)
        if( hydro_switch == 1 and (int(time) % 3600) == 0 ):
            f1012.write(time)
            f1013.write(time)
            for k in range( maxhydro ):
                f1012.write(",")
                f1012.write(qr_ave[hydro_i[k], hydro_j[k]])
                f1013.write(",")
                f1013.write(hr[hydro_i[k], hydro_j[k]]) # added by T.Sayama on July 1, 2021
            f1012.write("\n")
            f1013.write("\n")

        # open output files
        if( t == out_next ):
            print( "OUTPUT :", t, time)
            tt = tt + 1
            out_next = round((tt+1) * out_dt)
            t_char = int2char(tt)
            hs = np.where(domain == 0, -0.10, hs)
            if(riv_thresh >= 0):
                hr = np.where(domain == 0, -0.10, hr)
            if(riv_thresh >= 0):
                qr_ave = np.where(domain == 0, -0.10, qr_ave)
            gampt_ff = np.where(domain == 0, -0.10, gampt_ff)
            aevp = np.where(domain == 0, -0.10, aevp)
            if( evp_switch != 0 ):
                qe_t = np.where(domain == 0, -0.10, qe_t)
            hg = np.where(domain == 0, -0.10, hg)
            if(outswitch_hs == 1):
                ofile_hs = trim(outfile_hs) + trim(t_char) + ".out"
            if(outswitch_hs == 2):
                ofile_hs = trim(outfile_hs) + trim(t_char) + ".bin"
            if(outswitch_hr == 1):
                ofile_hr = trim(outfile_hr) + trim(t_char) + ".out"
            if(outswitch_hr == 2):
                ofile_hr = trim(outfile_hr) + trim(t_char) + ".bin"
            if(outswitch_hg == 1):
                ofile_hg = trim(outfile_hg) + trim(t_char) + ".out"
            if(outswitch_hg == 2):
                ofile_hg = trim(outfile_hg) + trim(t_char) + ".bin"
            if(outswitch_qr == 1):
                ofile_qr = trim(outfile_qr) + trim(t_char) + ".out"
            if(outswitch_qr == 2):
                ofile_qr = trim(outfile_qr) + trim(t_char) + ".bin"
            if(outswitch_qu == 1):
                ofile_qu = trim(outfile_qu) + trim(t_char) + ".out"
            if(outswitch_qu == 2):
                ofile_qu = trim(outfile_qu) + trim(t_char) + ".bin"
            if(outswitch_qv == 1):
                ofile_qv = trim(outfile_qv) + trim(t_char) + ".out"
            if(outswitch_qv == 2):
                ofile_qv = trim(outfile_qv) + trim(t_char) + ".bin"
            if(outswitch_gu == 1):
                ofile_gu = trim(outfile_gu) + trim(t_char) + ".out"
            if(outswitch_gu == 2):
                ofile_gu = trim(outfile_gu) + trim(t_char) + ".bin"
            if(outswitch_gv == 1):
                ofile_gv = trim(outfile_gv) + trim(t_char) + ".out"
            if(outswitch_gv == 2):
                ofile_gv = trim(outfile_gv) + trim(t_char) + ".bin"
            if(outswitch_gampt_ff == 1):
                ofile_gampt_ff = trim(outfile_gampt_ff) + trim(t_char) + ".out"
            if(outswitch_gampt_ff == 2):
                ofile_gampt_ff = trim(outfile_gampt_ff) + trim(t_char) + ".bin"

            if(outswitch_hs == 1):
                f100 = open( ofile_hs )
            if(outswitch_hr == 1):
                f101 = open( ofile_hr )
            if(outswitch_hg == 1):
                f102 = open( ofile_hg )
            if(outswitch_qr == 1):
                f103 = open( ofile_qr )
            if(outswitch_qu == 1):
                f104 = open( ofile_qu )
            if(outswitch_qv == 1):
                f105 = open( ofile_qv )
            if(outswitch_gu == 1):
                f106 = open( ofile_gu )
            if(outswitch_gv == 1):
                f107 = open( ofile_gv )
            if(outswitch_gampt_ff == 1):
                f108 = open( ofile_gampt_ff )

            # This should open in binary format
            if(outswitch_hs == 2):
                f100 = open( ofile_hs, "wb" )
            if(outswitch_hr == 2):
                f101 = open( ofile_hr, "wb" )
            if(outswitch_hg == 2):
                f102 = open( ofile_hr, "wb" )
            if(outswitch_qr == 2):
                f103 = open( ofile_qr, "wb" )
            if(outswitch_qu == 2):
                f104 = open( ofile_qu, "wb" )
            if(outswitch_qv == 2):
                f105 = open( ofile_qv, "wb" )
            if(outswitch_gu == 2):
                f106 = open( ofile_gu, "wb" )
            if(outswitch_gv == 2):
                f107 = open( ofile_gv, "wb" )
            if(outswitch_gampt_ff == 2):
                f108 = open( ofile_gampt_ff, "wb" )

            # TODO output (ascii)
            if(outswitch_hs == 1):
                f100.write(hs)
            if(outswitch_hr == 1):
                f101.write(hr)
            #if(outswitch_hr == 1):
                #f101.write(hr + zb_riv)
            if(outswitch_hg == 1):
                f102.write(hg)
            if(outswitch_qr == 1):
                f103.write(qr_ave) # [m3/s]
            if(outswitch_qu == 1):
                f104.write((qs_ave[1,i,j] + (qs_ave[3,i,j] - qs_ave[4,i,j]) / 2.0) * area)
            #if(outswitch_qv == 1):
                #f105.write((qs_ave[2,i,j] + (qs_ave[3,i,j] + qs_ave[4,i,j]) / 2.0) * area)
            if(outswitch_qv == 1):
                f105.write((qs_ave[2,i,j] + (qs_ave[3,i,j] + qs_ave[4,i,j]) / 2.0) * area)
            if(outswitch_gu == 1):
                f106.write((qg_ave[1,i,j] + (qg_ave[3,i,j] - qg_ave[4,i,j]) / 2.0) * area)
            if(outswitch_gv == 1):
                f107.write((qg_ave[2,i,j] + (qg_ave[3,i,j] + qg_ave[4,i,j]) / 2.0) * area)
            if(outswitch_gampt_ff == 1):
                f108.write(gampt_ff)
            #enddo

            # TODO output (binary)
            if(outswitch_hs == 2):
                f100.write(hs)
            if(outswitch_hr == 2):
                f101.write(hr)
            if(outswitch_hg == 2):
                f102.write(hg)
            if(outswitch_qr == 2):
                f103.write(qr_ave) # [m3/s]
            if(outswitch_qu == 2):
                f104.write((qs_ave[1,i,j] + (qs_ave[3,i,j] - qs_ave[4,i,j]) / 2.0) * area)
            if(outswitch_qv == 2):
                f105.write((qs_ave[2,i,j] + (qs_ave[3,i,j] + qs_ave[4,i,j]) / 2.0) * area)
            if(outswitch_gu == 2):
                f106.write((qg_ave[1,i,j] + (qg_ave[3,i,j] - qg_ave[4,i,j]) / 2.0) * area)
            if(outswitch_gv == 2):
                f107.write((qg_ave[2,i,j] + (qg_ave[3,i,j] + qg_ave[4,i,j]) / 2.0) * area)
            if(outswitch_gampt_ff == 2):
                f108.write(gampt_ff)

            if(outswitch_hs != 0): 
                f100.close()
            if(outswitch_hr != 0): 
                f101.close()
            if(outswitch_hg != 0): 
                f102.close()
            if(outswitch_qr != 0): 
                f103.close()
            if(outswitch_qu != 0): 
                f104.close()
            if(outswitch_qv != 0): 
                f105.close()
            if(outswitch_gu != 0): 
                f106.close()
            if(outswitch_gv != 0): 
                f107.close()
            if(outswitch_gampt_ff != 0):
                f108.close()

            if( tec_switch == 1 ):
                if (tt == 1):
                    #call Tecout_alloc(nx, ny, 4)
                    iMX = ny + 1
                    jMX = nx + 1
                    ValNum = ValNum_temp
                    X = np.zeros(jMX)
                    Y = np.zeros(iMX)
                    Z = np.zeros((iMX,jMX))
                    Z_buf = np.zeros((ny,nx))
                    SufHmax = np.zeros((iMX-1,jMX-1))
                    SufHmax.fill(-0.10)
                    f6 = open( tecfile )
                    #call Tecout_mkGrid(dx, dy, zs)
                    Z = Tecout_mkGrid(dx, dy, zs, X, Y, Y_temp, Z_buf, Z_temp)
                    #call Tecout_write_initialize(tt, width, depth, height, area_ratio)
                    err = Tecout_write_initialize( f6, tt, iMX, jMX, X, Y, Z, width, depth, height, area_ratio)
                #endif
                #call Tecout_write(tt, qp_t, hr, qr_ave, hs, area)
                err = Tecout_write( tt, f6, iMX, jMX, qp_t, hr, qr_ave, hs, area, SufHmax )
            #endif

            # For dt_check
            #call dt_check_riv(hr_idx, tt, ddt_chk_riv)
            #call dt_check_slo(hs_idx, tt, ddt_chk_slo)

            if( dam_switch == 1 ):
                if( tt == 1 ):
                    f1001.open("./out/dam_out.txt")
                    #call dam_write
                    err = dam_write(f1001, dam_num, dam_vol, dam_qin, dam_loc, dam_qout)
                    f1001.close()
                if(t == maxt ):
                    f1002.open("./out/damcnt_out.txt")
                    #call dam_write_cnt
                    err = dam_write_cnt(f1002, dam_num, dam_name, dam_iy, dam_ix, dam_volmax, dam_floodq, dam_maxfloodq, dam_rate, dam_vol)
                    f1002.close()
                #endif
            #endif

            #******* OUTPUT FOR UNSTEADY MODEL (DIR = -1) *********** added on Sep 15, 2019
            itemp = 0
            for k in range( riv_count ):
                i = riv_idx2i[k]
                j = riv_idx2j[k]
                if(dir[i,j] == -1):
                    itemp = 1
            #enddo

            if( itemp == 1):
                ofile_ro = './out/ro_' + trim(t_char) + ".out"
                f111 = open( ofile_ro )
                for k in range( riv_count ):
                    i = riv_idx2i[k]
                    j = riv_idx2j[k]
                    kk = down_riv_idx[k]
                    if(k == kk):
                        continue
                    ii = riv_idx2i[kk]
                    jj = riv_idx2j[kk]
                    if(dir[ii,jj] == -1):
                        f111.write(1, i, j, ii, jj, qr_ave_idx[k])
                    #endif
                f111.close()
                #enddo
            for k in range( slo_count ):
                i = slo_idx2i[k]
                j = slo_idx2j[k]
                if(dir[i,j] == -1):
                    continue
                for l in range( lmax ):
                    kk = down_slo_idx[l, k]
                    if(kk <= 0):
                        continue
                    ii = slo_idx2i[kk]
                    jj = slo_idx2j[kk]
                    if(dir[ii,jj] == -1):
                        f111.write(2, i, j, ii, jj, qs_ave[l,i,j] * area)
                    #endif
                #enddo
            #enddo

            for k in range( slo_count ):
                i = slo_idx2i[k]
                j = slo_idx2j[k]
                if(dir[i,j] == -1):
                    for l in range( lmax ):
                        kk = down_slo_idx[l, k]
                        if(kk <= 0):
                            continue
                        ii = slo_idx2i[kk]
                        jj = slo_idx2j[kk]
                        if(domain[ii,jj] != 1):
                            continue
                        f111.write( 3, ii, jj, i, j, -1.0 * qs_ave[l,ii,jj] * area)
                    #enddo
                #endif
            #enddo

            # added on Feb 14, 2021 to include rainfall on grid cells with dir = -1
            for k in range(1, slo_count):
                i = slo_idx2i[k]
                j = slo_idx2j[k]
                if(dir[i,j] == -1):
                    f111.write( 4, i, j, i, j, qp_t_idx[k] * area)
                #endif
            #enddo
            f111.close()
        #endif
    #endif
    # check water balance
    if(t % 1 == 0):
        #call storage_calc(hs, hr, hg, ss, sr, si, sg)
        ss, sr, si, sg = storage_calc(hs, hr, hg, domain, area, riv_thresh, riv, gampt_ff, gammag_idx, slo_ij2idx)
        print( rain_sum, pevp_sum, aevp_sum, sout, ss + sr + si + sg, (rain_sum - aevp_sum - sout - (ss + sr + si + sg) + sinit))
        f1000.write( rain_sum, pevp_sum, aevp_sum, sout, ss + sr + si + sg, (rain_sum - aevp_sum - sout - (ss + sr + si + sg) + sinit), ss, sr, si, sg)
    #endif
    #enddo
    #pause
    #end program RRI
