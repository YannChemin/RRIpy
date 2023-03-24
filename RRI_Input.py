version = "RRI_Input_Format_Ver1_4_2"
datadir = "solo30s"

rainfile = "./rain/rain_.dat"
demfile = "./topo/adem.txt"
accfile = "./topo/acc.txt"
dirfile = "./topo/adir.txt"

utm = 0				# utm(1) or latlon(0)
eight_dir = 1			# 4-direction (0), 8-direction(1)
lasth = 360			# lasth(hour)
dt = 600			# dt(second)
dt_riv = 60			# dt_riv
outnum = 60			# outnum [-]"
xllcorner_rain = 110.400000	# xllcorner_rain
yllcorner_rain = -8.158333,	# yllcorner_rain
cellsize_rain_x = 0.0083333333	# cellsize_rain"
cellsize_rain_y = 0.0083333333	# cellsize_rain"

ns_river = 3.000d-2		# ns_river
num_of_landuse 	= 1		# num_of_landuse
dif 		= [1]		# diffusion(1) or kinematic(0), one per LU
ns_slope 	= [4.000d-1]    # ns_slope, one per LU
soildepth 	= [1.0000] 	# soildepth, one per LU
gammaa 		= [4.750d-1]    # gammaa, one per LU

ksv 		= [0.0000]     	# kv (m/s), one per LU
faif 		= [3.163d-1]   	# Sf (m), one per LU

ka 		= [0.0000]     	# ka, one per LU
gammam 		= [0.0000]	# gammam, one per LU
beta 		= [8.0000]     	# beta, one per LU

ksg		= [0.0]         # ksg (m/s), 1/LU -- set zero for no bedrock gw
gammag		= [0.0370]      # gammag (-), one per LU
kg0		= [5.7d-5]      # kg0 (m/s), one per LU
fpg		= [0.10]        # fg (-), one per LU
rgl		= [0.0]         # rgl, one per LU

riv_thresh = 20			# riv_thresh
width_param_c = 5.000		# width_param_c (2.5)
width_param_s = 3.50d-1		# width_param_s (0.4)
depth_param_c = 9.50d-1		# depth_param_c (0.1)
depth_param_s = 2.00d-1		# depth_param_s (0.4)
height_param = 0.000		# height_param
height_limit_param = 20		# height_limit_param

rivfile_switch = 1
widthfile = "./riv/width.txt"
depthfile = "./riv/depth.txt"
heightfile = "./riv/height.txt"

init_slo_switch = 0
init_riv_switch = 0
init_gw_switch = 0
init_gampt_ff_switch = 0
initfile_slo = "./init/hs_init.out"
initfile_riv = "./init/hr_init.out"
initfile_gw = "./init/hg_init.out"
initfile_gampt_ff = "./init/gamptff_init.out"

bound_slo_wlev_switch = 0
bound_riv_wlev_switch = 0
boundfile_slo_wlev = "./bound/hs_wlev_bound.txt"
boundfile_riv_wlev = "./bound/hr_wlev_bound.txt"

bound_slo_disc_switch = 0
bound_riv_disc_switch = 0
boundfile_slo_disc = "./bound/qs_bound.txt"
boundfile_riv_disc = "./bound/qr_bound.txt"

land_switch = 0
landfile = "./topo/landuse.txt"

dam_switch = 0
damfile = "./damcnt.txt"

div_switch = 0
div_switch = "./div.txt"

evp_switch = 0
evpfile = "./rain/Evp.dat"
xllcorner_evp = 100.000000		# xllcorner_evp
yllcorner_evp = 10.000000		# yllcorner_evp
cellsize_evp_x = 0.0083333300		# cellsize x
cellsize_evp_y = 0.0083333300		# cellsize y

sec_length_switch = 0
sec_length_file = "./riv/length.txt"

sec_switch = 0
sec_map_file = "./riv/sec_map.txt"
sec_file = "./riv/section/sec_"

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
outfile_hs = "./out/hs_"
outfile_hr = "./out/hr_"
outfile_hg = "./out/hg_"
outfile_qr = "./out/qr_"
outfile_qu = "./out/qu_"
outfile_qv = "./out/qv_"
outfile_gu = "./out/gu_"
outfile_gv = "./out/gv_"
outfile_gampt_ff = "./out/gampt_ff_"
outfile_storage = "./out/storage.dat"

hydro_switch = 1
location_file = "./location.txt"
