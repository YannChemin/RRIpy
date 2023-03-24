# THIS IS RRI_Mod2.f90

#character*256 rainfile
#character*256 demfile
#character*256 accfile
#character*256 dirfile

int rivfile_switch
#character*256 widthfile
#character*256 depthfile
#character*256 heightfile

int init_slo_switch, init_riv_switch, init_gw_switch, init_gampt_ff_switch
#character*256 initfile_slo
#character*256 initfile_riv
#character*256 initfile_gw
#character*256 initfile_gampt_ff

int bound_slo_wlev_switch, bound_riv_wlev_switch
#character*256 boundfile_slo_wlev
#character*256 boundfile_riv_wlev

int bound_slo_disc_switch, bound_riv_disc_switch
#character*256 boundfile_slo_disc
#character*256 boundfile_riv_disc

int land_switch
#character*256 landfile

int div_switch
#character*256 divfile

int evp_switch
#character*256 evpfile

int sec_length_switch
#character*256 sec_length_file

int sec_switch
#character*256 sec_map_file
#character*256 sec_file

#int emb_switch
##character*256 embrfile, embbfile

int outswitch_hs
int outswitch_hr
int outswitch_hg
int outswitch_qr
int outswitch_qu
int outswitch_qv
int outswitch_gu
int outswitch_gv
int outswitch_gampt_ff
int outswitch_storage

#character*256 outfile_hs
#character*256 outfile_hr
#character*256 outfile_hg
#character*256 outfile_qr
#character*256 outfile_qu
#character*256 outfile_qv
#character*256 outfile_gu
#character*256 outfile_gv
#character*256 outfile_gampt_ff
#character*256 outfile_storage

int hydro_switch
#character*256 location_file

int lasth
int dt
int dt_riv
int outnum
float xllcorner_rain
float yllcorner_rain
float cellsize_rain_x, cellsize_rain_y

int utm
int eight_dir

float ns_river
int num_of_landuse
int dif[:]
float ns_slope[:]
float soildepth[:]
float gammaa[:]

float ksv[:]
float faif[:]
float infilt_limit[:]
# discharge per unit area

float ka[:]
float gammam[:]
float beta[:]
float da[:], dm[:]

float ksg[:]
float gammag[:]
float kg0[:]
float fpg[:]
float rgl[:]
int gw_switch

int riv_thresh
float width_param_c
float width_param_s
float depth_param_c
float depth_param_s
float height_param
int height_limit_param

int maxt, div
float ddt

int nx, ny, num_of_cell
float xllcorner, yllcorner
float cellsize
float length, area, dx, dy
int i4
parameter( i4 = 4 )

int domain[:,:], dir[:,:]
float zs[:,:], zb[:,:], zb_riv[:,:]
int riv[:,:], acc[:,:], land[:,:]
float width[:,:], depth[:,:], height[:,:], area_ratio[:,:], len_riv[:,:] # v1.4 add (len_riv[:,:]) 
float bound_slo_wlev[:,:], bound_riv_wlev[:,:]
float bound_slo_disc[:,:], bound_riv_disc[:,:]

float time
int tt_max_bound_slo_wlev, tt_max_bound_riv_wlev
int tt_max_bound_slo_disc, tt_max_bound_riv_disc
int :: t_bound_slo_wlev[:], t_bound_riv_wlev[:]
int :: t_bound_slo_disc[:], t_bound_riv_disc[:]
float gampt_ff[:,:], gampt_f[:,:]

int riv_count
int riv_idx2i[:], riv_idx2j[:], riv_ij2idx[:,:]
int down_riv_idx[:], domain_riv_idx[:]
float width_idx[:], depth_idx[:], height_idx[:], area_ratio_idx[:]
float zb_riv_idx[:], dis_riv_idx[:], len_riv_idx[:] # v1.4 add (len_riv_idx[:])
float bound_slo_wlev_idx[:,:], bound_riv_wlev_idx[:,:]
float bound_slo_disc_idx[:,:], bound_riv_disc_idx[:,:]
int dif_riv_idx[:]

int lmax, slo_count
int slo_idx2i[:], slo_idx2j[:], slo_ij2idx[:,:]
int down_slo_idx[:,:], domain_slo_idx[:], land_idx[:]
int down_slo_1d_idx[:]
float ns_slo_idx[:], soildepth_idx[:], gammaa_idx[:]
float ksv_idx[:], faif_idx[:], infilt_limit_idx[:]
float ka_idx[:], gammam_idx[:], beta_idx[:], da_idx[:], dm_idx[:]
float ksg_idx[:], gammag_idx[:], kg0_idx[:], fpg_idx[:], rgl_idx[:]
float zb_slo_idx[:], dis_slo_idx[:,:], len_slo_idx[:,:]
float dis_slo_1d_idx[:], len_slo_1d_idx[:]
int dif_slo_idx[:], acc_slo_idx[:]

#character*256 ofile_hs, ofile_hr, ofile_hg, ofile_qr, ofile_qu, ofile_qv, ofile_gu, ofile_gv, ofile_gampt_ff

int id_break

int div_id_max
int div_org_idx[:], div_dest_idx[:]
float div_rate[:]

float xllcorner_evp, yllcorner_evp
float cellsize_evp_x, cellsize_evp_y

int evp_i[:], evp_j[:]
int tt_max_evp
int t_evp[:]
int nx_evp, ny_evp
float qe[:,:,:], qe_t[:,:], qe_t_idx[:], aevp[:,:]
float aevp_sum, pevp_sum

float qrs[:,:]
int hs_id[:,:], hr_id[:,:]
int hs_id_idx[:], hr_id_idx[:]
float aevp_tsas[:], exfilt_hs_tsas[:], rech_hs_tsas[:]

int sec_id_max
int sec_map[:,:]
int sec_map_idx[:]
int sec_div[:]
float sec_length[:,:], sec_depth[:], sec_height[:]
float sec_hr[:,:], sec_area[:,:], sec_peri[:,:]
float sec_b[:,:], sec_ns_river[:,:], sec_length_idx[:]

#float emb_r[:,:], emb_b[:,:]
#float emb_r_idx[:], emb_b_idx[:]

#end module globals


# THIS IS RRI_Mod2.f90
# Variables used for "Adaptive Stepsize Control Runge-Kutta"

#module runge_mod

float eps, ddt_min_riv, ddt_min_slo
# for standard simulation
parameter( eps = 0.010 )
parameter( ddt_min_riv = 0.10 )
parameter( ddt_min_slo = 1.0 )

# for detailed simulation (new setting for detail)
#parameter( eps = 0.0050 )
#parameter( ddt_min_riv = 0.10 )
#parameter( ddt_min_slo = 1.0 )

# for rough simulation
#parameter( eps = 0.10 )
#parameter( ddt_min_riv = 1.0 )
#parameter( ddt_min_slo = 1.0 )

float safety, pgrow, pshrnk, errcon
parameter (safety=0.90,pgrow=-.20,pshrnk=-.250,errcon=1.89d-4)

float ks2[:], ks3[:], ks4[:], ks5[:], ks6[:]
float hs_temp[:], hs_err[:]

float kg2[:], kg3[:], kg4[:], kg5[:], kg6[:]
float hg_temp[:], hg_err[:]

float kr2[:], kr3[:], kr4[:], kr5[:], kr6[:]
#float hr_temp[:], hr_err[:]
float vr_temp[:], vr_err[:], hr_err[:] # v1.4

float a2,a3,a4,a5,a6,b21,b31,b32,b41,b42,b43,b51,b52,b53, &
        b54,b61,b62,b63,b64,b65,c1,c3,c4,c6,dc1,dc3,dc4,dc5,dc6
parameter (a2=.20,a3=.30,a4=.60,a5=1.0,a6=.8750,b21=.20,b31=3.0/40.0, &
     b32=9.0/40.0,b41=.30,b42=-.90,b43=1.20,b51=-11.0/54.0,b52=2.50, &
     b53=-70.0/27.0,b54=35.0/27.0,b61=1631.0/55296.0,b62=175.0/512.0, &
     b63=575.0/13824.0,b64=44275.0/110592.0,b65=253.0/4096.0,c1=37.0/378.0, &
     c3=250.0/621.0,c4=125.0/594.0,c6=512.0/1771.0,dc1=c1-2825.0/27648.0, &
     dc3=c3-18575.0/48384.0,dc4=c4-13525.0/55296.0,dc5=-277.0/14336.0, &
     dc6=c6-.250)

float errmax

#end module runge_mod


#THIS IS RRI_Mod_Dam.f90
#module dam_mod

int dam_switch
#character*256 damfile                            # dam control file
    
int dam_num                               # number of dam
#character*256 dam_name[:]        # dam name

int dam_kind[:]              # dam id number
int dam_ix[:], dam_iy[:]     # dam location (x, y)
int dam_loc[:]               # dam location (k)
int dam_state[:]             # dam status (full:1, not full:0)
int damflg[:]                # dam exist at each grid-cell (exist:1, not exist:0)

float dam_qin[:]               # dam inflow
float dam_vol[:]               # storage volume
float dam_vol_temp[:]          # storage volume (temporary)
float dam_volmax[:]            # maximum storage
float dam_qout[:]              # dam outflow
float dam_floodq[:]            # flood discharge
float dam_maxfloodq[:]         # max flood discharge (for non-const. cont.)
float dam_rate[:]              # discharge increasing rate

#end module dam_mod


#THIS IS RRI_Mod_Tecout.f90
#module tecout_mod
    int tec_switch
    #character*256 tecfile
    int iMX,jMX,ValNum
    float X[:],Y[:],Z[:,:],Z_buf[:,:]
    float SufHmax[:,:]
    #--------
    #contains
    #--------
    def alloc_Vals(nx,ny,ValNum_temp):
        #int :: nx,ny,ValNum_temp
        iMX = ny + 1
        jMX = nx + 1
        ValNum = ValNum_temp
        allocate( X[jMX],Y[iMX],Z[iMX,jMX],Z_buf[ny,nx] )
        allocate( SufHmax(iMX-1,jMX-1) )
        SufHmax[:,:] = -0.10
    #end subroutine alloc_Vals
#end module tecout_mod
