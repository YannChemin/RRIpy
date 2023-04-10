# RRI_Sub
import numpy as np

# river index setting
def riv_idx_setting(ny, nx, dx, dy, domain, riv, width, depth, height, area_ratio, zb_riv, dif, land, sec_map, len_riv, direc, sec_length_switch):
    """
    River index setting

    :param ny:
    :param nx:
    :param dx:
    :param dy:
    :param domain:
    :param riv:
    :param width:
    :param depth:
    :param height:
    :param area_ratio:
    :param zb_riv:
    :param dif:
    :param land:
    :param sec_map:
    :param len_riv:
    :param direc:
    :param sec_length_switch:

    :return: riv_idx2i
    :return: riv_idx2j
    :return: riv_ij2idx
    :return: down_riv_idx
    :return: domain_riv_idx
    :return: width_idx
    :return: depth_idx
    :return: height_idx
    :return: area_ratio_idx
    :return: zb_riv_idx
    :return: dis_riv_idx
    :return: dif_riv_idx
    :return: sec_map_idx
    :return: len_riv_idx
    :return: riv_count
    """
    riv_count = 0
    for i in range( ny ):
        for j in range( nx ):
            if( domain[i,j] != 0 and riv[i,j] == 1 ):
                riv_count = riv_count + 1
        #enddo
    #enddo

    riv_idx2i = np.zeros(riv_count)
    riv_idx2j = np.zeros(riv_count)
    riv_ij2idx = np.zeros((ny, nx))

    down_riv_idx = np.zeros(riv_count)
    domain_riv_idx = np.zeros(riv_count)

    width_idx = np.zeros(riv_count)
    depth_idx = np.zeros(riv_count)
    height_idx = np.zeros(riv_count)
    area_ratio_idx = np.zeros(riv_count)

    zb_riv_idx = np.zeros(riv_count)
    dis_riv_idx = np.zeros(riv_count)
    
    dif_riv_idx = np.zeros(riv_count)
    sec_map_idx = np.zeros(riv_count) # add v1.4
    len_riv_idx = np.zeros(riv_count) # add v1.4

    riv_count = 0
    for i in range( ny ):
        for j in range( nx ):
            if(domain[i,j] == 0 or riv[i,j] != 1):
                continue
            # domain[i, j] = 1 or 2 and riv[i, j] = 1
            riv_idx2i[riv_count] = i
            riv_idx2j[riv_count] = j
            domain_riv_idx[riv_count] = domain[i, j]
            width_idx[riv_count] = width[i, j]
            depth_idx[riv_count] = depth[i, j]
            height_idx[riv_count] = height[i, j]
            area_ratio_idx[riv_count] = area_ratio[i, j]
            zb_riv_idx[riv_count] = zb_riv[i, j]
            riv_ij2idx[i, j] = riv_count
            dif_riv_idx[riv_count] = dif[int(land[i, j])-1]
            sec_map_idx[riv_count] = sec_map[i, j]
            len_riv_idx[riv_count] = len_riv[i, j] # add v1.4
            riv_count = riv_count + 1
        #enddo
    #enddo

    # search for downstream gridcell (down_idx)
    riv_count = 0
    for i in range( ny ):
        for j in range( nx ):
            if(domain[i,j] == 0 or riv[i,j] != 1):
                continue
            # domain(i, j) = 1 or 2 and riv(i, j) = 1
            # right
            if( direc[i,j] == 1 ):
                ii = i
                jj = j + 1
                distance = dx
            # right down
            elif( direc[i,j] == 2 ):
                ii = i + 1
                jj = j + 1
                distance = np.sqrt( dx*dx + dy*dy )
            # down
            elif( direc[i,j] == 4 ):
                ii = i + 1
                jj = j
                distance = dy
            # left down
            elif( direc[i,j] == 8 ):
                ii = i + 1
                jj = j - 1
                distance = np.sqrt( dx*dx + dy*dy )
            # left
            elif( direc[i,j] == 16 ):
                ii = i
                jj = j - 1
                distance = dx
            # left up
            elif( direc[i,j] == 32 ):
                ii = i - 1
                jj = j - 1
                distance = np.sqrt( dx*dx + dy*dy )
            # up
            elif( direc[i,j] == 64 ):
                ii = i - 1
                jj = j
                distance = dy
            # right up
            elif( direc[i,j] == 128 ):
                ii = i - 1
                jj = j + 1
                distance = np.sqrt( dx*dx + dy*dy )
            elif( direc[i,j] == 0 or direc[i,j] == -1 ):
                ii = i
                jj = j
            else:
                raise Exception ( "dir[%d, %d] error %f" % ( i, j, direc[i, j]) )
            #endif

            # If the downstream cell is outside the domain, set domain[i, j] = 2
            if( ii < 0 or ii > ny-1 or jj < 0 or jj > nx-1 ):
                domain[i, j] = 2
                direc[i, j] = 0
                ii = i
                jj = j
            #endif
            if( domain[ii, jj] == 0 ):
                domain[i, j] = 2
                direc[i, j] = 0
                ii = i
                jj = j
            #endif

            if( riv[ii,jj] == 0 ):
                raise Exception ("riv[%d,%d] should be 1 [%d %d]" % (i, j, ii, jj) )
            #endif
            dis_riv_idx[riv_count] = distance
            down_riv_idx[riv_count] = riv_ij2idx[ii, jj]

            riv_count = riv_count + 1
        #enddo
    #enddo

    # add v1.4
    if( sec_length_switch == 1 ):
        for k in range(riv_count): 
            kk = down_riv_idx[k]
            dis_riv_idx[k] = ( len_riv_idx[k] + len_riv_idx[kk] ) / 2.0
        #enddo
    #endif

    return(riv_idx2i, riv_idx2j, riv_ij2idx, down_riv_idx, domain_riv_idx, width_idx, depth_idx, height_idx, area_ratio_idx, zb_riv_idx, dis_riv_idx, dif_riv_idx, sec_map_idx, len_riv_idx, riv_count)

#end def riv_idx_setting


def sub_riv_ij2idx( riv_count, a, riv_idx2i, riv_idx2j, a_idx ):
    """
    # 2D -> 1D (ij2idx)
    
    :param riv_count:
    :param a:
    :param riv_idx2i:
    :param riv_idx2j:

    :return: a_idx updated
    """
    #real(8) a(ny, nx), a_idx(riv_count)
    for k in range( riv_count ):
        a_idx[k] = a[riv_idx2i[k], riv_idx2j[k]]
    #enddo

    return(a_idx)

#end def sub_riv_ij2idx


def sub_riv_idx2ij( a_idx, a ):
    """
    # 1D -> 2D (idx2ij)

    :param a:
    :param riv_count:
    :param riv_idx2i:
    :param riv_idx2j:
    :param a_idx:

    :return: Updated a
    """
    #real(8) a_idx(riv_count), a(ny, nx)
    a.fill(0.0)
    for k in range( riv_count ):
        a[riv_idx2i[k], riv_idx2j[k]] = a_idx[k]
    #enddo
    
    return(a)

#end def sub_riv_idx2ij


def slo_idx_setting(ny, nx, domain, zb, acc, land, dif, ns_slope, soildepth, gammaa, ksv, faif, infilt_limit, ka, gammam, beta, da, dm, ksg, gammag, kg0, fpg, rgl, eight_dir, dy, dx, direc ):
    """
    Slope index setting

    :param ny:
    :param nx:
    :param domain:
    :param zb:
    :param acc:
    :param land:
    :param dif:
    :param ns_slope:
    :param soildepth:
    :param gammaa:
    :param ksv:
    :param faif:
    :param infilt_limit:
    :param ka:
    :param gammam:
    :param beta:
    :param da:
    :param dm:
    :param ksg:
    :param gammag:
    :param kg0:
    :param fpg:
    :param rgl:
    :param eight_dir: switch 0=4directions 1=8directions
    :param dy:
    :param dx:
    :param direc:

    :return: slo_idx2i
    :return: slo_idx2j
    :return: slo_ij2idx
    :return: down_slo_idx
    :return: domain_slo_idx
    :return: zb_slo_idx
    :return: dis_slo_idx
    :return: len_slo_idx
    :return: acc_slo_idx
    :return: down_slo_1d_idx
    :return: dis_slo_1d_idx
    :return: len_slo_1d_idx
    :return: land_idx
    :return: dif_slo_idx
    :return: ns_slo_idx
    :return: soildepth_idx
    :return: gammaa_idx
    :return: ksv_idx
    :return: faif_idx
    :return: infilt_limit_idx
    :return: ka_idx
    :return: gammam_idx
    :return: beta_idx
    :return: da_idx
    :return: dm_idx
    :return: ksg_idx
    :return: gammag_idx
    :return: kg0_idx
    :return: fpg_idx
    :return: rgl_idx
    :return: slo_count
    """
    #real(8) distance, len, l1, l2, l3
    #real(8) l1_kin, l2_kin, l3_kin
    slo_count = 0
    for i in range( ny ):
        for j in range( nx ):
            if(domain[i,j] != 0):
                slo_count = slo_count + 1
        #enddo
    #enddo

    slo_idx2i = np.zeros(slo_count)
    slo_idx2j = np.zeros(slo_count)
    slo_ij2idx = np.zeros((ny, nx))
    down_slo_idx = np.zeros((4, slo_count))
    domain_slo_idx = np.zeros(slo_count)
    zb_slo_idx = np.zeros(slo_count)
    dis_slo_idx = np.zeros((4, slo_count))
    len_slo_idx = np.zeros((4, slo_count))
    acc_slo_idx = np.zeros(slo_count)
    down_slo_1d_idx = np.zeros(slo_count)
    dis_slo_1d_idx = np.zeros(slo_count)
    len_slo_1d_idx = np.zeros(slo_count)
    land_idx = np.zeros(slo_count)

    dif_slo_idx = np.zeros(slo_count)
    ns_slo_idx = np.zeros(slo_count)
    soildepth_idx = np.zeros(slo_count)
    gammaa_idx = np.zeros(slo_count)

    ksv_idx = np.zeros(slo_count)
    faif_idx = np.zeros(slo_count)
    infilt_limit_idx = np.zeros(slo_count)
    ka_idx = np.zeros(slo_count)
    gammam_idx = np.zeros(slo_count)
    beta_idx = np.zeros(slo_count)
    da_idx = np.zeros(slo_count)
    dm_idx = np.zeros(slo_count)
    ksg_idx = np.zeros(slo_count)
    gammag_idx = np.zeros(slo_count)
    kg0_idx = np.zeros(slo_count)
    fpg_idx = np.zeros(slo_count)
    rgl_idx = np.zeros(slo_count)

    slo_count = 0
    for i in range( ny ):
        for j in range( nx ):
            if(domain[i,j] == 0):
                continue
            slo_idx2i[slo_count] = i
            slo_idx2j[slo_count] = j
            domain_slo_idx[slo_count] = domain[i, j]
            zb_slo_idx[slo_count] = zb[i, j]
            acc_slo_idx[slo_count] = acc[i, j]
            slo_ij2idx[i, j] = slo_count
            land_idx[slo_count] = land[i, j]

            dif_slo_idx[slo_count] = dif[int(land[i, j])-1]
            ns_slo_idx[slo_count] = ns_slope[int(land[i,j])-1]
            soildepth_idx[slo_count] = soildepth[int(land[i, j])-1]
            gammaa_idx[slo_count] = gammaa[int(land[i, j])-1]

            ksv_idx[slo_count] = ksv[int(land[i,j])-1]
            faif_idx[slo_count] = faif[int(land[i,j])-1]
            infilt_limit_idx[slo_count] = infilt_limit[int(land[i,j])-1]
  
            ka_idx[slo_count] = ka[int(land[i,j])-1]
            gammam_idx[slo_count] = gammam[int(land[i,j])-1]
            beta_idx[slo_count] = beta[int(land[i,j])-1]
            da_idx[slo_count] = da[int(land[i,j])-1]
            dm_idx[slo_count] = dm[int(land[i,j])-1]

            ksg_idx[slo_count] = ksg[int(land[i, j])-1]
            gammag_idx[slo_count] = gammag[int(land[i, j])-1]
            kg0_idx[slo_count] = kg0[int(land[i, j])-1]
            fpg_idx[slo_count] = fpg[int(land[i, j])-1]
            rgl_idx[slo_count] = rgl[int(land[i, j])-1]

            slo_count = slo_count + 1
        #enddo
    #enddo

    if( eight_dir == 1 ):
        # Hromadka etal, JAIH2006 (USE THIS AS A DEFAULT)
        # 8-direction
        lmax = 4
        l1 = dy / 2.0
        l2 = dx / 2.0
        l3 = np.sqrt(dx ** 2.0 + dy ** 2.0) / 4.0
    elif( eight_dir == 0 ):
        # 4-direction
        lmax = 2
        l1 = dy
        l2 = dx
        l3 = 0.0
    else:
        raise Exception ( "error: eight_dir should be 0 or 1.")
    #endif

    # search for downstream gridcell (down_slo_idx)
    slo_count = 0
    down_slo_idx.fill(-1)

    for i in range( ny ):
        for j in range( nx ):
            if(domain[i,j] == 0):
                continue
            # domain(i, j) = 1 or 2
            slo_count = slo_count + 1
            # 8-direction: lmax = 4, 4-direction: lmax = 2
            for l in range( lmax ): # (1: rightC2: down, 3: right down, 4: left down)
                if( l == 1 ):
                    ii = i
                    jj = j + 1
                    distance = dx
                    length = l1
                elif( l == 2 ):
                    ii = i + 1
                    jj = j
                    distance = dy
                    length = l2
                elif( l == 3 ):
                    ii = i + 1
                    jj = j + 1
                    distance = np.sqrt( dx * dx + dy * dy )
                    length = l3
                else:
                    ii = i + 1
                    jj = j - 1
                    distance = np.sqrt( dx * dx + dy * dy )
                    length = l3
                #endif

                if( ii > ny ):
                    continue
                if( jj > nx ):
                    continue
                if( ii < 0 ):
                    continue
                if( jj < 0 ):
                    continue
                if( domain[ii,jj] == 0 ):
                    continue

                down_slo_idx[l, slo_count] = slo_ij2idx[ii, jj]
                dis_slo_idx[l, slo_count] = distance
                len_slo_idx[l, slo_count] = length
            #enddo
        #enddo
    #enddo

    # search for downstream gridcell (down_slo_1d_idx) (used only for kinematic with 1-direction)
    slo_count = 0

    down_slo_1d_idx.fill(-1)
    dis_slo_1d_idx.fill(l1)
    len_slo_1d_idx.fill(dx)

    l1_kin = dy
    l2_kin = dx
    l3_kin = dx * dy / np.sqrt( dx ** 2.0 + dy ** 2.0 )

    for i in range( ny ):
        for j in range( nx ):
            if(domain[i,j] == 0):
                continue
            # domain(i, j) = 1 or 2
            # right
            if( direc[i,j] == 1 ):
                ii = i
                jj = j + 1
                distance = dx
                length = l1_kin
            # right down
            elif( direc[i,j] == 2 ):
                ii = i + 1
                jj = j + 1
                distance = np.sqrt( dx*dx + dy*dy )
                length = l3_kin
            # down
            elif( direc[i,j] == 4 ):
                ii = i + 1
                jj = j
                distance = dy
                length = l2_kin
            # left down
            elif( direc[i,j] == 8 ):
                ii = i + 1
                jj = j - 1
                distance = np.sqrt( dx*dx + dy*dy )
                length = l3_kin
            # left    
            elif( direc[i,j] == 16 ):
                ii = i
                jj = j - 1
                distance = dx
                length = l1_kin
            # left up
            elif( direc[i,j] == 32 ):
                ii = i - 1
                jj = j - 1
                distance = np.sqrt( dx*dx + dy*dy )
                length = l3_kin
            # up
            elif( direc[i,j] == 64 ):
                ii = i - 1
                jj = j
                distance = dy
                length = l2_kin
            # right up
            elif( direc[i,j] == 128 ):
                ii = i - 1
                jj = j + 1
                distance = np.sqrt( dx*dx + dy*dy )
                length = l3_kin
            elif( direc[i,j] == 0 or direc[i,j] == -1 ):
                ii = i
                jj = j
                distance = dx
                length = l1_kin
            else:
                raise Exception ("dir(%d, %d) error %d %d" % ( i, j, direc[i, j]))
            #endif

            if( ii > ny ):
                continue
            if( jj > nx ):
                continue
            if( ii < 0 ):
                continue
            if( jj < 0 ):
                continue
            if( domain[ii,jj] == 0 ): 
                continue

            down_slo_1d_idx[slo_count] = slo_ij2idx[ii, jj]
            dis_slo_1d_idx[slo_count] = distance
            len_slo_1d_idx[slo_count] = length

            slo_count = slo_count + 1
        #enddo
    #enddo

    return(slo_idx2i, slo_idx2j, slo_ij2idx, down_slo_idx, domain_slo_idx, zb_slo_idx, dis_slo_idx, len_slo_idx, acc_slo_idx, down_slo_1d_idx, dis_slo_1d_idx, len_slo_1d_idx, land_idx, dif_slo_idx, ns_slo_idx, soildepth_idx, gammaa_idx, ksv_idx, faif_idx, infilt_limit_idx, ka_idx, gammam_idx, beta_idx, da_idx, dm_idx, ksg_idx, gammag_idx, kg0_idx, fpg_idx, rgl_idx, slo_count)

#end def slo_idx_setting


def sub_slo_ij2idx( k, slo_count, a_idx, a, slo_idx2i, slo_idx2j ):
    """
    2D -> 1D (ij2idx)
    
    :param k:
    :param slo_count:
    :param a_idx:
    :param a:
    :param slo_idx2i:
    :param slo_idx2j:

    :return: Updated a_idx
    """
    #real(8) a(ny, nx), a_idx(slo_count)
    for k in range ( slo_count ):
        a_idx[k] = a[slo_idx2i[k], slo_idx2j[k]]
    #enddo

    return(a_idx)

#end def sub_slo_ij2idx


def sub_slo_idx2ij( a, slo_count, slo_idx2i, slo_idx2j, a_idx ):
    """
    1D -> 2D (idx2ij)

    :param a:
    :param slo_count:
    :param slo_idx2i:
    :param slo_idx2j:
    :param a_idx:

    :return: updated a
    """
    #real(8) a_idx(slo_count), a(ny, nx)
    a.fill(0.0)
    for k in range( slo_count ):
        a[slo_idx2i[k], slo_idx2j[k]] = a_idx[k]
    #enddo

    return(a)

#end def sub_slo_idx2ij


def sub_slo_idx2ij4( a, slo_count, slo_idx2i, slo_idx2j, a_idx ):
    """
    2D -> 1D (ij2idx)

    :param a: input array 
    :param slo_count:
    :param slo_idx2i:
    :param slo_idx2j:
    :param a_idx:

    :return: Updated a
    """
    #real(8) a_idx(i4, slo_count), a(i4, ny, nx)
    a.fill(0.0)
    for i in range( 4 ):
        for k in range( slo_count ):
            a[i, slo_idx2i[k], slo_idx2j[k]] = a_idx[i, k]
        #enddo
    #enddo
    
    return(a)

#end def sub_slo_idx2ij4


def storage_calc(hs, hr, hg, domain, area, riv_thresh, riv, gampt_ff, gammag_idx, slo_ij2idx):
    """
    # storage calculation

    :param hs:
    :param hr:
    :param hg:
    :param domain:
    :param aera:
    :param riv_thresh:
    :param riv:
    :param gampt_ff:
    :param gammag_idx:
    :param slo_ij2idx:

    :return: ss
    :return: sr
    :return: si
    :return: sg
    """
    ss = 0.0
    sr = 0.0
    si = 0.0
    sg = 0.0
    for i in range( ny ):
        for j in range( nx ):
            if( domain[i,j] == 0 ): 
                continue 
            ss = ss + hs[i,j] * area
            #if(riv_thresh.ge.0 .and. riv(i,j).eq.1) sr = sr + hr(i,j) * area * area_ratio(i,j)
            # modified v1.4
            if( riv_thresh == 0 and riv[i,j] == 1 ):
                vr_temp = hr2vr(hr[i, j], riv_ij2idx[i,j])
                sr = sr + vr_temp
            #endif
            si = si + gampt_ff[i,j] * area
            sg = sg - hg[i,j] * gammag_idx[slo_ij2idx[i,j]] * area # storage deficit
        #enddo
    #enddo
    return(ss, sr, si, sg)
#end def storage_calc



def int2char( num ):
    """
    numbers to characters
    
    Temporary wrapper
    """
    return(str(num))

#end def int2char



def hubeny_sub( x1_deg, y1_deg, x2_deg, y2_deg ):
    """
    Hubeny_sub.f90

    :param x1_deg:
    :param x2_deg:
    :param x3_deg:
    :param x4_deg:

    :return: d
    """
    #real(8) x1_deg, y1_deg, x2_deg, y2_deg
    #real(8) x1, y1, x2, y2
    #real(8) pi, dx, dy, mu, a, b, e, W, N, M, d
    pi = 3.1415926535897
    x1 = x1_deg * pi / 180.0
    y1 = y1_deg * pi / 180.0
    x2 = x2_deg * pi / 180.0
    y2 = y2_deg * pi / 180.0

    dy = y1 - y2
    dx = x1 - x2
    mu = (y1 + y2) / 2.

    a = 6378137.0000 # Semi-Major Axis
    b = 6356752.3140 # Semi-Minor Axis

    e = np.sqrt((a**2.0 - b**2.0) / (a**2.0))
    W = np.sqrt(1. - e**2.0 * (np.sin(mu))**2.0)
    N = a / W
    M = a * (1. - e ** 2.0) / W**(3.0)

    d = np.sqrt((dy * M) ** 2.0 + (dx * N * np.cos(mu)) ** 2.0)
    
    return(d)

#end def

