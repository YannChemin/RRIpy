# RRI_Sub

# river index setting
def riv_idx_setting():
    """
    River index setting
    """
    #use globals
    #implicit none
    #integer i, j, ii, jj, k, kk # add v1.4(k, kk)
    #real(8) distance

    riv_count = 0
    for i in range(ny):
        for j in range (nx):
            if( domain[i,j] != 0 and riv[i,j] == 1):
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
    for i in range(ny):
        for j in range(nx):
            if(domain[i,j] == 0 or riv[i,j] != 1):
                continue
            # domain[i, j] = 1 or 2 and riv(i, j) = 1
            riv_count = riv_count + 1
            riv_idx2i[riv_count] = i
            riv_idx2j[riv_count] = j
            domain_riv_idx[riv_count] = domain[i, j]
            width_idx[riv_count] = width[i, j]
            depth_idx[riv_count] = depth[i, j]
            height_idx[riv_count] = height[i, j]
            area_ratio_idx[riv_count] = area_ratio[i, j]
            zb_riv_idx[riv_count] = zb_riv[i, j]
            riv_ij2idx[i, j] = riv_count
            dif_riv_idx[riv_count] = dif[land[i, j]]
            sec_map_idx[riv_count] = sec_map[i, j]
            len_riv_idx[riv_count] = len_riv[i, j] # add v1.4
        #enddo
    #enddo

    # search for downstream gridcell (down_idx)
    riv_count = 0
    for i in range(ny):
        for j in range(nx):
            if(domain[i,j] == 0 or riv[i,j] == 1):
                continue
            # domain(i, j) = 1 or 2 and riv(i, j) = 1
            riv_count = riv_count + 1
            # right
            if( direc[i,j] == 1 ):
                ii = i
                jj = j + 1
                distance = dx
            # right down
            else if( direc[i,j] == 2 ):
                ii = i + 1
                jj = j + 1
                distance = np.sqrt( dx*dx + dy*dy )
            # down
            else if( direc[i,j] == 4 ):
                ii = i + 1
                jj = j
                distance = dy
            # left down
            else if( direc[i,j] == 8 ):
                ii = i + 1
                jj = j - 1
                distance = np.sqrt( dx*dx + dy*dy )
            # left
            else if( direc[i,j] == 16 ):
                ii = i
                jj = j - 1
                distance = dx
            # left up
            else if( direc[i,j] == 32 ):
                ii = i - 1
                jj = j - 1
                distance = np.sqrt( dx*dx + dy*dy )
            # up
            else if( direc[i,j] == 64 ):
                ii = i - 1
                jj = j
                distance = dy
            # right up
            else if( direc[i,j] ==.128 ):
                ii = i - 1
                jj = j + 1
                distance = np.sqrt( dx*dx + dy*dy )
            else if( direc[i,j] == 0 or direc[i,j] == -1 ):
                ii = i
                jj = j
            else
                raise Exception ( "dir[%d, %d] error %f" % ( i, j, direc[i, j]) )
            #endif

            # If the downstream cell is outside the domain, set domain(i, j) = 2
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
                raise Exception ("riv(%d, %d) should be 1 %d %d" % (i, j, ii, jj) )
            #endif

            dis_riv_idx[riv_count] = distance
            down_riv_idx[riv_count] = riv_ij2idx[ii, jj]

        #enddo
    #enddo

    # add v1.4
    if( sec_length_switch == 1 ):
        for k in range(riv_count): 
            kk = down_riv_idx[k]
            dis_riv_idx[k] = ( len_riv_idx[k] + len_riv_idx[kk] ) / 2.0
        #enddo
    #endif

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


def slo_idx_setting
    """
    Slope index setting


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
            slo_count = slo_count + 1
            slo_idx2i[slo_count] = i
            slo_idx2j[slo_count] = j
            domain_slo_idx[slo_count] = domain[i, j]
            zb_slo_idx[slo_count] = zb[i, j]
            acc_slo_idx[slo_count] = acc[i, j]
            slo_ij2idx[i, j] = slo_count
            land_idx[slo_count] = land[i, j]

            dif_slo_idx[slo_count] = dif[land[i, j]]
            ns_slo_idx[slo_count] = ns_slope[land[i,j]]
            soildepth_idx[slo_count] = soildepth[land[i, j]]
            gammaa_idx[slo_count] = gammaa[land[i, j]]

            ksv_idx[slo_count] = ksv[land[i,j]]
            faif_idx[slo_count] = faif[land[i,j]]
            infilt_limit_idx[slo_count] = infilt_limit[land[i,j]]
  
            ka_idx[slo_count] = ka[land[i,j]]
            gammam_idx[slo_count] = gammam[land[i,j]]
            beta_idx[slo_count] = beta[land[i,j]]
            da_idx[slo_count] = da[land[i,j]]
            dm_idx[slo_count] = dm[land[i,j]]

            ksg_idx[slo_count] = ksg[land[i, j]]
            gammag_idx[slo_count] = gammag[land[i, j]]
            kg0_idx[slo_count] = kg0[land[i, j]]
            fpg_idx[slo_count] = fpg[land[i, j]]
            rgl_idx[slo_count] = rgl[land[i, j]]

        #enddo
    #enddo

    if( eight_dir == 1 ):
        # Hromadka etal, JAIH2006 (USE THIS AS A DEFAULT)
        # 8-direction
        lmax = 4
        l1 = dy / 2.0
        l2 = dx / 2.0
        l3 = np.sqrt(dx ** 2.0 + dy ** 2.0) / 4.0
    else if( eight_dir == 0 ):
        # 4-direction
        lmax = 2
        l1 = dy
        l2 = dx
        l3 = 0.0
    else
        raise Exception ( "error: eight_dir should be 0 or 1.")
    #endif

    # search for downstream gridcell (down_slo_idx)
    slo_count = 0
    down_slo_idx.fill(-1)

    for i in range( ny ):
        do j in range( nx ):
            if(domain[i,j] == 0):
                continue
            # domain(i, j) = 1 or 2
            slo_count = slo_count + 1
            # 8-direction: lmax = 4, 4-direction: lmax = 2
            for l in range( lmax ): # (1: right�C2: down, 3: right down, 4: left down)
                if( l == 1 ):
                    ii = i
                    jj = j + 1
                    distance = dx
                    length = l1
                else if( l == 2 ):
                    ii = i + 1
                    jj = j
                    distance = dy
                    length = l2
                else if( l == 3 ):
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
            slo_count = slo_count + 1
            # right
            if( direc[i,j] == 1 ):
                ii = i
                jj = j + 1
                distance = dx
                length = l1_kin
            # right down
            else if( direc[i,j] == 2 ):
                ii = i + 1
                jj = j + 1
                distance = np.sqrt( dx*dx + dy*dy )
                length = l3_kin
            # down
            else if( direc[i,j] == 4 ):
                ii = i + 1
                jj = j
                distance = dy
                length = l2_kin
            # left down
            else if( direc[i,j] == 8 ):
                ii = i + 1
                jj = j - 1
                distance = np.sqrt( dx*dx + dy*dy )
                length = l3_kin
            # left    
            else if( direc[i,j] == 16 ):
                ii = i
                jj = j - 1
                distance = dx
                length = l1_kin
            # left up
            else if( direc[i,j] == 32 ):
                ii = i - 1
                jj = j - 1
                distance = np.sqrt( dx*dx + dy*dy )
                length = l3_kin
            # up
            else if( direc[i,j] == 64 ):
                ii = i - 1
                jj = j
                distance = dy
                length = l2_kin
            # right up
            else if( direc[i,j] == 128 ):
                ii = i - 1
                jj = j + 1
                distance = np.sqrt( dx*dx + dy*dy )
                length = l3_kin
            else if( direc[i,j] == 0 or direc[i,j] == -1 ):
                ii = i
                jj = j
                distance = dx
                length = l1_kin
            else:
                raise Exception ("dir(%d, %d) error %d %d" %( i, j, direc[i, j])
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

    #enddo
#enddo

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


def storage_calc(hs, hr, hg, domain, area, ):
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
    #real(8) hs(ny, nx), hr(ny, nx), hg(ny, nx)
    #real(8) ss, sr, si, sg
    #real(8) vr_temp # add v1.4

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



def int2char( num )
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

    d = np.sqrt((dy * M) ** 2.0 + (dx * N * cos(mu)) ** 2.0)
    
    return(d)

#end def

