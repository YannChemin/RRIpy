# RRI_Slope

def funcs(hs_idx, qp_t_idx, fs_idx, qs_idx ):
    """
    Variable definition (slope)

    :param :
    :param :
    :param :
    :param :
    """
    #real(8) hs_idx(slo_count), qp_t_idx(slo_count), fs_idx(slo_count)
    #real(8) qs_idx(i4, slo_count)
    fs_idx.fill(0.0)
    qs_idx.fill(0.0)
    qs_idx = qs_calc(hs_idx)
    # boundary condition for slope (discharge boundary)
    if( bound_slo_disc_switch >= 1 ):
        #itemp = time / dt_bound_slo + 1
        itemp = -1
        for jtemp in range( tt_max_bound_slo_disc ):
            if( t_bound_slo_disc[jtemp-1] < (time + ddt) and (time + ddt) <= t_bound_slo_disc[jtemp] ):
                itemp = jtemp
        #enddo
        for k in range( slo_count ):
            if( bound_slo_disc_idx[itemp,k] <= -100.0 ):
                continue # not boundary
            # right
            if( direc(slo_idx2i[k], slo_idx2j[k])==1 ):
                qs_idx[1,k] = bound_slo_disc_idx[itemp,k] / area
            # right down
            else if( direc(slo_idx2i[k], slo_idx2j[k])==2 ):
                qs_idx[3,k] = bound_slo_disc_idx[itemp,k] / area
            # down
            else if( direc(slo_idx2i[k], slo_idx2j[k])==4 ):
                qs_idx[2,k] = bound_slo_disc_idx[itemp,k] / area
            # left down
            else if( direc(slo_idx2i[k], slo_idx2j[k])==8 ):
                qs_idx[4,k] = bound_slo_disc_idx[itemp,k] / area
            # left
            else if( direc(slo_idx2i[k], slo_idx2j[k])==16 ):
                qs_idx[1,k] = - bound_slo_disc_idx[itemp,k] / area
            # left up
            else if( direc(slo_idx2i[k], slo_idx2j[k])==32 ):
                qs_idx[3,k] = - bound_slo_disc_idx[itemp,k] / area
            # up
            else if( direc(slo_idx2i[k], slo_idx2j[k])==64 ):
                qs_idx[2,k] = - bound_slo_disc_idx[itemp,k] / area
            # right up
            else if( direc(slo_idx2i[k], slo_idx2j[k])==128 ):
                qs_idx[4,k] = - bound_slo_disc_idx[itemp,k] / area
            #endif
        #enddo
    #endif
    # qs_idx > 0 --> discharge flowing out from a cell

    #$omp parallel do
    do k = 1, slo_count
        fs_idx[k] = qp_t_idx[k] - (qs_idx(1,k) + qs_idx(2,k) + qs_idx(3,k) + qs_idx(4,k))
    #enddo
    #$omp end parallel do

    for k in range( slo_count ):
        for l in range( lmax ):
            if( dif_slo_idx[k] == 0 and l == 2 ):
                break # kinematic -> 1-direction
            kk = down_slo_idx[l, k]
            if( dif_slo_idx[k] == 0 ):
                kk = down_slo_1d_idx[k]
            if( kk == -1 ):
                continue
            fs_idx[kk] = fs_idx[kk] + qs_idx[l, k]
        #enddo
    #enddo
    return(fs_idx)
#enddef funcs


def qs_calc(slo_count, zb_slo_idx, hs_idx, ns_slo_idx, ka_idx, da_idx, dm_idx, beta_idx, dif_slo_idx, lmax, down_slo_idx, down_slo_1d_idx, dis_slo_idx, dis_slo_1d_idx, len_slo_idx, len_slo_1d_idx, lev_p, lev_n, area, qs_idx):
    """
    Lateral discharge (slope)

    :param slo_count: 
    :param zb_slo_idx:
    :param hs_idx:
    :param ns_slo_idx:
    :param ka_idx:
    :param da_idx:
    :param dm_idx:
    :param beta_idx:
    :param dif_slo_idx:
    :param lmax:
    :param down_slo_idx:
    :param down_slo_1d_idx:
    :param dis_slo_idx:
    :param dis_slo_1d_idx:
    :param len_slo_idx:
    :param len_slo_1d_idx:
    :param lev_p:
    :param lev_n:
    :param area:
    :param qs_idx:


    :return: qs_idx 
    """
    #$omp parallel do private(kk,zb_p,hs_p,ns_p,ka_p,da_p,dm_p,b_p,dif_p,l,distance,len, &
    #$omp                     zb_n,hs_n,ns_n,ka_n,da_n,dm_n,b_n,dif_n,lev_p,lev_n,dh,hw,q)
    for k in range( slo_count ):
        zb_p = zb_slo_idx[k]
        hs_p = hs_idx[k]
        ns_p = ns_slo_idx[k]
        ka_p = ka_idx[k]
        da_p = da_idx[k]
        dm_p = dm_idx[k]
        b_p  = beta_idx[k]
        dif_p = dif_slo_idx[k]
        # 8-direction: lmax = 4, 4-direction: lmax = 2
        for l in range( lmax ): # (1: right, 2: down, 3: right down, 4: left down)
            if( dif_p == 0 and l == 2 ):
                break # kinematic -> 1-direction
            kk = down_slo_idx[l,k]
            if( dif_p == 0 ):
                kk = down_slo_1d_idx[k]
            if( kk == -1 ):
                continue
            distance = dis_slo_idx[l,k]
            length = len_slo_idx[l,k]
            if( dif_p == 0 ):
                distance = dis_slo_1d_idx[k]
            if( dif_p == 0 ):
                length = len_slo_1d_idx[k]

            # information of the destination cell
            zb_n = zb_slo_idx[kk]
            hs_n = hs_idx[kk]
            ns_n = ns_slo_idx[kk]
            ka_n = ka_idx[kk]
            da_n = da_idx[kk]
            dm_n = dm_idx[kk]
            b_n = beta_idx[kk]
            dif_n = dif_slo_idx[kk]

            call h2lev(hs_p, k, lev_p)
            call h2lev(hs_n, kk, lev_n)

            # diffusion wave
            dh = ((zb_p + lev_p) - (zb_n + lev_n)) / distance

            # 1-direction : kinematic wave
            if( dif_p == 0 ):
                dh = max( (zb_p - zb_n) / distance, 0.001 )

            # embankment
            #if(emb_switch==1):
            # if(l==1) emb = emb_r_idx[k]
            # if(l==2) emb = emb_b_idx[k]
            # if(l==3) emb = max( emb_r_idx[k], emb_b_idx[k] )
            # if(l==4) emb = max( emb_r_idx[kk], emb_b_idx[k] )
            ##endif

            # water coming in or going out?
            if( dh >= 0.0 ):
                # going out
                hw = hs_p
                #if(emb > 0.0) hw = max(hs_p - emb, 0.0)
                if( zb_p < zb_n ):
                    hw = max(0.0, zb_p + hs_p - zb_n)
                    q = hq(ns_p, ka_p, da_p, dm_p, b_p, hw, dh, length, area)
                    qs_idx(l,k) = q
                else:
                    # coming in
                    hw = hs_n
                    #if(emb > 0.0) hw = max(hs_n - emb, 0.0)
                    dh = abs(dh)
                    if( zb_n < zb_p ):
                        hw = max(0.0, zb_n + hs_n - zb_p)
                        q = hq(ns_n, ka_n, da_n, dm_n, b_n, hw, dh, length, area)
                        qs_idx(l,k) = -q
            #endif
        #enddo
    #enddo
    return(qs_idx)
    #$omp end parallel do

#enddef qs_calc

def hq(ns_p, ka_p, da_p, dm_p, b_p, h, dh, length, area):
    """
    Water depth and discharge relationship

    :param b_p:
    :param ka_p:
    :param dh:
    :param da_p:
    :param ns_p:
    :param h:
    :param dm_p:
    :param length:
    :param area:

    :return: q
    """
    #real(8) ns_p, da_p, dm_p, ka_p, b_p, h, dh, len, q
    #real(8) km, vm, va, al, m
    if( b_p > 0.0 ):
        km = ka_p / b_p
    else:
        km = 0.0
    #endif
    vm = km * dh

    if( da_p > 0.0 ):
        va = ka_p * dh
    else:
        va = 0.0
    #endif

    if( dh < 0 ):
        dh = 0.0
        al = sqrt(dh) / ns_p
        m = 5.0 / 3.0

    if( h < dm_p ):
        q = vm * dm_p * (h / dm_p) ** b_p
    else if( h < da_p ):
        q = vm * dm_p + va * (h - dm_p)
    else:
        q = vm * dm_p + va * (h - dm_p) + al * (h - da_p) ** m
    #endif

    # discharge per unit area
    # (q multiply by width and divide by area)
    q = q * length / area

    # water depth limitter (1 mm)
    # note: it can be set to zero
    #if( h<=0.001 ) q = 0.0
    
    return(q)

#enddef hq


def h2lev(h, k, lev)
    """
    Water depth (h) to actual water level (lev)
    """
    #real(8) h, lev
    #real(8) rho
    #real(8) da_temp

    da_temp = soildepth_idx[k] * gammaa_idx[k]

    if( soildepth_idx[k] == 0.0 ):
        lev = h
    else if( h >= da_temp ): # including da = 0
        lev = soildepth_idx[k] + (h - da_temp) # surface water
    else:
        if(soildepth_idx[k] > 0.0 ):
            rho = da_temp / soildepth_idx[k]
        lev = h / rho
#endif

#enddef h2lev
