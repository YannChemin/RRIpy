# RRI_Riv.f90

def funcr( riv_count, vr_idx, hr_idx, bound_riv_wlev_switch, tt_max_bound_riv_wlev, t_bound_riv_wlev, time, ddt, qr_div_idx ):
    """
    Variable definition (river)

    :param riv_count:
    :param vr_idx:
    :param hr_idx:
    :param bound_riv_wlev_switch: 
    :param tt_max_bound_riv_wlev:
    :param t_bound_riv_wlev:
    :param time:
    :param ddt:
    :param qr_div_idx:


    :return: fr_idx
    :return: qr_idx
    """
    #real(8) hr_idx(riv_count), vr_idx(riv_count), fr_idx(riv_count), qr_idx(riv_count)
    #real(8) qr_sum_idx(riv_count), qr_div_idx(riv_count)

    fr_idx = np.zeros(riv_count)
    qr_idx = np.zeros(riv_count)
    qr_sum_idx = np.zeros(riv_count)
    qr_div_idx = np.zeros(riv_count)

    # add v1.4
    for k in range( riv_count ):
        hr_idx[k] = vr2hr(vr_idx[k], k) 
    #enddo

    # boundary condition for river (water depth boundary)
    if( bound_riv_wlev_switch >= 1 ):
        itemp = -1
        for jtemp in range( tt_max_bound_riv_wlev ):
            if( t_bound_riv_wlev(jtemp-1) < (time + ddt) and (time + ddt) <= t_bound_riv_wlev(jtemp) ):
                itemp = jtemp
        #enddo
        for k in range( riv_count ):
            if( bound_riv_wlev_idx(itemp, k) <= -100.0 ):
                continue # not boundary
            hr_idx[k] = bound_riv_wlev_idx(itemp, k)
            vr_idx[k] = hr2vr(hr_idx[k], k) # add v1.4
        #enddo
    #endif

    qr_idx = qr_calc(hr_idx)
    if( div_switch == 1 ):
        call RRI_Div(qr_idx, hr_idx, qr_div_idx)

    # dam control
    if( dam_switch == 1 ):
        call dam_prepare(qr_idx) # calculate inflow to dam
        for i in range( dam_num ):
            if( dam_volmax[i] > 0.0 ):
                # dam
                call dam_operation( dam_loc[i] )
                qr_idx[ int(dam_loc[i]) ] = dam_qout[i]
                qr_sum_idx[ int(dam_loc[i]) ] = qr_sum_idx[ int(dam_loc[i]) ] + dam_qin[ int(dam_loc[i]) ] - dam_qout[i]
            else if( dam_volmax[i] == 0.0 ):
                # barrage
                if( hr_idx[ int(dam_loc[i]) ] <= dam_floodq[i] ):
                    qr_idx[ int(dam_loc[i]) ] = 0.0
                #endif
            else:
                # water gate (dam_volmax[i] < 0) # added on Aug 7, 2021 by T.Sayama
                k = dam_loc[i]
                kk = down_riv_idx[k]
                call gate_operation( dam_loc[i], hr_idx[k], hr_idx[kk] )
                qr_idx( int(dam_loc[i]) ) = dam_qout[i]
                #qr_sum_idx( dam_loc[i] ) = qr_sum_idx( dam_loc[i] ) + dam_qin( dam_loc[i] ) - dam_qout[i]
            #endif
        #enddo
    #endif

    # boundary condition for river (discharge boundary)
    if( bound_riv_disc_switch >= 1 ):
        for jtemp in range( tt_max_bound_riv_disc ):
            if( t_bound_riv_disc(jtemp-1) < (time + ddt) and (time + ddt) <= t_bound_riv_disc(jtemp) ):
                itemp = jtemp
        #enddo
        for k in range( riv_count ):
            if( bound_riv_disc_idx[itemp, k] <= -100.0 ):
                continue # not boundary
            #qr_idx[k] = bound_riv_disc_idx[itemp, k] / area  # qr_idx: discharge per unit area
            qr_idx[k] = bound_riv_disc_idx[itemp, k] # modified v1.4
            # linear interpolation of the boundary condition
            #qr_idx[k] = bound_riv_disc_idx(itemp-1, k) * (t_bound_riv_disc(itemp) - (time + ddt)) &
            #           + bound_riv_disc_idx(itemp, k) * ((time + ddt) - t_bound_riv_disc(itemp-1))
            #qr_idx[k] = qr_idx[k] / (t_bound_riv_disc(itemp) - t_bound_riv_disc(itemp-1))
            hr_idx[k] = 0.0
            vr_idx[k] = 0.0 # add v1.4
        #enddo
    #endif

    # qr_sum > 0 --> discharge flowing out from a cell
    for k in range( riv_count ):
        # outflow from [k]
        qr_sum_idx[k] = qr_sum_idx[k] + qr_idx[k]
        kk = down_riv_idx[k]
        if(domain_riv_idx[kk]==0):
            continue
        # qr_sum minus (flowing into) discharge at the destination cell
        qr_sum_idx[kk] = qr_sum_idx[kk] - qr_idx[k]
    #enddo

    # diversion
    if( div_switch == 1 ):
        for l in range( div_id_max ):
            # outflow from [k]
            k = div_org_idx(l)
            kk = down_riv_idx[k]
            if( div_rate(l) <= 0.0 ): # modified by T.Sayama on Dec. 7, 2022 v1.4.2.7
                qr_sum_idx[k] = qr_sum_idx[k] + qr_div_idx[k]
            else:
                qr_sum_idx[kk] = qr_sum_idx[kk] + qr_div_idx[k]
            #endif
            kk = div_dest_idx(l)
            if(domain_riv_idx[kk]==0):
                continue
            # qr_sum minus (flowing into) discharge at the destination cell
            qr_sum_idx[kk] = qr_sum_idx[kk] - qr_div_idx[k] 
        #enddo
    #endif

    # qr_sum divide by area_ratio
    # (area_ratio = ratio of river area against total cell area) 
    #fr_idx = -qr_sum_idx / area_ratio_idx
    fr_idx = - qr_sum_idx # modified v1.4
    return(fr_idx, qr_idx)
#enddef funcr



def qr_calc(hr_idx, qr_idx):
    """
    Lateral discharge (river)

    :param hr_idx:

    :return: qr_idx
    """
    #real(8) hr_idx(riv_count), qr_idx(riv_count), qr_div_idx(riv_count)
    #real(8) zb_p, hr_p
    #real(8) zb_n, hr_n
    #real(8) dh, distance
    #real(8) qr_temp, hw

    qr_idx(:) = 0.0
    qr_div_idx(:) = 0.0

    #$omp parallel do private(kk,zb_p,hr_p,distance,zb_n,hr_n,dh,hw,qr_temp,dif_p,dif_n)
    for k in range( riv_count):
        if(domain_riv_idx[k] == 2):
            continue
        zb_p = zb_riv_idx[k]
        hr_p = hr_idx[k]
        dif_p = dif_riv_idx[k]
        distance = dis_riv_idx[k]
        # information of the destination cell
        kk = down_riv_idx[k]
        zb_n = zb_riv_idx[kk]
        hr_n = hr_idx[kk]
        dif_n = dif_riv_idx[kk]
        # diffusion wave
        dh = ((zb_p + hr_p) - (zb_n + hr_n)) / distance # diffussion
        # kinematic wave
        if( dif_p == 0 ):
            dh = max( (zb_p - zb_n) / distance, 0.001 )
        # the destination cell is outlet (domain = 2)
        if( domain_riv_idx[kk] == 2 ):
            dh = (zb_p + hr_p - zb_n) / distance # kinematic wave (+hr_p)
        # the destination cell is outlet (domain = 2 ) with water depth boundary water table
        if( domain_riv_idx[kk] == 2 and bound_riv_wlev_switch >= 1 ):
        # ver 1.4.2 mod by T.Sayama on June 24, 2015
        #if( bound_riv_wlev_idx(1, kk) <= -100.0 ) continue # not boundary
        #dh = ((zb_p + hr_p) - (zb_n + hr_n)) / distance # diffussion
            if( bound_riv_wlev_idx(1, kk) > -100.0 ):
                dh = ((zb_p + hr_p) - (zb_n + hr_n)) / distance # diffussion 
        #endif
        # kinematic wave (for dam and water gate) modified by T.Sayama on Aug 9, 2021
        if ( damflg[k] > 0 ):
            dh = max((zb_p - zb_n) / distance, 0.001)
        #endif
        if ( damflg[kk] > 0 ):
            if( dam_floodq(damflg[kk]) > 0 ):
                dh = max((zb_p - zb_n) / distance, 0.001)
            #endif
        #endif
        # from a tributary (levee height =< 0) to a main river (levee height > 0) : no reverse flow
        #if( height_idx[k] <= 0.0 and height_idx[kk] > 0.0 ) dh = max( (zb_p - zb_n) / distance, 0.001 )

        if( dh >= 0.0 ):
            hw = hr_p
            if( zb_p < zb_n ):
                hw = max(0.0, zb_p + hr_p - zb_n)
            #call hq_riv(hw, dh, width_idx[k], qr_temp)
            call hq_riv(hw, dh, k, width_idx[k], qr_temp)
            qr_idx[k] = qr_temp
        else:
            # reverse flow
            hw = hr_n
            if( zb_n < zb_p ):
                hw = max(0.0, zb_n + hr_n - zb_p)
            dh = abs(dh)
            #call hq_riv(hw, dh, width_idx[k], qr_temp)
            call hq_riv(hw, dh, kk, width_idx[k], qr_temp)
            qr_idx[k] = -qr_temp
        #endif
    #enddo
    #$omp end parallel do
    return(qr_idx)

#enddef qr_calc


def hq_riv(dh, ns_river, w, h, k):
    """
    Water depth and discharge relationship

    :param dh:
    :param ns_river:
    :param w:
    :param h:
    :param k:

    :return: q
    """
    a = sqrt(abs(dh)) / ns_river
    #m = 5.0 / 3.0
    r = (w * h) / (w + 2.0 * h)
    # modified v1.4
    #q = a * h ** m * w
    # modified v1.4.2.4
    q = a * r ** (2.0 / 3.0) * w * h

    if( sec_map_idx[k] > 0 ):
        call sec_hq_riv(h, dh, k, q)
    #endif

    # discharge per unit area
    # (q multiply by width and divide by area)
    #q = q * w / area
    return(q)
#enddef hq_riv
