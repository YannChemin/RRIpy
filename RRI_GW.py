from globals_vars import *

# RRI_GW

def funcg(slo_count, hg_idx, fg_idx, qg_idx):
    """
    :Description: variable definition (slope)

    :param slo_count:
    """
    #real(8) hg_idx(slo_count), fg_idx(slo_count)
    #real(8) qg_idx(i4, slo_count)
    fg_idx = np.zeros(slo_count)
    qg_idx = np.zeros(slo_count)
    qg_idx = qg_calc(hg_idx)
    # qg_idx > 0 --> discharge flowing out from a cell

    #$omp parallel do
    for k in range(1, slo_count):
        if(gammag_idx[k] > 0.0):
            fg_idx[k] = (qg_idx[0,k] + qg_idx[1,k] + qg_idx[2,k] + qg_idx[3,k]) / gammag_idx[k]
    #enddo
    #$omp end parallel do

    for k in range(1, slo_count):
        for l in range(1, lmax):
            if( dif_slo_idx[k] == 0 and l == 2 ):
                exit # kinematic -> 1-direction
            kk = down_slo_idx[l, k]
            if( dif_slo_idx[k] == 0 ):
                kk = down_slo_1d_idx[k]
            if( kk == -1 ):
                continue
            if(gammag_idx[k] > 0.0):
                fg_idx[kk] = fg_idx[kk] - qg_idx[l, k] / gammag_idx[k]
        #enddo
    #enddo

#end def funcg


def qg_calc(slo_count, hg_idx, qg_idx):
    """
    :Description: lateral gw discharge (slope)

    :param slo_count:
    """
    #real(8) hg_idx(slo_count)
    #real(8) qg_idx(i4, slo_count), qg
    #real(8) zb_p, hg_p, gammag_p, kg0_p, ksg_p, fpg_p
    #real(8) zb_n, hg_n, gammag_n, kg0_n, ksg_n, fpg_n
    #real(8) dh, distance
    #real(8) len, hw
    #integer dif_p, dif_n

    qg_idx = 0.0

    #$omp parallel do private(kk,zb_p,hg_p,gammag_p,kg0_p,ksg_p,fpg_p,l,distance,len,zb_n,hg_n,gammag_n,kg0_n,ksg_n,fpg_n,dh,dif_p,dif_n,qg)
    for k in range(1, slo_count):
        zb_p = zb_slo_idx[k]
        hg_p = hg_idx[k]
        ksg_p = ksg_idx[k]
        gammag_p = gammag_idx[k]
        kg0_p = kg0_idx[k]
        fpg_p = fpg_idx[k]
        dif_p = dif_slo_idx[k]
        if( ksg_p <= 0.0 ):
            continue

        # 8-direction: lmax = 4, 4-direction: lmax = 2
        for l in range(1, lmax): # (1: rightC2: down, 3: right down, 4: left down)
            if( dif_p == 0 and l == 2 ):
                exit # kinematic -> 1-direction
            kk = down_slo_idx[l, k]
            if( dif_p == 0 ):
                kk = down_slo_1d_idx[k]
            if( kk == -1 ):
                continue

            distance = dis_slo_idx[l, k]
            len = len_slo_idx[l, k]
            if( dif_p == 0 ):
                distance = dis_slo_1d_idx[k]
            if( dif_p == 0 ):
                len = len_slo_1d_idx[k]

            # information of the destination cell
            zb_n = zb_slo_idx[kk]
            hg_n = hg_idx[kk]
            gammag_n = gammag_idx[kk]
            kg0_n = kg0_idx[kk]
            ksg_n = ksg_idx[kk]
            fpg_n = fpg_idx[kk]
            dif_n = dif_slo_idx[kk]
            if( ksg_n <= 0.0 ):
                continue

            # diffusion wave
            dh = ((zb_p - hg_p) - (zb_n - hg_n)) / distance

            # 1-direction : kinematic wave
            if( dif_p == 0 ):
                dh = max( (zb_p - zb_n) / distance, 0.001 )

            # water coming in or going out?
            if( dh >= 0.0 ):
                qg = hg_calc(gammag_p, kg0_p, ksg_p, fpg_p, hg_p, dh, len, area)
                qg_idx[l,k] = qg
            else:
                # coming in
                dh = abs(dh)
                qg = hg_calc(gammag_n, kg0_n, ksg_n, fpg_n, hg_n, dh, len, area)
                qg_idx[l,k] = -qg
            #endif
    #enddo
#enddo
#$omp end parallel do

#end def qg_calc


def hg_calc(gammag_p, kg0_p, ksg_p, fpg_p, hg_p, dh, length, area)
    """
    :Description: water depth and gw discharge relationship
    
    :param gammag_p:
    :param kg0_p:
    :param ksg_p:
    :param fpg_p:
    :param hg_p:
    :param dh:
    :param length:
    :return: qg the discharge per unit area
    """
    #real(8) gammag_p, kg0_p, ksg_p, fpg_p, hg_p, dh, len, qg
    qg = dh * kg0_p / fpg_p * exp( - fpg_p * hg_p )
    # discharge per unit area
    # (qg multiply by width and divide by area)
    qg = qg * length / area
    return(qg)
#end def hg_calc


def gw_recharge(hs_idx, gampt_ff_idx, hg_idx):
    """
    # gw recharge
    """
    #real(8) hs_idx(slo_count), gampt_ff_idx(slo_count), hg_idx(slo_count)
    rech = 0.0
    rech_ga = 0.0
    rech_hs = 0.0
    for k in range(1, slo_count):
        rech = ksg_idx[k] * dt
        if( hg_idx[k] < 0.0 ):
            continue

        if( hg_idx[k] * gammag_idx[k] < rech ):
            rech = hg_idx[k] * gammag_idx[k]

        if( rech < gampt_ff_idx[k] and ksv_idx[k] < 0.0 ):
            rech_ga = rech
            rech_hs = 0.0
        else:
            if( ksv_idx[k] > 0.0):
                rech_ga = gampt_ff_idx[k]

            if( rech - rech_ga .lt. hs_idx[k] ):
                rech_hs = rech - rech_ga
            else:
                rech_hs = hs_idx[k]
            #endif
        #endif
        hs_idx[k] = hs_idx[k] - rech_hs
        if( ksv_idx[k] > 0.0 ):
            gampt_ff_idx[k] = gampt_ff_idx[k] - rech_ga

        rech = rech_ga + rech_hs
        if(gammag_idx[k] > 0.0):
            hg_idx[k] = hg_idx[k] - rech / gammag_idx[k]

        rech_hs_tsas[k] = rech_hs / float(dt)
    #enddo

#end def gw_recharge

def gw_lose(hg_idx):
    """
    # gw lose
    """
    #real(8) hg_idx(slo_count)
    for k in range(1, slo_count):
        if( rgl_idx[k] > 0.0):
            hg_idx[k] = hg_idx[k] + rgl_idx[k] / gammag_idx[k] * dt
    #enddo

#end def gw_lose


def gw_exfilt(hs_idx, gampt_ff_idx, hg_idx)
    """
    # gw exfilt
    """
    #real(8) hs_idx(slo_count), gampt_ff_idx(slo_count), hg_idx(slo_count)
    #real(8) exfilt, exfilt_ga, exfilt_hs
    for k in range(1, slo_count):
        exfilt = 0.0
        exfilt_ga = 0.0
        exfilt_hs = 0.0
        if( hg_idx[k] >= 0.0 ):
            continue

        exfilt = - hg_idx[k] * gammag_idx[k]
        if( infilt_limit_idx[k] > gampt_ff_idx[k] and infilt_limit_idx[k] - gampt_ff_idx[k] >= exfilt and ksv_idx[k] > 0.0 and infilt_limit_idx[k] > 0.0 ):
            exfilt_ga = exfilt
            exfilt_hs = 0.0
        else
            if( ksv_idx[k] > 0.0 and infilt_limit_idx[k] > 0.0 ):
                exfilt_ga = infilt_limit_idx[k] - gampt_ff_idx[k]
            #endif
            exfilt_hs = exfilt - exfilt_ga
        #endif
        hg_idx[k] = 0.0
        if( ksg_idx[k] > 0.0 and infilt_limit_idx[k] > 0.0 ):
            gampt_ff_idx[k] = gampt_ff_idx[k] + exfilt_ga
        
        hs_idx[k] = hs_idx[k] + exfilt_hs
        exfilt_hs_tsas[k] = exfilt_hs / float(dt)
    #enddo
#end def gw_exfilt


# initial setting for hg
def hg_init(slo_count):
    return(np.zeros(slo_count))
#end def hg_init
