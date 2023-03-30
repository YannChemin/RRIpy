from globals_vars import *

# Evapotranspiration
def evp( evp_switch, itemp, jtemp, t_evp, evp_i, evp_j, qe_t, slo_count, hs_idx, gampt_ff_idx, aevp, aevp_tsas, aevp_sum, tevp_sum):
    """
    :Description: RRI_Evp 

    :param evp_switch: evp_switch = 1 : Allow ET from gampt_ff_idx = 2 : Do not take ET from gampt_ff_idx
    
    :param itemp: ???
    :param jtemp: ???
    :param t_evp: ??? seems to be an evap temporal index array
    :param evp_i: ???
    :param evp_j: ???
    :param qe_t:  ???
    :param slo_count: ???
    :param hs_idx:
    :param gampt_ff_idx:
    :return: aevp
    :return: aevp_tsas
    :return: hs_idx
    :return: gampt_ff_idx
    :return: aevp_sum (Added by Yann)
    :return: pevp_sum (Added by Yann)
    """
    #real(8) hs_idx(slo_count), gampt_ff_idx(slo_count), qe_t_temp

    itemp = -1
    for jtemp in range(1, tt_max_evp):
        if( t_evp[jtemp-1] < time and time <= t_evp[jtemp] ):
            itemp = jtemp
    #enddo
    for i in range(1, ny):
        for j in range(1, nx):
            qe_t[i, j] = qe(itemp, evp_i[i], evp_j[j])
        #enddo
    #enddo
    qe_t_idx = sub_slo_ij2idx(qe_t)

    for k in range(1, slo_count):
        i = slo_idx2i[k]
        j = slo_idx2j[k]
        qe_t_temp = qe_t_idx[k]
 
    # from ver 1.3.3
    #if( dm_idx[k] .gt. 0.0 .and. hs_idx[k] .lt. dm_idx[k] ):
    # qe_t_temp = qe_t_idx[k] * hs_idx[k] / dm_idx[k]
    ##endif

    if( evp_switch == 1 ):
        aevp[i, j] = min( qe_t_temp, (hs_idx[k] + gampt_ff_idx[k]) / dt )
    else:
        aevp[i, j] = min( qe_t_temp, hs_idx[k] / dt )
    #endif
    aevp_tsas[k] = min( qe_t_temp, hs_idx[k] / dt )

    hs_idx[k] = hs_idx[k] - aevp[i, j] * dt
    if( hs_idx[k] < 0.0 ):
        if( evp_switch == 1 ):
            gampt_ff_idx[k] = gampt_ff_idx[k] + hs_idx[k]

            hs_idx[k] = 0.0
            if( gampt_ff_idx[k] < 0.0 ):
                gampt_ff_idx[k] = 0.0
            #endif
        #endif
        aevp_sum = aevp_sum + aevp[i, j] * dt * area
        pevp_sum = pevp_sum + qe_t_idx[k] * dt * area
    #enddo
    return(aevp, aevp_tsas, hs_idx, gampt_ff_idx, aevp_sum, pevp_sum)
#end def evp
