from globals_vars import *

def infilt(hs_idx, gampt_ff_idx, gampt_f_idx):
    """
    :Description: infiltration

    :param slo_count:
    :param gampt_f_idx: infiltration capacity [m/s]
    :param gampt_ff_idx: accumulated infiltration depth [m]
    :param ksv_idx: 
    :param faif_idx:
    :param gammaa_idx:
    :param hs_idx:
    :return: hs_idx
    :return: gampt_ff_idx the accumulated infiltration depth [m]
    :return: gampt_f_idx the infiltration capacity [m/s]
    """
    #real(8) hs_idx(slo_count), gampt_ff_idx(slo_count), gampt_f_idx(slo_count)
    #real(8) gampt_ff_temp

    for k in range(1, slo_count):
        gampt_f_idx[k] = 0.0
        gampt_ff_temp = gampt_ff_idx[k]
        if( gampt_ff_temp <= 0.01 ):
            gampt_ff_temp = 0.01

        # gampt_f_idx[k] : infiltration capacity [m/s]
        # gampt_ff : accumulated infiltration depth [m]
        gampt_f_idx[k] = ksv_idx[k] * (1.0 + faif_idx[k] * gammaa_idx[k] / gampt_ff_temp)

        # gampt_f_idx[k] : infiltration capacity -> infiltration rate [m/s]
        if( gampt_f_idx[k] >= hs_idx[k] / dt ):
            gampt_f_idx[k] = hs_idx[k] / dt

        # gampt_ff should not exceeds a certain level
        if( infilt_limit_idx[k] >= 0.0 and gampt_ff_idx[k] >= infilt_limit_idx[k] ):
            gampt_f_idx[k] = 0.0

        # update gampt_ff [m]
        gampt_ff_idx[k] = gampt_ff_idx[k] + gampt_f_idx[k] * dt

        # hs : hs - infiltration rate * dt [m]
        hs_idx[k] = hs_idx[k] - gampt_f_idx[k] * dt
        if( hs_idx[k] <= 0.0 ):
            hs_idx[k] = 0.0
    #enddo
    return(hs_idx, gampt_ff_idx, gampt_f_idx)
#end def infilt
