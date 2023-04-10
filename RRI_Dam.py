# RRI_Dam.f90
# modified by T.Sayama on Dec. 2, 2020

def dam_read(riv_count, dam_switch, damfile):
    """
    Reading dam control file

    :param riv_count:
    :param dam_switch:
    :param damfile:

    :return: damflg
    :return: dam_qin
    :return: dam_num
    :return: dam_name
    :return: dam_kind
    :return: dam_ix
    :return: dam_iy
    :return: dam_vol
    :return: dam_vol_temp
    :return: dam_volmax
    :return: dam_state
    :return: dam_qout
    :return: dam_loc
    :return: dam_floodq
    :return: dam_maxfloodq
    :return: dam_rate
    """
    damflg = np.zeros(riv_count)
    dam_qin = np.zeros(riv_count)
    if( dam_switch == 1 ): 
        f99 = open(damfile)
        lines_list = f99.readlines()
        dam_num = int(lines_list[0]) #only one value in the first line
    dam_name = np.zeros(dam_num)
    dam_kind = np.zeros(dam_num)
    dam_ix = np.zeros(dam_num)
    dam_iy = np.zeros(dam_num)
    dam_vol = np.zeros(dam_num)
    dam_vol_temp = np.zeros(dam_num)
    dam_volmax = np.zeros(dam_num)
    dam_state = np.zeros(dam_num)
    dam_qout = np.zeros(dam_num)
    dam_loc = np.zeros(dam_num)
    dam_floodq = np.zeros(dam_num)
    dam_maxfloodq = np.zeros(dam_num)
    dam_rate = np.zeros(dam_num)
    # Scroll through the dam descriptions lines (one dam per line)
    for i in range( len(line_list) - 1 ):
        dam_name[i], dam_iy[i], dam_ix[i], dam_volmax[i], dam_floodq[i], dam_maxfloodq[i], dam_rate[i], dam_vol[i] = lines_list[i+1].split("")
        dam_loc[i] = riv_ij2idx( dam_iy[i], dam_ix[i] )
        damflg[dam_loc[i]] = i
    #enddo
    f99.close()
    
    return(damflg, dam_qin, dam_num, dam_name, dam_kind, dam_ix, dam_iy, dam_vol, dam_vol_temp, dam_volmax, dam_state, dam_qout, dam_loc, dam_floodq, dam_maxfloodq, dam_rate)
#end def dam_read

#######################################################################

def dam_prepare(qr_idx):
    """
    Calculating inflow to dam

    This is a passthrough function...
    """
    # modified by TS on June 16, 2016
    dam_qin[:] = qr_idx[:]
    
    return(qr_idx)
#end def dam_prepare

#######################################################################

def dam_operation(dam_qout, damflg, k, dam_maxfloodq, dam_rate, dam_qin, dam_floodq, dam_vol, ddt, dam_state, dam_vol_temp):
    """
    Operation of the dam

    :param dam_qout:
    :param damflg:
    :param k:
    :param dam_maxfloodq:
    :param dam_rate:
    :param dam_qin:
    :param dam_floodq:
    :param dam_vol:
    :param ddt:
    :param dam_state:
    :param dam_vol_temp:

    :return: dam_qout
    :return: dam_vol_temp
    """
    dam_qout[damflg[k]] = 0.0
    if ( dam_maxfloodq[damflg[k]] == 0 and dam_rate[damflg[k]] == 0 ):
        # constant flood peak cut
        if ( dam_qin[k] < dam_floodq[damflg[k]] ):
            if( dam_vol(damflg[k]) <= 0 ): # 
                dam_qout[damflg[k]] = dam_qin[k]
            else #
                if( dam_qin[k] < 0.25 * dam_floodq[damflg[k]] ):
                    dam_qout[damflg[k]] = 0.25 * dam_floodq[damflg[k]] #
                    qdiff = (dam_qin[k] - dam_qout[damflg[k]]) #
                    dam_vol_temp[damflg[k]] = dam_vol_temp[damflg[k]] + qdiff * ddt #
                else
                    dam_qout[damflg[k]] = dam_qin[k]
                #endif
            #endif
        else
            if (dam_state[damflg[k]] == 0):
                # still have space
                dam_qout[damflg[k]] = dam_floodq[damflg[k]]
                qdiff = (dam_qin[k] - dam_qout[damflg[k]])
                dam_vol_temp[damflg[k]] = dam_vol_temp[damflg[k]] + qdiff * ddt
            else:
                # no more space
                dam_qout[damflg[k]] = dam_qin[k]
            ##endif
        ##endif
    else:
        # non-constant flood peak cut
        if ( dam_qin[k] < dam_floodq[damflg[k]] ):
            if( dam_vol[damflg[k]] <= 0 ): #
                dam_qout[damflg[k]] = dam_qin[k]
            else: #
                if( dam_qin[k] < 0.25 * dam_floodq[damflg[k]] ):
                    dam_qout[damflg[k]] = 0.25 * dam_floodq[damflg[k]] #
                    qdiff = (dam_qin[k] - dam_qout[damflg[k]]) #
                    dam_vol_temp[damflg[k]] = dam_vol_temp[damflg[k]] + qdiff * ddt #
                else:
                    dam_qout[damflg[k]] = dam_qin[k]
                #endif
            #endif
        else:
            if (dam_state[damflg[k]] == 0):
                # still have space
                dam_qout[damflg[k]] = dam_floodq[damflg[k]] + dam_rate[damflg[k]] * (dam_qin[k] - dam_floodq[damflg[k]])
                if( dam_qout[damflg[k]] > dam_maxfloodq[damflg[k]] ):
                    dam_qout[damflg[k]]  = dam_maxfloodq[damflg[k]]
                    qdiff = (dam_qin[k] - dam_qout[damflg[k]])
                    dam_vol_temp[damflg[k]] = dam_vol_temp[damflg[k]] + qdiff * ddt
                else:
                    # no more space
                    dam_qout[damflg[k]] = dam_qin[k]
                ##endif
            ##endif
        ##endif

        return(dam_qout, dam_vol_temp)

#end def dam_operation

#######################################################################

def gate_operation(dam_qout, damflg, k, down_riv_idx, dam_qin, dam_floodq):
    """
    Dam gates operation

    :param dam_qout:
    :param damflg:
    :param k:
    :param down_riv_idx:
    :param dam_qin:
    :param damfloodq:

    :return: dam_qout
    """
    dam_qout[damflg[k]] = 0.0
    kk = down_riv_idx[k]
    dam_qout[damflg[k]] = dam_qin[k]
    if( dam_qout[damflg[k]] >= dam_floodq[damflg[k]] ):
        dam_qout[damflg[k]] = dam_floodq[damflg[k]]
    
    return(dam_qout)
#end def gate_operation

#######################################################################

def dam_checkstate(dam_num, dam_vol_temp, dam_vol, dam_state, dam_vol_max):
    """
    Update the state of the dam

    :param dam_num:
    :param dam_vol_temp:
    :param dam_vol:
    :param dam_state:
    :param dam_vol_max:

    :return: dam_state (0 or 1)
    """
    for i in range( dam_num ):
        dam_vol_temp[i] = dam_vol_temp[i] / 6.0
        dam_vol[i] = dam_vol[i] + dam_vol_temp[i]
        dam_state[i] = 0 # added on January 8, 2021
        if (dam_vol[i] > dam_volmax[i]):
            dam_state[i] = 1
    #enddo

    return(dam_state)
    
#end def dam_checkstate

#######################################################################

def dam_write(f1001, dam_num, dam_vol, dam_qin, dam_loc, dam_qout):
    """
    Write Dam info

    :param f1001: Open file handle
    :param dam_num: Number of dams
    :param dam_vol:
    :param dam_qin:
    :param dam_loc:
    :param dam_qout:

    :return: 0
    """
    for i in range( dam_num ):
        f1001.write(dam_vol[i], dam_qin(dam_loc[i]), dam_qout[i]) # v1.4

    return(0)

#end def dam_write

#######################################################################

def dam_write_cnt(f1002, dam_num, dam_name, dam_iy, dam_ix, dam_volmax, dam_floodq, dam_maxfloodq, dam_rate, dam_vol):
    """
    Write dam detailed info

    :param f1002: Open file 1002 handle
    :param dam_num:
    :param dam_name:
    :param dam_iy:
    :param dam_ix:
    :param dam_volmax:
    :param dam_floodq:
    :param dam_maxfloodq:
    :param dam_rate:
    :param dam_vol:

    :return: 0
    """
    f1002.write( dam_num )
    for i in range( dam_num ):
        f1002.write( dam_name[i], dam_iy[i], dam_ix[i], dam_volmax[i], dam_floodq[i], dam_maxfloodq[i], dam_rate[i], dam_vol[i])
    #enddo

    return(0)

#end def dam_write_cnt
