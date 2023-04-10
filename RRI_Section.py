# RRI_Section

def set_section(sec_id_max, sec_map, secfile, depth, width, height, riv):
    """
    Set section

    :param sec_id_max: The maximum value of sec_map
    :param sec_map:
    :param sec_div:
    :param secfile: section filename
    :param depth:
    :param height:
    :param width:
    :param riv:

    :return: width
    :return: depth
    :return: sec_depth
    :return: height
    :return: sec_height
    :return: sec_hr
    :return: sec_area
    :return: sec_peri
    :return: sec_b
    :return: sec_ns_river
    :return: riv
    """
    #integer i, j, k, id, div_max
    #real(8) rdummy
    #character*256 sec_file_name
    #character*6 id_char

    sec_div = np.zeros(sec_id_max)
    sec_depth = np.zeros(sec_id_max)
    sec_height = np.zeros(sec_id_max)

    div_max = 0
    for i in range(ny):
        for j in range(nx):
            ids = sec_map[i, j]
            if( ids > 0 ):
                id_char = int2char( ids )
                sec_file_name = trim(sec_file) + trim(id_char) + ".txt"
                f10 = open(sec_file_name)
                sec_div[ids], sec_depth[ids], sec_height[ids] = f10.read().split(" ")
                if(div_max <= sec_div[ids]):
                    div_max = sec_div[ids]
                close(10)
            #endif
        #enddo
    #enddo

    sec_hr = np.zeros((sec_id_max, div_max))
    sec_area = np.zeros((sec_id_max, div_max))
    sec_peri = np.zeros((sec_id_max, div_max))
    sec_b = np.zeros((sec_id_max, div_max))
    sec_ns_river = np.zeros((sec_id_max, div_max))

    for i in range(ny):
        for j in range(nx):
            ids = sec_map[i, j]
            if( ids == 0 ):
                id_char = int2char( ids )
                sec_file_name = trim(sec_file) + trim(id_char) + ".txt"
                f10 = open(sec_file_name)
                # TODO this is not sure as no example input file is provided
                sec_div[ids], sec_depth[ids], sec_height[ids] = f10.read().split(" ")
                depth[i, j] = sec_depth[ids]
                height[i, j] = sec_height[ids]
                if(height[i, j] < 0.0):
                    height[i, j] = 0.0 # added ver 1.4.2.4
                for k in range(sec_div[ids]):
                    sec_hr[ids, k], sec_peri[ids, k], sec_b[ids, k], sec_ns_river[ids, k] = f10.read().split(" ")
                #enddo
                sec_area[ids, 1] = sec_b[ids, 1] * sec_hr[ids, 1]
                for k in range(1, sec_div[ids]):
                    sec_area[ids, k] = sec_area[ids, k-1] + sec_b[ids, k] * (sec_hr[ids, k] - sec_hr[ids, k-1])
                #enddo
                #width[i, j] = sec_b[ids, div_max]
                width[i, j] = sec_b[ids, sec_div[ids]] # modified on Feb 14, 2021
                f10.close()
                riv[i, j] = 1 # added by T.Sayama on Dec. 27, 2021
            #endif
        #enddo
    #enddo
    return(width, depth, sec_depth, height, sec_height, sec_div, sec_hr, sec_area, sec_peri, sec_b, sec_ns_river, riv)
#end


def sec_hq_riv(h, dh, k, q):
    """
    Section h q river

    :param sec_map_idx:
    :param k:
    :param sec_div:
    :param sec_hr:
    :param h:
    :param sec_peri:
    :param sec_ns_river:
    :param sec_b:
    :param sec_area:
    
    :return: q the discharge
    """
    #real(8) h, dh, q, p, n, a, b
    ids = sec_map_idx[k]
    div_max = sec_div[ids]

    if( h <= sec_hr[ids, 1] ):
        p = sec_peri[ids, 1]
        n = sec_ns_river[ids, 1]
        a = h * sec_b[ids, 1]

    elif( h > sec_hr[ids, div_max] ):
        p = sec_peri[ids, div_max]
        n = sec_ns_river[ids, div_max]
        a = sec_area[ids, div_max] + (h - sec_hr[ids, div_max]) * sec_b[ids, div_max]
    else:
        for i in range( 1, div_max ):
            if( h <= sec_hr[ids, i] ):
                p = sec_peri[ids, i]
                n = sec_ns_river[ids, i]
                a = sec_area[ids, i - 1] + (h - sec_hr[ids, i - 1]) * sec_b[ids, i]
                break
            #endif
        #enddo
    #endif

    q = 1.0 / n  * ( a / p ) ** (2.0 / 3.0) * sqrt( abs(dh) ) * a

    return(q)
    #end


def hr2vr(sec_map_idx, k, hr, area, area_ratio_idx, sec_div, sec_hr, sec_b, sec_area, len_riv_idx):
    """
    Convert hr to vr
    
    :param sec_map_idx:
    :param k:
    :param hr:
    :param area:
    :param area_ratio_idx:
    :param sec_div:
    :param sec_hr:
    :param sec_b:
    :param sec_area:
    :param len_riv_idx:

    :return: vr
    """
    #real(8) hr, vr, a
    ids = sec_map_idx[k]
    if( ids <= 0 ):
        vr = hr * area * area_ratio_idx[k]
    else:
        div_max = sec_div[ids]
        if( hr <= sec_hr[ids, 1] ):
            a = hr * sec_b[ids, 1]
        elif( hr > sec_hr[ids, div_max] ):
            a = sec_area[ids, div_max] + (hr - sec_hr[ids, div_max]) * sec_b[ids, div_max]
        else:
            for i in range( 1, div_max ):
                if( hr <= sec_hr[ids, i] ):
                    a = sec_area[ids, i - 1] + (hr - sec_hr[ids, i - 1]) * sec_b[ids, i]
                    break
                #endif
            #enddo
        #endif
        vr = a * len_riv_idx[k]
    #endif
    return(vr)
#end


def vr2hr( sec_map_idx, k, vr, area, area_ratio_idx, sec_div, len_riv_idx, sec_area, sec_b, sec_hr ):
    """
    Convert hr to vr

    :param sec_map_idx:
    :param k:
    :param vr:
    :param area:
    :param area_ratio_idx:
    :param sec_div:
    :param len_riv_idx:
    :param sec_area:
    :param sec_b:
    :param sec_hr:
    
    :return: hr
    """ 
    #real(8) hr, vr, a
    ids = sec_map_idx[k]
    if( ids <= 0 ):
        hr = vr / ( area * area_ratio_idx[k] )
    else:
        div_max = sec_div[ids]
        a = vr / len_riv_idx[k]
        if( a <= sec_area[ids, 1] ):
            hr = a / sec_b[ids, 1]
        elif( a > sec_area[ids, div_max] ):
            hr = (a - sec_area[ids, div_max]) / sec_b[ids, div_max] + sec_hr[ids, div_max]
        else:
            for i in range(1, div_max):
                if( a <= sec_area[ids, i] ):
                        hr = (a - sec_area[ids, i-1]) / sec_b[ids, i] + sec_hr[ids, i - 1]
                        break
                #endif
            #enddo
        #endif
    #endif
    return(hr)
#end


def hr_update(hr_org, vr_inc, k):
    """
    Update hr values

    :param hr_org:
    :param k:
    :param vr_inc:

    :return: hr
    """
    vr_org = hr2vr(hr_org, k)
    vr_new = vr_org + vr_inc
    hr_new = vr2hr(vr_new, k)

    return(hr_new)
#end


def sec_h2b(sec_map_idx, k, width_idx, sec_div, sec_hr, sec_b):
    """
    Section h2b
    
    :param sec_map_idx:
    :param k:
    :param width_idx:
    :param sec_div:
    :param sec_hr:
    :param sec_b:

    :return: b
    """
    ids = sec_map_idx[k]
    if( ids <= 0 ):
        b = width_idx[k]
    else:
        div_max = sec_div[ids]
        for i in range( div_max ):
            if( h <= sec_hr[ids, i] ):
                b = sec_b[ids, i]
                break
            #endif
            b = sec_b[ids, div_max]
        #enddo
    #endif
    return(b)
#end
