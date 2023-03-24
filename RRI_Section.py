# RRI_Section

def set_section(sec_id_max):
    """
    Set section

    :param sec_id_max: The maximum value of sec_map

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
            ids = sec_map(i, j)
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

return
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

    else if( h > sec_hr[ids, div_max] ):
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


def hr2vr(hr, k, vr):
#use globals
#implicit none

#real(8) hr, vr, a
#integer k, id, div_max, i

ids = sec_map_idx(k)

if( ids .le. 0 ) then

 vr = hr * area * area_ratio_idx(k)

else

 div_max = sec_div(ids)

 if( hr .le. sec_hr(ids, 1) ) then

  a = hr * sec_b(ids, 1)

 elseif( hr .gt. sec_hr(ids, div_max) ) then

  a = sec_area(ids, div_max) + (hr - sec_hr(ids, div_max)) * sec_b(ids, div_max)

 else

  do i = 2, div_max

   if( hr .le. sec_hr(ids, i) ) then
    a = sec_area(ids, i - 1) + (hr - sec_hr(ids, i - 1)) * sec_b(ids, i)
    exit
   #endif

  #enddo
 #endif

 vr = a * len_riv_idx(k)

#endif

return
#end


def vr2hr( vr, k, hr ):
#use globals
#implicit none

#real(8) hr, vr, a
#integer k, id, div_max, i

ids = sec_map_idx(k)

if( ids .le. 0 ) then

 hr = vr / ( area * area_ratio_idx(k) )

else

 div_max = sec_div(ids)

 a = vr / len_riv_idx(k)

 if( a .le. sec_area(ids, 1) ) then

  hr = a / sec_b(ids, 1)

 elseif( a .gt. sec_area(ids, div_max) ) then

  hr = (a - sec_area(ids, div_max)) / sec_b(ids, div_max) + sec_hr(ids, div_max)

 else

  do i = 2, div_max

   if( a .le. sec_area(ids, i) ) then
    hr = (a - sec_area(ids, i-1)) / sec_b(ids, i) + sec_hr(ids, i - 1)
    exit
   #endif

  #enddo
 #endif

#endif

return
#end


def hr_update(hr_org, vr_inc, k, hr_new):
#use globals
#implicit none

#real(8) hr_org, vr_inc, hr_new, vr_org, vr_new
#integer k

vr_org = hr2vr(hr_org, k)
vr_new = vr_org + vr_inc
hr_new = vr2hr(vr_new, k)

return(hr_new)
#end


def sec_h2b(h, k, b):
#use globals
#implicit none

#real(8) h, b
#integer k, id, div_max, i

ids = sec_map_idx(k)
if( ids .le. 0 ) then

 b = width_idx(k)

else

 div_max = sec_div(ids)

 do i = 1, div_max
  if( h .le. sec_hr(ids, i) ) then
   b = sec_b(ids, i)
   exit
  #endif
  b = sec_b(ids, div_max)
 #enddo

#endif

return
#end
