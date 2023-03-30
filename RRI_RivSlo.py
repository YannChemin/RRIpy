# RRI_RivSlo

def funcrs( hr, hs )
    """
    River and slope interaction
    """
    #real(8) hr(ny, nx), hs(ny, nx)
    #real(8) hrs # discharge amount from slope to river [m/s]
    #real(8) hr_top, hs_top, h1, h2
    #real(8) mu1, mu2, mu3
    #real(8) hr_new, len, b, ar # add v1.4
    mu1 = (2.0 / 3.0) ** (3.0 / 2.0)
    mu2 = 0.350
    mu3 = 0.910
    for i in range( ny ):
        for j in range( nx ):
            if( domain[i, j] == 0 ):
                continue
            if( riv[i, j] == 0 ):
                continue
            hs_top = hs[i, j]
            hr_top = hr[i, j] - depth[i, j]
            k = riv_ij2idx[i, j]
            length = len_riv_idx[k]
            # (Case a) : (height = 0 and hr_top < 0) or (height > 0 and hr_top < 0 and hs_top <= height)
            # -> From slope to river : step fall (hrs : positive)
            if(( height[i, j] == 0.0 and hr_top < 0.0 ) or ( height[i, j] > 0.0 and hr_top < 0.0 and hs_top <= height[i, j] )):
                #hrs = mu1 * hs_top * sqrt( 9.810 * hs_top ) * dt * length / area
                hrs = mu1 * hs_top * sqrt( 9.810 * hs_top ) * dt * length * 2.0 / area # modified 1.4.2.4 
                if( hrs > hs[i, j] ):
                    hrs = hs[i, j]
                hs[i, j] = hs[i, j] - hrs
                #hr[i, j] = hr[i, j] + hrs / area_ratio[i, j]
                call hr_update(hr[i, j], hrs * area, k, hr_new) # add v1.4
                hr[i, j] = hr_new # add v1.4
                qrs[i, j] = hrs
                # avoid the situation of hr_top > hs_top
                hs_top = hs[i, j]
                hr_top = hr[i, j] - depth[i, j]
                if( hr_top >= -0.000010 and hr_top > hs_top ):
                    for count in range( 10 ):
                        call sec_h2b(hr[i, j], k, b)
                        ar = len * b / area
                        hrs = ( hs_top - hr_top ) / ( 1.0 + 1.0 / ar )
                        hs[i, j] = hs[i, j] - hrs
                        #hr[i, j] = hr[i, j] + hrs / area_ratio[i, j]
                        call hr_update(hr[i, j], hrs * area, k, hr_new) # add v1.4
                        hr[i, j] = hr_new # add v1.4
                        qrs[i, j] = qrs[i, j] + hrs
                        if( abs(hs[i, j] - (hr[i, j] - depth[i, j])) < 0.000010 ):
                            break
                        hs_top = hs[i, j]
                        hr_top = hr[i, j] - depth[i, j]
                    #enddo
                    hr[i, j] = hs[i, j] + depth[i, j] # 
                #endif
            # (Case b) : height > 0 and hr_top <= height and hr_top >= 0
            # -> No exchange
            else if( height[i, j] > 0.0 and hs_top <= height[i, j] and hr_top <= height[i, j] and hr_top >= 0.0 ):
                qrs[i, j] = 0.0
                continue
            # (Case c) : hs <= hrt & hrt >= height
            # -> From river to slope : overtopping (hrs : negative)
            # (incl. hs = 0 and hrt > 0)
            else if( hs_top <= hr_top and hr_top >= height[i, j] ):
                h1 = hr_top - height[i, j]
                h2 = hs_top - height[i, j]
                if( h2 / h1  <=  2.0 / 3.0 ):
                    #hrs = - mu2 * h1 * sqrt( 2.0 * 9.810 * h1 ) * dt * length / area # modified v1.4
                    hrs = - mu2 * h1 * sqrt( 2.0 * 9.810 * h1 ) * dt * length * 2.0 / area # modified v1.4.2.4
                else:
                    #hrs = - mu3 * h2 * sqrt( 2.0 * 9.810 * (h1 - h2) ) * dt * length / area # modified v1.4
                    hrs = - mu3 * h2 * sqrt( 2.0 * 9.810 * (h1 - h2) ) * dt * length * 2.0 / area # modified v1.4.2.4
                #endif
                call sec_h2b(hr[i, j], k, b)
                ar = length * b / area
                if( abs(hrs / ar) > (hr_top-height[i, j]) ):
                    hrs = - (hr_top-height[i, j]) * ar
                    #if( abs(hrs / area_ratio[i, j]) .gt. (hr_top-height[i, j]) ) hrs = - (hr_top-height[i, j]) * area_ratio[i, j]
                    qrs[i, j] = hrs
                    hs[i, j] = hs[i, j] - hrs
                    #hr[i, j] = hr[i, j] + hrs / area_ratio[i, j]
                    call hr_update(hr[i, j], hrs * area, k, hr_new) # add v1.4
                    hr[i, j] = hr_new # add v1.4
                    # avoid the situation of hs_top > hr_top
                    hs_top = hs[i, j]
                    hr_top = hr[i, j] - depth[i, j]
                    #if( hr_top > -0.000010 and hs_top > hr_top ): # modified from v1.4 (‚±‚±‚Í‚Ç‚Á‚¿‚©H)
                    if( hs_top > hr_top ):
                        for count in range( 10 ):
                            call sec_h2b(hr[i, j], k, b)
                            ar = len * b / area
                            hrs = ( hs_top - hr_top ) / ( 1.0 + 1.0 / ar )
                            hs[i, j] = hs[i, j] - hrs
                            #hr[i, j] = hr[i, j] + hrs / area_ratio[i, j]
                            call hr_update(hr[i, j], hrs * area, k, hr_new) # add v1.4
                            hr[i, j] = hr_new # add v1.4
                            qrs[i, j] = qrs[i, j] + hrs
                            if( abs(hs[i, j] - (hr[i, j] - depth[i, j])) .lt. 0.000010 ):
                                break
                            hs_top = hs[i, j]
                            hr_top = hr[i, j] - depth[i, j]
                        #enddo
                        hr[i, j] = hs[i, j] + depth[i, j] # ÅIŽè’i
                    #endif
            # (Case d) : hs > hrt & hs >= height
            # -> From slope to river : overtopping (hrs : positive)
            # (incl. hs = 0 and hrt > 0)
            else if( hs_top >= hr_top and hs_top >= height[i, j] ):
                h1 = hs_top - height[i, j]
                h2 = hr_top - height[i, j]
                if( h2 / h1  <=  2.0 / 3.0 ):
                    #hrs = mu2 * h1 * sqrt( 2.0 * 9.810 * h1 ) * dt * length / area # v1.4
                    hrs = mu2 * h1 * sqrt( 2.0 * 9.810 * h1 ) * dt * length * 2.0 / area # v1.4.2.4
                else:
                    #hrs = mu3 * h2 * sqrt( 2.0 * 9.810 * (h1 - h2) ) * dt * length / area # v1.4
                    hrs = mu3 * h2 * sqrt( 2.0 * 9.810 * (h1 - h2) ) * dt * length * 2.0 / area # v1.4.2.4
                #endif
                if( hrs > (hs_top - height[i, j]) ):
                    hrs = hs[i, j] - height[i, j]
                    qrs[i, j] = hrs
                    hs[i, j] = hs[i, j] - hrs
                    #hr[i, j] = hr[i, j] + hrs / area_ratio[i, j]
                    call hr_update(hr[i, j], hrs * area, k, hr_new) # add v1.4
                    hr[i, j] = hr_new # add v1.4
                    # avoid the situation of hr_top > hs_top
                    hs_top = hs[i, j]
                    hr_top = hr[i, j] - depth[i, j]
                    if( hr_top >= -0.000010 and hr_top > hs_top ):
                        for count in range( 10 ):
                            call sec_h2b(hr[i, j], k, b)
                            ar = length * b / area
                            hrs = ( hs_top - hr_top ) / ( 1.0 + 1.0 / ar )
                            hs[i, j] = hs[i, j] - hrs
                            #hr[i, j] = hr[i, j] + hrs / area_ratio[i, j]
                            call hr_update(hr[i, j], hrs * area, k, hr_new) # add v1.4
                            hr[i, j] = hr_new # add v1.4
                            qrs[i, j] = qrs[i, j] + hrs
                            if( abs(hs[i, j] - (hr[i, j] - depth[i, j])) < 0.000010 ):
                                break
                            hs_top = hs[i, j]
                            hr_top = hr[i, j] - depth[i, j]
                        #enddo
                        hr[i, j] = hs[i, j] + depth[i, j] #
                    #endif
            else:
                # Condition not considered above
                raise Exception ( "Error : RivSlo" )
            #endif
            qrs[i, j] = qrs[i, j] / float(dt) # [m/s]
        #enddo
    #enddo
#end def funcrs
