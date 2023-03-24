from globals_vars import *
from dam_mod import *
import numpy

def RRI_TSAS(tt, hs_idx, hr_idx, hg_idx, qs_ave_idx, qr_ave_idx, qg_ave_idx, qp_t_idx):
    """
    :Description: Output for T-SAS program by RRI
    
    :param tt:
    :param hs_idx:
    :param hr_idx:
    :param hg_idx:
    :param qs_ave_idx:
    :param qr_ave_idx:
    :param qe_ave_idx:
    :param qp_t_idx:
    :result: output files for T-SAS program 
    """
    tsasfolder = "./tsas/"
    tsasdata = "./tsas/data/"

    t_char=int2char(tt)

    # Initial Output
    if(tt == 0):
        ofile_tsas_hs_id = trim(tsasfolder) + "hs_id.txt"
        ofile_tsas_hr_id = trim(tsasfolder) + "hr_id.txt"
        ofile_tsas_outlet = trim(tsasfolder) + "outlet.txt"
        ofile_tsas_hg_param = trim(tsasfolder) + "param.txt"

        f1101 = open(ofile_tsas_hs_id)
        f1102 = open(ofile_tsas_hr_id)
        f1103 = open(ofile_tsas_outlet)
        f1104 = open(ofile_tsas_hg_param)

        #allocate( hs_id(ny, nx), hr_id(ny, nx) )
        hs_id = np.zeros((ny, nx))
        hr_id = np.zeros((ny, nx))
        #allocate( hs_id_idx(slo_count), hr_id_idx(riv_count) )
        hs_id_idx = np.zeros(slo_count)
        hr_id_idx = np.zeros(riv_count)
        
        hs_id = np.fill(-999)
        hr_id = np.fill(-999)

        tsasid = 0

        for k in range(1, slo_count):
            tsasid = tsasid + 1
            i = slo_idx2i[k]
            j = slo_idx2j[k]
            hs_id[i, j] = tsasid
            hs_id_idx[k] = tsasid
            if( domain_slo_idx[k] == 2 ):
                f1103.write(tsasid)
        #enddo

        for k in range(1, riv_count):
            tsasid = tsasid + 1
            i = riv_idx2i[k]
            j = riv_idx2j[k]
            hr_id[i, j] = tsasid
            hr_id_idx[k] = tsasid
            if( domain_riv_idx[k] == 2 ):
                f1103.write(tsasid)
        #enddo

        f1101.write(("ncols %d" % (nx))
        f1101.write("nrows %d" % (ny))
        f1101.write("xllcorner %22.12f" % (xllcorner))
        f1101.write("yllcorner %22.12f" % (yllcorner))
        f1101.write("cellsize %20.12f" % (cellsize))
        f1101.write("NODATA_value %d" % (-999))

        f1102.write("ncols %d" % (nx))
        f1102.write("nrows %d" % (ny))
        f1102.write("xllcorner  %22.12f" % (xllcorner))
        f1102.write("yllcorner %22.12f" % (yllcorner))
        f1102.write("cellsize %20.12f" % (cellsize))
        f1102.write("NODATA_value %d" % (-999))

        for i in range(1, ny):
            # These may have to be cleaned from "[" and "]"
            f1101.write(hs_id[i])
            f1102.write(hr_id[i])
        #enddo

        f1104.write(ns_slope[1])
        f1104.write(soildepth[1])
        f1104.write(gammaa[1])
        f1104.write(ka[1])
        f1104.write(gammam[1])
        f1104.write(beta[1])
        f1104.write(ksg[1])
        f1104.write(gammag[1])
        f1104.write(kg0[1])
        f1104.write(fpg[1])
        f1104.write(rgl[1])
        f1104.write(area)

        f1101.close()
        f1102.close()
        f1103.close()
        f1104.close()
    #endif

    # Data Output
    t_char = int2char(tt)

    ofile_tsas_sto = trim(tsasdata) + "tsas_sto_" + trim(t_char) + ".txt"
    ofile_tsas_flux = trim(tsasdata) + "tsas_flux_" + trim(t_char) + ".txt"
    ofile_tsas_rain = trim(tsasdata) + "tsas_rain_" + trim(t_char) + ".txt"

    f1103 = open(ofile_tsas_sto)
    f1104 = open(ofile_tsas_flux)
    f1105 = open(ofile_tsas_rain)

    # storage
    for k in range(1, slo_count):
        f1103.write(hs_idx[k], hg_idx[k], 1)
    #enddo
    for k in range(1, riv_count):
        vr_temp = hr2vr(hr_idx[k], k)
        f1103.write(vr_temp, 0.0, 2)
    #enddo

    # flux [m3/s]

    # qs
    for k in range(1, slo_count):
        for l in range(1, lmax):
            dif_p = dif_slo_idx[k]
            if( dif_p == 0 and l == 2 ):
                exit # kinematic -> 1-direction
            if( dif_p == 0 ):
                kk = down_slo_1d_idx[k]
            kk = down_slo_idx[l, k]
            if( kk == -1 ):
                continue

            # calc dh/dx
            distance = dis_slo_idx[l, k]
            if( dif_p == 0 ):
                distance = dis_slo_1d_idx[k]
            zb_p = zb_slo_idx[k]
            hs_p = hs_idx[k]
            zb_n = zb_slo_idx[kk]
            hs_n = hs_idx[kk]
            call h2lev(hs_p, k, lev_p)
            call h2lev(hs_n, kk, lev_n)
            dh = ((zb_p + lev_p) - (zb_n + lev_n)) / distance 
            if( dif_p == 0 ):
                dh = max( (zb_p - zb_n) / distance, 0.001 )

            f1104.write(hs_id_idx[k], hs_id_idx[kk],(qs_ave_idx[l, k] * area), dh, (qg_ave_idx[l, k] * area), 1)
        #enddo
    #enddo

    # qr
    for k in range(1, riv_count):
        if( domain_riv_idx[k] == 2) continue
        kk = down_riv_idx[k]
        f1104.write(hr_id_idx[k], hr_id_idx[kk], qr_ave_idx[k], 0.0, 0.0, 2) # v1.4
    #enddo

    # qrs
    for i in range (1, ny):
        for j in range (1, nx):
            if( domain[i, j] == 0 ):
                continue
            if( riv[i, j] == 0 ):
                continue
            f1104.write(hs_id[i, j], hr_id[i, j], (qrs[i, j] * area), 0.0, 0.0, 3)
        #enddo
    #enddo

    # id, rainfall, aevp, infilt
    for k in range(1, slo_count):
        i = slo_idx2i[k]
        j = slo_idx2j[k]
        f1105.write(hs_id_idx[k], (qp_t_idx[k] * area), (aevp_tsas[k] * area), (gampt_f[i, j] * area), (rech_hs_tsas[k] * area), (exfilt_hs_tsas[k] * area))
    #enddo

    f1103.close()
    f1104.close()
    f1105.close()

#end def RRI_TSAS
