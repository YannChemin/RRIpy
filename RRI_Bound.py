# RRI_Bound.f90

def read_bound(bound_slo_wlev_switch, bound_slo_wlev ):
    """
    :Description: Read Boundary conditions for slope cells

    Two options ?

    """
    #float bound_opt2[:]
    #int bound_opt2_loci[:], bound_opt2_locj[:]
    #float rdummy
    #float rdummy_dim[:]
    # boundary condition (wlev)
    # for slope cells (wlev)
    if( bound_slo_wlev_switch == 1 or bound_slo_wlev_switch == 2 ):
        f17 = open(boundfile_slo_wlev)
        lines_list = f17.readlines()
        tt = 0
        #if( bound_slo_wlev_switch == 1 ):
        if( bound_slo_wlev_switch == 2 ):
            t, itemp, jtemp = lines_list[0].split(" ")
            if(nx != itemp or ny != jtemp):
                raise Exception ( "error in boundary file (slo, wlev)" )

            for i in range(1, ny):
                rdummy = f17.read()
            #enddo
            tt = tt + 1
        #enddo
        tt_max_bound_slo_wlev = tt - 1
    else # option 1
        num_of_bound_point = f17.read()
        ctemp = f17.read()
        ctemp = f17.read()
        tt = tt + 1
    #enddo
    tt_max_bound_slo_wlev = tt - 1

    #allocate( bound_opt2(num_of_bound_point) )
    bound_opt2 = np.zeros(num_of_bound_point)
    #allocate( bound_opt2_loci(num_of_bound_point), bound_opt2_locj(num_of_bound_point) )
    bound_opt2_loci = np.zeros(num_of_bound_point)
    bound_opt2_locj = np.zeros(num_of_bound_point)

    #endif

    allocate( t_bound_slo_wlev(0:tt_max_bound_slo_wlev), bound_slo_wlev(ny, nx), &
    bound_slo_wlev_idx(0:tt_max_bound_slo_wlev, slo_count), rdummy_dim(slo_count) )
    rewind(17)

 bound_slo_wlev = -999.90
 bound_slo_wlev_idx = -999.90
 rdummy_dim = -999.90

 #if( bound_slo_wlev_switch .eq. 1 ):
 if( bound_slo_wlev_switch .eq. 2 ):

  do tt = 0, tt_max_bound_slo_wlev
   read(17, *) t_bound_slo_wlev(tt), itemp, jtemp
   do i = 1, ny
    read(17,*) (bound_slo_wlev(i, j), j = 1, nx)
   #enddo
   call sub_slo_ij2idx( bound_slo_wlev, rdummy_dim )
   do k = 1, slo_count
    bound_slo_wlev_idx(tt, k) = rdummy_dim[k]
   #enddo
  #enddo
  deallocate( bound_slo_wlev, rdummy_dim )
  print( "done: reading boundary file for slope cells (wlev)"
  close(17)

 else # option 1

  read(17, *) itemp
  read(17, *) ctemp, (bound_opt2_loci[k], k = 1, num_of_bound_point)
  read(17, *) ctemp, (bound_opt2_locj[k], k = 1, num_of_bound_point)

  do tt = 0, tt_max_bound_slo_wlev
   read(17, *) t_bound_slo_wlev(tt), (bound_opt2[k], k = 1, num_of_bound_point)
   do k = 1, num_of_bound_point
    bound_slo_wlev( bound_opt2_loci[k], bound_opt2_locj[k] ) = bound_opt2[k]
   #enddo
   call sub_slo_ij2idx( bound_slo_wlev, rdummy_dim )
   do k = 1, slo_count
    bound_slo_wlev_idx(tt, k) = rdummy_dim[k]
   #enddo
  #enddo
  deallocate( bound_slo_wlev, rdummy_dim )
  print( "done: reading boundary file for slope cells (wlev)"
  close(17)
  deallocate( bound_opt2 )
  deallocate( bound_opt2_loci, bound_opt2_locj )

 #endif
#endif

# for river cells (wlev)

if( bound_riv_wlev_switch .eq. 1 .or. bound_riv_wlev_switch .eq. 2 ):
 open( 18, file = boundfile_riv_wlev, status = 'old' )

 tt = 0
 #if( bound_riv_wlev_switch .eq. 1 ):
 if( bound_riv_wlev_switch .eq. 2 ):

  do
   read(18, *, iostat = ios) t, itemp, jtemp
   if(nx.ne.itemp .or. ny.ne.jtemp) stop "error in boundary file (riv, wlev)"
   do i = 1, ny
    read(18, *, iostat = ios) (rdummy, j = 1, nx)
   #enddo
   if( ios.lt.0 ) exit
   tt = tt + 1
  #enddo
  tt_max_bound_riv_wlev = tt - 1

 else # option 1

  read(18, *) num_of_bound_point
  read(18, '(a)') ctemp
  read(18, '(a)') ctemp
  do
   read(18, *, iostat = ios) t
   if( ios.lt.0 ) exit
   tt = tt + 1
  #enddo
  tt_max_bound_riv_wlev = tt - 1

  allocate( bound_opt2(num_of_bound_point) )
  allocate( bound_opt2_loci(num_of_bound_point), bound_opt2_locj(num_of_bound_point) )

 #endif

 allocate( t_bound_riv_wlev(0:tt_max_bound_riv_wlev), bound_riv_wlev(ny, nx), &
  bound_riv_wlev_idx(0:tt_max_bound_riv_wlev, riv_count), rdummy_dim(riv_count) )
 rewind(18)

 bound_riv_wlev = -999.90
 bound_riv_wlev_idx = -999.90
 rdummy_dim = -999.90

 #if( bound_riv_wlev_switch .eq. 1 ):
 if( bound_riv_wlev_switch .eq. 2 ):

  do tt = 0, tt_max_bound_riv_wlev
   read(18, *) t_bound_riv_wlev(tt), itemp, jtemp
   do i = 1, ny
    read(18,*) (bound_riv_wlev(i, j), j = 1, nx)
   #enddo
   call sub_riv_ij2idx( bound_riv_wlev, rdummy_dim )
   do k = 1, riv_count
    bound_riv_wlev_idx(tt, k) = rdummy_dim[k]
   #enddo
  #enddo
  deallocate( bound_riv_wlev, rdummy_dim )
  print( "done: reading boundary file for river cells (wlev)"
  close(18)

 else # option 1

  read(18, *) itemp
  read(18, *) ctemp, (bound_opt2_loci[k], k = 1, num_of_bound_point)
  read(18, *) ctemp, (bound_opt2_locj[k], k = 1, num_of_bound_point)

  do tt = 0, tt_max_bound_riv_wlev

   read(18, *) t_bound_riv_wlev(tt), (bound_opt2[k], k = 1, num_of_bound_point)
   do k = 1, num_of_bound_point
    bound_riv_wlev( bound_opt2_loci[k], bound_opt2_locj[k] ) = bound_opt2[k]
   #enddo
   call sub_riv_ij2idx( bound_riv_wlev, rdummy_dim )
   do k = 1, riv_count
    bound_riv_wlev_idx(tt, k) = rdummy_dim[k]
   #enddo
  #enddo
  deallocate( bound_riv_wlev, rdummy_dim )
  print( "done: reading boundary file for river cells (wlev)"
  close(18)
  deallocate( bound_opt2 )
  deallocate( bound_opt2_loci, bound_opt2_locj )

 #endif
#endif

# boundary condition (disc)

# for slope cells (disc)

if( bound_slo_disc_switch .eq. 1 .or. bound_slo_disc_switch .eq. 2 ):
 open( 17, file = boundfile_slo_disc, status = 'old' )

 tt = 0
 #if( bound_slo_disc_switch .eq. 1 ):
 if( bound_slo_disc_switch .eq. 2 ):

  do
   read(17, *, iostat = ios) t, itemp, jtemp
   if(nx.ne.itemp .or. ny.ne.jtemp) stop "error in boundary file (slo, disc)"
   do i = 1, ny
    read(17, *, iostat = ios) (rdummy, j = 1, nx)
   #enddo
   if( ios.lt.0 ) exit
   tt = tt + 1
  #enddo
  tt_max_bound_slo_disc = tt - 1

 else # option 1

  read(17, *) num_of_bound_point
  read(17, '(a)') ctemp
  read(17, '(a)') ctemp
  do
   read(17, *, iostat = ios) t
   if( ios.lt.0 ) exit
   tt = tt + 1
  #enddo
  tt_max_bound_slo_disc = tt - 1

  allocate( bound_opt2(num_of_bound_point) )
  allocate( bound_opt2_loci(num_of_bound_point), bound_opt2_locj(num_of_bound_point) )

 #endif

 allocate( t_bound_slo_disc(0:tt_max_bound_slo_disc), bound_slo_disc(ny, nx), &
  bound_slo_disc_idx(0:tt_max_bound_slo_disc, slo_count), rdummy_dim(slo_count) )
 rewind(17)

 bound_slo_disc = -999.90
 bound_slo_disc_idx = -999.90
 rdummy_dim = -999.90

 #if( bound_slo_disc_switch .eq. 1 ):
 if( bound_slo_disc_switch .eq. 2 ):

  do tt = 0, tt_max_bound_slo_disc
   read(17, *) t_bound_slo_disc(tt), itemp, jtemp
   do i = 1, ny
    read(17,*) (bound_slo_disc(i, j), j = 1, nx)
   #enddo
   call sub_slo_ij2idx( bound_slo_disc, rdummy_dim )
   do k = 1, slo_count
    bound_slo_disc_idx(tt, k) = rdummy_dim[k]
   #enddo
  #enddo
  deallocate( bound_slo_disc, rdummy_dim )
  print( "done: reading boundary file for slope cells (disc)"
  close(17)

 else # option 1

  read(17, *) itemp
  read(17, *) ctemp, (bound_opt2_loci[k], k = 1, num_of_bound_point)
  read(17, *) ctemp, (bound_opt2_locj[k], k = 1, num_of_bound_point)

  do tt = 0, tt_max_bound_slo_disc
   read(17, *) t_bound_slo_disc(tt), (bound_opt2[k], k = 1, num_of_bound_point)
   do k = 1, num_of_bound_point
    bound_slo_disc( bound_opt2_loci[k], bound_opt2_locj[k] ) = bound_opt2[k]
   #enddo
   call sub_slo_ij2idx( bound_slo_disc, rdummy_dim )
   do k = 1, slo_count
    bound_slo_disc_idx(tt, k) = rdummy_dim[k]
   #enddo
  #enddo
  deallocate( bound_slo_disc, rdummy_dim )
  print( "done: reading boundary file for slope cells (wlev)"
  close(17)
  deallocate( bound_opt2 )
  deallocate( bound_opt2_loci, bound_opt2_locj )

 #endif
#endif

# for river cells (disc)

if( bound_riv_disc_switch .eq. 1 .or. bound_riv_disc_switch .eq. 2 ):
 open( 18, file = boundfile_riv_disc, status = 'old' )

 tt = 0
 #if( bound_riv_disc_switch .eq. 1 ):
 if( bound_riv_disc_switch .eq. 2 ):

  do
   read(18, *, iostat = ios) t, itemp, jtemp
   if(nx.ne.itemp .or. ny.ne.jtemp) stop "error in boundary file (riv, disc)"
   do i = 1, ny
    read(18, *, iostat = ios) (rdummy, j = 1, nx)
   #enddo
   if( ios.lt.0 ) exit
   tt = tt + 1
  #enddo
  tt_max_bound_riv_disc = tt - 1

 else # option 1

  read(18, *) num_of_bound_point
  read(18, '(a)') ctemp
  read(18, '(a)') ctemp
  do
   read(18, *, iostat = ios) t
   if( ios.lt.0 ) exit
   tt = tt + 1
  #enddo
  tt_max_bound_riv_disc = tt - 1

  allocate( bound_opt2(num_of_bound_point) )
  allocate( bound_opt2_loci(num_of_bound_point), bound_opt2_locj(num_of_bound_point) )

 #endif

 allocate( t_bound_riv_disc(0:tt_max_bound_riv_disc), bound_riv_disc(ny, nx), &
  bound_riv_disc_idx(0:tt_max_bound_riv_disc, riv_count), rdummy_dim(riv_count) )
 rewind(18)

 bound_riv_disc = -999.90
 bound_riv_disc_idx = -999.90
 rdummy_dim = -999.90

 #if( bound_riv_disc_switch .eq. 1 ):
 if( bound_riv_disc_switch .eq. 2 ):

  do tt = 0, tt_max_bound_riv_disc
   read(18, *) t_bound_riv_disc(tt), itemp, jtemp
   do i = 1, ny
    read(18,*) (bound_riv_disc(i, j), j = 1, nx)
   #enddo
   call sub_riv_ij2idx( bound_riv_disc, rdummy_dim )
   do k = 1, riv_count
    bound_riv_disc_idx(tt, k) = rdummy_dim[k]
   #enddo
  #enddo
  deallocate( bound_riv_disc, rdummy_dim )
  print( "done: reading boundary file for river cells (disc)" )
  close(18)

 else # option 1

  read(18, *) itemp
  read(18, *) ctemp, (bound_opt2_loci[k], k = 1, num_of_bound_point)
  read(18, *) ctemp, (bound_opt2_locj[k], k = 1, num_of_bound_point)

  do tt = 0, tt_max_bound_riv_disc
   read(18, *) t_bound_riv_disc(tt), (bound_opt2[k], k = 1, num_of_bound_point)
   do k = 1, num_of_bound_point
    bound_riv_disc( bound_opt2_loci[k], bound_opt2_locj[k] ) = bound_opt2[k]
   #enddo
   call sub_riv_ij2idx( bound_riv_disc, rdummy_dim )
   do k = 1, riv_count
    bound_riv_disc_idx(tt, k) = rdummy_dim[k]
   #enddo
  #enddo
  deallocate( bound_riv_disc, rdummy_dim )
  print( "done: reading boundary file for river cells (disc)"
  close(18)
  deallocate( bound_opt2 )
  deallocate( bound_opt2_loci, bound_opt2_locj )

 #endif
#endif

#end def read_bound
