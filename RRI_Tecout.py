#*********************************************************************
def Tecout_alloc(nx,ny,ValNum_temp):
    """
    Old Memory allocation function (NOT USED)
    """
    iMX = ny + 1
    jMX = nx + 1
    ValNum = ValNum_temp
    X = np.zeros(jMX)
    Y = np.zeros(iMX)
    Z = np.zeros((iMX,jMX))
    Z_buf = np.zeros((ny,nx))
    SufHmax = np.zeros((iMX-1,jMX-1))
    SufHmax.fill(-0.10)
    f6 = open( tecfile )
#enddef Tecout_alloc

#*********************************************************************

def Tecout_mkGrid(DX, DY, Z_temp, X, Y, Y_temp, Z_buf, Z_temp):
    """
    """
    X[0] = 0.
    for j in range(1, jMX):
        X[j] = X[j-1] + DX
    #enddo
    #Yç¿ïWçÏê¨
    Y_temp[0] = 0.
    for i in range(1, iMX):
        Y_temp[i] = Y_temp[i-1] + DY
    #enddo
    for i in range(iMX, 0, -1):
        Y[iMX-i] = Y_temp[i]
    #enddo
    #Zç¿ïWçÏê¨-Cnt2Node
    for i in range( iMX ):
        for j in range(jMX-1):
            Z_buf[i,j] = Z_temp[i,j]
            if (Z_buf(i,j) < 0.):
                Z_buf[i,j] = 0.
            #enddo
        #enddo
    for i in range( iMX ):
        for j in range( jMX ):
            if (i==0):
                if (j==0):
                    if (Z_buf[i,j] > 0.):
                        Z[i,j] = Z_buf[i,j]
                    else:
                        Z[i,j] = 0.
                    #endif
                else if (j == jMX):
                    if (Z_buf[i,j] > 0.):
                          Z[i,j] = Z_buf[i,j]
                      else:
                          Z[i,j] = 0.
                      #endif
                else:
                    divN = 0.
                    if (Z_buf(i,j-1) > 0.):
                        divN = divN + 1.
                    if (Z_buf[i,j] > 0.):
                        divN = divN + 1.
                    if (divN > 0.):
                        Z[i,j] = (Z_buf[i,j-1] + Z_buf[i,j])/divN
                    else
                        Z(i,j) = 0.
                    #endif
                #endif
            else if (i == iMX):
                if (j == 1):
                      if (Z_buf[i-1,j] > 0.):
                          Z[i,j] = Z_buf[i-1,j]
                      else:
                          Z[i,j] = 0.
                      #endif
                else if (j==jMX):#TODO
                      if (Z_buf(i-1,j-1)>0.) then
                          Z(i,j)=Z_buf(i-1,j-1)
                      else
                          Z(i,j)=0.
                      #endif
                  else
                      divN=0.
                      if (Z_buf(i-1,j-1)>0.) divN=divN+1.
                      if (Z_buf(i-1,j  )>0.) divN=divN+1.
                      if (divN>0.) then
                          Z(i,j)=(Z_buf(i-1,j-1)+Z_buf(i-1,j))/divN
                      else
                          Z(i,j)=0.
                      #endif
                  #endif
              else
                  if (j==1) then
                      divN=0.
                      if (Z_buf(i-1,j)>0.) divN=divN+1.
                      if (Z_buf(i  ,j)>0.) divN=divN+1.
                      if (divN>0.) then
                          Z(i,j)=(Z_buf(i-1,j)+Z_buf(i,j))/divN
                      else
                          Z(i,j)=0.
                      #endif
                  elseif (j==jMX) then
                      divN=0.
                      if (Z_buf(i-1,j-1)>0.) divN=divN+1.
                      if (Z_buf(i  ,j-1)>0.) divN=divN+1.
                      if (divN>0.) then
                          Z(i,j)=(Z_buf(i-1,j-1)+Z_buf(i,j-1))/divN
                      else
                          Z(i,j)=0.
                      #endif
                  else
                      divN=0.
                      if (Z_buf(i-1,j-1)>0.) divN=divN+1.
                      if (Z_buf(i  ,j-1)>0.) divN=divN+1.
                      if (Z_buf(i-1,j  )>0.) divN=divN+1.
                      if (Z_buf(i  ,j  )>0.) divN=divN+1.
                      if (divN>0.) then
                          Z(i,j)= &
                           ( Z_buf(i-1,j-1) + Z_buf(i  ,j  ) + &
                             Z_buf(i-1,j  ) + Z_buf(i  ,j-1) ) / divN
                      else
                          Z(i,j)=0.
                      #endif
                  #endif
              #endif
          #enddo
      #enddo
#enddef  Tecout_mkGrid
#*********************************************************************
def Tecout_write_initialize( k, val1, val2, val3, val4)
    use Tecout_Mod
    implicit none
#---define
    integer :: i,j,k
    real(8) :: val1(iMX-1,jMX-1), val2(iMX-1,jMX-1), val3(iMX-1,jMX-1), val4(iMX-1,jMX-1)

    #ç≈èâÇ…ã§í ïîï™ÇèoóÕ
        write(6,'(a)') &
          ' VARIABLES = "X","Y","Z","width","depth","height","area_ratio" &
                        "Rain","River_H","River_Q","Surface_H","Surface_H_Max" '
        write(6,'(a,i5,a)') 'ZONE T="',k,'"'
        write(6,'(2(a,i5),a)') &
          'I=' ,jMX, ',J=' ,iMX, ',DATAPACKING=BLOCK'
        write(6,'(a,i10)') 'SOLUTIONTIME=',k
        write(6,'(a)') 'VARLOCATION=([4-12]=CELLCENTERED)'
        
         for i in range( iMX
             for j=1,jMX
                 write(6,'(f12.3,",",$)') X(j)
             #enddo
             write(6,*)
         #enddo
         for i in range( iMX
             for j=1,jMX
                 write(6,'(f12.3,",",$)') Y(i)
             #enddo
             write(6,*)
         #enddo
         for i in range( iMX
             for j=1,jMX
                 write(6,'(f12.3,",",$)') Z(i,j)
             #enddo
             write(6,*)
         #enddo
         for i in range( iMX-1
             for j=1,jMX-1
                 write(6,'(f12.5,",",$)') val1(i,j)
             #enddo
             write(6,*)
         #enddo
         for i in range( iMX-1
             for j=1,jMX-1
                 write(6,'(f12.5,",",$)') val2(i,j)
             #enddo
             write(6,*)
         #enddo
         for i in range( iMX-1
             for j=1,jMX-1
                 write(6,'(f12.5,",",$)') val3(i,j)
             #enddo
             write(6,*)
         #enddo
         for i in range( iMX-1
             for j=1,jMX-1
                 write(6,'(f12.5,",",$)') val4(i,j)
             #enddo
             write(6,*)
         #endfor    
#enddef Tecout_write_initialize    
#*********************************************************************
def Tecout_write( k, val1, val2, val3, val4, area )
    use Tecout_Mod
    implicit none
#---define
    integer :: i,j,k
    real(8) :: val1(iMX-1,jMX-1), val2(iMX-1,jMX-1), val3(iMX-1,jMX-1), val4(iMX-1,jMX-1), area
                #"Rain","River_H","River_Q","Surface_H"
    if (k==1) then
    #ÇPâÒñ⁄
#        write(6,'(a)') &
#          'VARIABLES = "X","Y","Z","Rain","River_H","River_Q",Surface'
#        write(6,'(a,i5,a)') 'ZONE T="',k,'"'
#        write(6,'(2(a,i5),a)') &
#          'I=' ,jMX, ',J=' ,iMX, ',DATAPACKING=BLOCK'
#        write(6,'(a,i10)') 'SOLUTIONTIME=',k
#        write(6,'(a)') 'VARLOCATION=([4-7]=CELLCENTERED)'
#        
#         for i in range( iMX
#         #for i=iMX,1,-1
#             for j=1,jMX
#                 write(6,'(f12.3,",",$)') X(j)
#             #enddo
#             write(6,*)
#         #enddo
#         for i in range( iMX
#         #for i=iMX,1,-1
#             for j=1,jMX
#                 write(6,'(f12.3,",",$)') Y(i)
#             #enddo
#             write(6,*)
#         #enddo
#         for i in range( iMX
#         #for i=iMX,1,-1
#             for j=1,jMX
#                 write(6,'(f12.3,",",$)') Z(i,j)
#             #enddo
#             write(6,*)
#         #enddo
         for i in range( iMX-1
             for j=1,jMX-1
                 write(6,'(f12.5,",",$)') val1(i,j) * 3600.0 * 1000.0
             #enddo
             write(6,*)
         #enddo
         for i in range( iMX-1
             for j=1,jMX-1
                 write(6,'(f12.5,",",$)') val2(i,j)
             #enddo
             write(6,*)
         #enddo
         for i in range( iMX-1
             for j=1,jMX-1
                 write(6,'(f12.5,",",$)') val3(i,j) * area
             #enddo
             write(6,*)
         #enddo
         for i in range( iMX-1
             for j=1,jMX-1
                 write(6,'(f12.5,",",$)') val4(i,j)
             #enddo
             write(6,*)
         #enddo
         where(SufHmax < val4 ) SufHmax = val4
         for i in range( iMX-1
             for j=1,jMX-1
                 write(6,'(f12.5,",",$)') SufHmax(i,j)
             #enddo
             write(6,*)
         #enddo
         
      else
      #2âÒñ⁄à»ç~
        #write(6,'(a)') &
        #  'VARIABLES = "X","Y","Z","Rain","River_H","River_Q",Surface'
        write(6,'(a,i5,a)') 'ZONE T="',k,'"'
        write(6,'(2(a,i5),a)') &
          'I=' ,jMX, ',J=' ,iMX, ',DATAPACKING=BLOCK'
        write(6,'(a,i10)') 'SOLUTIONTIME=',k
        write(6,'(a)') 'VARLOCATION=([4-12]=CELLCENTERED)'
        write(6,'(a)') 'VARSHARELIST=([1-7]=1)'
        
         for i in range( iMX-1
             for j=1,jMX-1
                 write(6,'(f12.5,",",$)') val1(i,j) * 3600.0 * 1000.0
             #enddo
             write(6,*)
         #enddo
         for i in range( iMX-1
             for j=1,jMX-1
                 write(6,'(f12.5,",",$)') val2(i,j)
             #enddo
             write(6,*)
         #enddo
         for i in range( iMX-1
             for j=1,jMX-1
                 write(6,'(f12.5,",",$)') val3(i,j) * area
             #enddo
             write(6,*)
         #enddo
         for i in range( iMX-1
             for j=1,jMX-1
                 write(6,'(f12.5,",",$)') val4(i,j)
             #enddo
             write(6,*)
         #enddo
         where(SufHmax < val4 ) SufHmax = val4
         for i in range( iMX-1
             for j=1,jMX-1
                 write(6,'(f12.5,",",$)') SufHmax(i,j)
             #enddo
             write(6,*)
         #enddo
         
      #endif
#enddef Tecout_write


