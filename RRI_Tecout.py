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
    Initialize a grid

    :param DX:
    :param DY:
    :param Z_temp:
    :param X:
    :param Y:
    :param Y_temp:
    :param Z_buf:
    :param Z_temp:

    :return: Z
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
                else if (j == jMX):#TODO
                      if (Z_buf[i-1,j-1) > 0.):
                          Z[i,j] = Z_buf[i-1,j-1]
                      else:
                          Z[i,j] = 0.
                      #endif
                  else:
                      divN = 0.
                      if (Z_buf[i-1,j-1] > 0.):
                        divN = divN + 1.
                      if (Z_buf[i-1,j) > 0.):
                        divN = divN + 1.
                      if (divN > 0.):
                          Z[i,j] = (Z_buf[i-1,j-1] + Z_buf[i-1,j])/divN
                      else:
                          Z(i,j) = 0.
                      #endif
                  #endif
              else:
                  if (j == 1):
                      divN = 0.
                      if (Z_buf[i-1,j] > 0.):
                        divN = divN + 1.
                      if (Z_buf[i,j] > 0.):
                        divN = divN + 1.
                      if (divN > 0.):
                          Z[i,j] = (Z_buf[i-1,j]+Z_buf[i,j])/divN
                      else:
                          Z[i,j] = 0.
                      #endif
                  else if (j == jMX):
                      divN = 0.
                      if (Z_buf[i-1,j-1] > 0.):
                        divN = divN + 1.
                      if (Z_buf[i,j-1] > 0.):
                        divN = divN + 1.
                      if (divN > 0.):
                          Z[i,j] = (Z_buf[i-1,j-1]+Z_buf[i,j-1])/divN
                      else:
                          Z[i,j] = 0.
                      #endif
                  else:
                      divN = 0.
                      if (Z_buf[i-1,j-1] > 0.):
                        divN = divN + 1.
                      if (Z_buf[i,j-1] > 0.):
                        divN = divN + 1.
                      if (Z_buf[i-1,j] > 0.):
                        divN = divN + 1.
                      if (Z_buf[i,j] > 0.):
                        divN = divN + 1.
                      if (divN > 0.):
                          Z[i,j] = (Z_buf[i-1,j-1] + Z_buf[i,j] + Z_buf[i-1,j] + Z_buf[i,j-1]) / divN
                      else:
                          Z[i,j] = 0.
                      #endif
                  #endif
              #endif
          #enddo
      #enddo
      return(Z)
#enddef  Tecout_mkGrid

#*********************************************************************

def Tecout_write_initialize( f6, k, iMX, jMX, X, Y, Z, val1, val2, val3, val4):
    """
    Start writing Tecout

    :param f6:
    :param k:
    :param iMX:
    :param jMX:
    :param X:
    :param Y:
    :param Z:
    :param val1:
    :param val2:
    :param val3:
    :param val4:

    :return: Filling the f6 file with data
    """
    #real(8) :: val1(iMX-1,jMX-1), val2(iMX-1,jMX-1), val3(iMX-1,jMX-1), val4(iMX-1,jMX-1)

    #ç≈èâÇ…ã§í ïîï™ÇèoóÕ
    f6.write(' VARIABLES = "X","Y","Z","width","depth","height","area_ratio","Rain","River_H","River_Q","Surface_H","Surface_H_Max" ')
    f6.write( 'ZONE T="',k,'"')
    f6.write( 'I=' ,jMX, ',J=' ,iMX, ',DATAPACKING=BLOCK')
    f6.write( 'SOLUTIONTIME=',k)
    f6.write( 'VARLOCATION=([4-12]=CELLCENTERED)')
    for i in range( iMX ):
        for j in range( jMX ):
            f6.write( X[j] )
        #enddo
    #enddo
    for i in range( iMX ):
        for j in range( jMX ):
            f6.write( Y[i] )
        #enddo
    #enddo
    for i in range( iMX ):
        for j in range( jMX ):
            f6.write( Z[i,j] )
        #enddo
    #enddo
    for i in range( iMX-1 ):
        for j in range( jMX-1 ):
            f6.write( val1[i,j] )
        #enddo
    #enddo
    for i in range( iMX-1 ):
        for j in range( jMX-1 ):
            f6.write( val2[i,j] )
        #enddo
    #enddo
    for i in range( iMX-1 ):
        for j in range( jMX-1 ):
            f6.write( val3[i,j] )
        #enddo
    #enddo
    for i in range( iMX-1 ):
        for j in range( jMX-1 ):
            f6.write( val4[i,j] )
        #enddo
    #endfor    
    return(0)#to return something
#enddef Tecout_write_initialize   

#*********************************************************************

def Tecout_write( k, f6, iMX, jMX, val1, val2, val3, val4, area, SufHmax )
    """
    Write data in Tecout

    :param k:
    :param f6:
    :param iMX:
    :param jMX:
    :param val1:
    :param val2:
    :param val3:
    :param val4:
    :param area:
    :param SufHmax:

    :return: Write into f6 Tecout file
    """
    #real(8) :: val1(iMX-1,jMX-1), val2(iMX-1,jMX-1), val3(iMX-1,jMX-1), val4(iMX-1,jMX-1), area
                #"Rain","River_H","River_Q","Surface_H"
    if (k == 1):
        #ÇPâÒñ⁄
        for i in range( iMX-1 ):
            for j in range( jMX-1 ):
                f6.write( val1[i,j] * 3600.0 * 1000.0 )
            #enddo
        #enddo
        for i in range( iMX-1 ):
            for j in range( jMX-1 ):
                f6.write( val2[i,j] )
            #enddo
        #enddo
        for i in range( iMX-1 ):
            for j in range( jMX-1 ):
                f6.write( val3[i,j] * area )
            #enddo
        #enddo
        for i in range( iMX-1 ):
            for j in range( jMX-1 ):
                f6.write( val4[i,j] )
            #enddo
        #enddo
        SufHmax = np.where(SufHmax < val4, val4, SufHmax)
        for i in range( iMX-1 ):
            for j in range( jMX-1 ):
                f6.write( SufHmax[i,j] )
            #enddo
        #enddo
    else:
        #2âÒñ⁄à»ç~
        #write(6,'(a)')  'VARIABLES = "X","Y","Z","Rain","River_H","River_Q",Surface'
        f6.write( 'ZONE T="',k,'"' )
        f6.write( 'I=' ,jMX, ',J=' ,iMX, ',DATAPACKING=BLOCK' )
        f6.write( 'SOLUTIONTIME=',k )
        f6.write( 'VARLOCATION=([4-12]=CELLCENTERED)' )
        f6.write( 'VARSHARELIST=([1-7]=1)' )
        for i in range( iMX-1 ):
             for j in range( jMX-1 ):
                f6.write( val1[i,j] * 3600.0 * 1000.0 )
             #enddo
        #enddo
        for i in range( iMX-1 ):
            for j in range( jMX-1 ):
                f6.write( val2[i,j] )
            #enddo
        #enddo
        for i in range( iMX-1 ):
            for j in  range( jMX-1 ):
                f6.write( val3[i,j] * area )
            #enddo
        #enddo
        for i in range( iMX-1 ):
            for j=1,jMX-1 ):
                f6.write( val4[i,j] )
            #enddo
        #enddo
        SufHmax = np.where(SufHmax < val4, val4, SufHmax )
        for i in range( iMX-1 ):
            for in range( jMX-1 ):
                f6.write( SufHmax[i,j] )
            #enddo
        #enddo
      #endif
      return(0)#to return something
#enddef Tecout_write


