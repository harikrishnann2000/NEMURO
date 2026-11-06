c_ ---------------------------------------------------------------------
c_ RCS lines preceded by "c_ "
c_ ---------------------------------------------------------------------
c_
c_ $Source: /home/orr/WWW/Biotic/boundcond/RCS/scco2.f,v $ 
c_ $Revision: 1.1 $
c_ $Date: 1999/03/22 12:57:48 $   ;  $State: Exp $
c_ $Author: orr $ ;  $Locker:  $
c_
c_ ---------------------------------------------------------------------
c_ $Log: scco2.f,v $
c_ Revision 1.1  1999/03/22 12:57:48  orr
c_ Initial revision
c_
c_ ---------------------------------------------------------------------
c_ 
      function scco2(t)
c
c*********************************************************************
c
c  Computes the Schmidt number of CO2 in seawater using the 
c  formulation presented by Wanninkhof (1992, J. Geophys. Res., 97,
c  7373-7382).  Input is temperature in deg C.
c
c*********************************************************************
c
      scco2 = 2073.1 - 125.62*t + 3.6276*t**2 - 0.043219*t**3
c
      return
      end
