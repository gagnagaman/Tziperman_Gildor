      subroutine landice(rtime )

      implicit none

      include 'ocean.h'
      include 'couple.h'
      include 'atmosphere.h'
      include 'seaice.h'
      include 'landice.h'



c compute the extension of the ice-sheets on land. ususally in the 
c following names, N stand for North pole and S for the South pole (Antartica)


         InewTemp(1)=0.0
         InewTemp(4)=0.0
         LIdiff(1)=0.0
         LIdiff(4)=0.0

*dbgz:
c      if (NMSTEP2.eq.MID_JUNE) THEN
c        SINK_PER=gamaLI*(solar(4)-502.555)*rhowLice
c      endif
           SINK_PER=0.0  ! seasonal run      

c start with South pole, i.e Antartica:

c check for melting/creating of Land Ice

c$$$       if (teta(1,k).ge.T1melt) then                 ! need to melt Land ice
c$$$
c$$$         Iflow(1)=0.009d+6*rhowLice
c$$$         Im(1)=0.03*(teta(1,k)-T1melt)*rhowLice
c$$$
c$$$c add Iflow to Icemelt
c$$$         Im(1)=Im(1)+Iflow(1)
c$$$
c$$$c make sure it doesn't shrink less the zero size
c$$$         IceVolS(kp)=DMAX1( (IceVolS(km)-tdt*im(1)),  0.d0  )
c$$$
c$$$         LIdiff(1)=IceVolS(kp)-IceVolS(km)
c$$$
c$$$         ILenS(kp)=( IceVolS(kp)/LIc1 )**(2./3.)
c$$$
c$$$
c$$$c           ======================================
c$$$       else        ! may need to create Land ice
c$$$c           ======================================
c$$$
c$$$         Im(1)=0.009d+6*rhowLice ! actually ice flow.  in m^3/s
c$$$
c$$$c calculate how much more ice is wanted to be created
c$$$
c$$$         InewTemp(1)=DMAX1(
c$$$     &     DMAX1(0.25d0,ILenS(k)/L(1))*fl(1)*(FW21+FW1)*rhowLice, 0.d0)
c$$$
c$$$c make sure it is a positive size
c$$$         IceVolS(kp)=DMAX1(IceVolS(km)+tdt*(InewTemp(1)-im(1)),0.d0)
c$$$
c$$$c make sure it does not exceed the whole box
c$$$
c$$$         IceVolS(kp)=DMIN1(IceVolS(kp),ISmax)
c$$$
c$$$         if (IceVolS(kp).gt.0.0) then
c$$$           LIdiff(1)=IceVolS(kp)-IceVolS(km)
c$$$         else
c$$$           LIdiff(1)=DMIN1(IceVolS(kp)-IceVolS(km),im(1))
c$$$         endif
c$$$
c$$$         ILenS(kp)=( IceVolS(kp)/LIc1 )**(2./3.)
c$$$
c$$$       endif
c$$$ 


c do the north pole, i.e Arctic
      

c chack for melting/creating of Land Ice

      if (teta(4,k).ge.T4melt) then    ! need to melt Land ice
      
         Iflow(4)=0.085d+6*rhowLice   ! ice flow 
         Im(4)=0.200762*(teta(4,k)-T4melt)*rhowLice


c add Iflow to Icemelt
         Im(4)=Im(4)+Iflow(4)

c make sure it doesn't shrink less the zero size
         IceVolN(kp)=DMAX1( (IceVolN(km)-tdt*im(4)),  0.d0  )

         LIdiff(4)=IceVolN(kp)-IceVolN(km)

c  calculate new half-length of the meridional extension of the sheet

         ILenN(kp)=( IceVolN(kp)/(2.*LIc4) )**(2./3.)
         LIdiff(4)=IceVolN(kp)-IceVolN(km)

      else

        Im(4)=0.085d+6*rhowLice

c calculate how much more ice is wanted to be created
        InewTemp(4)=DMAX1(
     &  DMAX1(0.25d0,ILenN(k)/L(4))*fl(4)*(FW34+FW4)*rhowLice,0.d0)

c make sure it is a positive size
        IceVolN(kp)=DMAX1(IceVolN(km)+tdt*(InewTemp(4)-im(4)),0.d0)

c make sure it does not exceed the whole box

        IceVolN(kp)=DMIN1(IceVolN(kp),INmax)
        if (IceVolN(kp).gt.0.0) then
           LIdiff(4)=IceVolN(kp)-IceVolN(km)
        else
           LIdiff(4)=DMIN1(IceVolN(kp)-IceVolN(km),im(4))
        endif
c  calculate new half-length of the meridional extension of the sheet


        ILenN(kp)=(IceVolN(kp)/(2.*LIc4) )**(2./3.)

      endif

      
c      fi(1,kp)=ILenS(kp)/L(1)
      fi(4,kp)=(2*ILenN(kp))/L(4)

      return
      end 
