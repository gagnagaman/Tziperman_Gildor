      subroutine humidity

      implicit none

      include 'atmosphere.h'
      include 'couple.h'

      real*8 A,B

c calculate the saturation vapor pressure using approximate
c Clausius-Clapeyron equation.

      
      A=2.53d+11  ! Pascall
      B=5.42d+3   ! K
      
      do i=1,4
      qsat(i)=A*exp(-B/Teta(i,k))
      enddo

c calculate new specific humidity base on constant relative
c humidity RH in each box

      do i=1,4
         q(i)=RH(i)*epsilon*qsat(i)/Ps
      enddo


      return
      end
