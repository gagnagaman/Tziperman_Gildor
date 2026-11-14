
      subroutine Asm2dt1

      implicit none

      include 'atmosphere.h'
      include 'couple.h'
      include 'biochemical.h'


c  this routine averages fields every nsmth time steps to eliminate 
c  time-step splitting caused by the leap-frog scheme.

      do  i=1,4

       teta(i,k) = teta(i,k)+0.25*( teta(i,kp)-2*teta(i,k )+teta(i,km))


       end do

       Catmos(k) = Catmos(k)+
     &          0.25*( Catmos(kp)-2*Catmos(k )+Catmos(km))

       return
       end
