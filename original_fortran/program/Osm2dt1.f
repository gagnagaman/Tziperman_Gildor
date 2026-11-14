
      subroutine Osm2dt1

      implicit none

      include 'ocean.h'
      include 'couple.h'
      include 'seaice.h'
      include 'biochemical.h'
      include 'landice.h'
      include 'atmosphere.h'

c  this routine averages fields every nsmth time steps to eliminate 
c  time-step splitting caused by the leap-frog scheme.



      do  i=1,8

         t(i,k)=t(i,k)+0.25*(t(i,kp)-2*t(i,k)+t(i,km))
         s(i,k)=s(i,k)+0.25*(s(i,kp)-2*s(i,k)+s(i,km))

      end do


      do i=1,4,3
         sivol(i,k)=sivol(i,k)+
     &        0.25*( sivol(i,kp)-2*sivol(i,k)+sivol(i,km))
         fsi(i,k)=fsi(i,k)+
     &        0.25*( fsi(i,kp)-2*fsi(i,k)+fsi(i,km))
         fi(i,k)=fi(i,k)+
     &        0.25*( fi(i,kp)-2*fi(i,k)+fi(i,km))

      enddo

        IceVolS(k)=IceVolS(k)+
     &      0.25*(IceVolS(kp)-2*IceVolS(k)+IceVolS(km))
        IceVolN(k)=IceVolN(k)+
     &      0.25*(IceVolN(kp)-2*IceVolN(k)+IceVolN(km))
        IlenS(k)=IlenS(k)+0.25*(IlenS(kp)-2*IlenS(k)+IlenS(km))  
        IlenN(k)=IlenN(k)+0.25*(IlenN(kp)-2*IlenN(k)+IlenN(km))  

      return
      end

