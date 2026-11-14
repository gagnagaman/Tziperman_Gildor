      subroutine seaice
      implicit none
c this subroutine calculate the amount of sea ice in each box.
c more specificaly it claculate the fraction of the box coverd by 
c specified thickness as long as there are still open water in the box
c and calculate the ice thickness when all the box is cover by ice

      include 'couple.h'
      include 'atmosphere.h'
      include 'ocean.h'
      include 'seaice.h'


c check if need to create more ice in order to restore the ocean 
c temperature to Tice. 
      do i=1,4,3
         Inewf(i)=0.0
         Imelt(i)=0.0
         Havail(i)=0.0
         Hneed(i)=0.0
         LSIdiff(i)=0.

         if (t(i,k).lt.Tice) then                  ! create sea-ice

             Havail(i)=3.8d-7*CwrhoW*vol(i)*(Tice-t(i,k))   ! check amount of heat in J needed

c consider also the part of precipitation falls on the sea-ice

             sivol(i,kp)=tdt*( (Havail(i)/rhoILice)+ipnew(i)/tdt ) 
     &                                   + sivol(i,km)

             inewf(i)=Havail(i)*tdt/rhoILice

c calculate new fsi
             if ((SIvol(i,kp)/Di(i)).lt.Oarea(i)) then
              fsi(i,kp)=SIvol(i,kp)/(Di(i)*Oarea(i))
              Idepth(i)=Di(i) 
             else
              fsi(i,kp)=1.
              Idepth(i)=SIvol(i,kp)/Oarea(i)
             endif

       else  ! ocean tepereture > TiceC. may need to melt ice

              Havail(i)=3.8d-7*CwrhoW*vol(i)*(Tice-t(i,k))   ! check amount of heat in J 
              Hneed(i)=(SIvol(i,k)+ipnew(i))*rhoILice        ! heat need for melt all the ice

                if (Hneed(i).gt.abs(Havail(i))) then           ! more ice than we can melt
c make sure we do not get negative values
                   sivol(i,kp)=DMAX1(
     &             tdt*( (Havail(i)/rhoILice)+ipnew(i)/tdt ) 
     &                          + sivol(i,km),                0.d0)

                else   ! we don't need so much heat

                   Havail(i)=-Hneed(i)

                   sivol(i,kp)=DMAX1(
     &             tdt*( (Havail(i)/rhoILice)+(ipnew(i)/tdt) ) 
     &                          + sivol(i,km),               0.d0)
               
                endif

           if ((SIvol(i,kp)/Di(i)).lt.Oarea(i)) then
              fsi(i,kp)=SIvol(i,kp)/(Di(i)*Oarea(i))
              Idepth(i)=Di(i) 
           else
              fsi(i,kp)=1.
              Idepth(i)=SIvol(i,kp)/Oarea(i)
           endif


            imelt(i)=Havail(i)*tdt/rhoILice
       endif

       LSIdiff(i)=sivol(i,kp)-sivol(i,km)

       enddo

       return
       end
    
