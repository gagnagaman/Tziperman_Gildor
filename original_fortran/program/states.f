      real*8 function states(tt,ss)

      implicit none

      include 'ocean.h'
      include 'couple.h'

      real*8 tt,ss

      states=rhoW*(1.+beta*(ss-Sref) )
      

      return
      end
