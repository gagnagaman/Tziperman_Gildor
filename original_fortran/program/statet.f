      real*8 function statet(tt,ss)

      implicit none

      include 'ocean.h'
      include 'couple.h'

      real*8 tt,ss

      statet=rhoW*(1.-alpha*(tt-Tref) )
      

      return
      end
