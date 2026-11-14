      subroutine biochem_coef

      implicit none

      include 'ocean.h'
      include 'biochem_coef.h'
      include 'couple.h'


c this subroutine calculate the dissociation coefficients needed
c for calculation of pCO2 in biochemical.f

c it also calculate solubility of oxygen according to Weiss 1970
c this is done for the upper boxes only in a separate loop

      do n=1,8


c solubility of CO2 in seawater         [CO2*]
c                                   K0=--------
c                                        [CO2g]
c is from Weiss (1974)

c      ln(K0(n))=93.4517*(100./t(n,k))-60.2409+23.3585*ln(t(n,k)/100.)
c     $       + s(n,k)*( 0.023517-0.023656*(t(n,k)/100.) + 
c     $                  0.0047036*(t(n,k)/100.)*(t(n,k)/100.) )

      K0(n)=exp( 
     $   (9345.17/t(n,k))       -60.2409+23.3585*dlog(t(n,k)/100.)
     $  +s(n,k)*( 0.023517-0.00023656*t(n,k)  + 
     $   0.00000047036*t(n,k)*t(n,k) ))

c      K0(n)=exp(
c     $   93.4517*( 100./t(n,k) )-60.2409+23.3585*dlog(t(n,k)/100.)
c     $  +s(n,k)*( 0.023517-0.023656*( t(n,k)/100. ) + 
c     $   0.0047036*( t(n,k)/100. )*( t(n,k)/100. ) )  )


c       [H+][HCO3-]           [H+][CO3--]
c    K1=-----------    and K2=-----------  are from Millero (1987)
c         [CO2*]                [HCO3-] 
c


      K1(n)=10.**(-((6320.18/t(n,k))-126.3405+19.568*dlog(t(n,k))+
     $    (19.894-(840.39/t(n,k))-3.0189*dlog(t(n,k)))*
     $    sqrt(s(n,k))+0.00678*s(n,k)) )


c      pK2=5143.69*(1./t(n,k))-90.1833+14.613*ln(t(n,k)+
c     $    (17.176-690.59*(1./t(n,k))-2.6719*ln(t(n,k))*sqrt(s(n,k))
c     $    +0.02169*s(n,k)

      K2(n)=10.**(-((5143.69/t(n,k))-90.1833+14.613*dlog(t(n,k))+
     $    (17.176-(690.59/t(n,k))-2.6719*dlog(t(n,k)))*
     $    sqrt(s(n,k))+0.02169*s(n,k)) )


c for Boric acid,        [H+][B(OH4)-]
c                    Kb=-------------- is from Johansson & Wedborg (1982)
c                          [B(OH)3 

c      pKb=1030.5*(1./t(n,k))+5.5076-0.015469*s(n,k)+0.00015339*s(n,k)*s(n,k)

      Kb(n)=10.**(-((1030.5/t(n,k))+5.5076-0.015469*s(n,k)+
     $           1.5339d-4*s(n,k)*s(n,k) ) )



c for water. Kw=[h+]*[OH-], from Dickson & Riley (1979)

c      pKw=3441.0*(1./t(n,k))+2.241-0.09415*sqrt(s(n,k))

      Kw(n)=10.**(-((3441.0/t(n,k))+2.241-0.09415*sqrt(s(n,k)) ))

      enddo


c      do n=1,4
c loop for oxygen solubility, from Weiss (1974)
c

c      KO2(n)=exp(
c     $   255.5907*( 100./t(n,k) )-177.7888+146.4813*dlog(t(n,k)/100.)
c     $  +s(n,k)*(-0.037362+0.016504*( t(n,k)/100. ) - 
c     $   0.0020564*( t(n,k)/100. )*( t(n,k)/100. ) )  )

c      enddo


      return
      end
