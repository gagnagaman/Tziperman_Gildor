
c constants and numerical values specified for the run and stay unchange:

      real*8 Gama0,albedoL,albedoI,albedoSI,Solar(4),KT,KWa,
     &       Plw0(4),Ac1,Kapa,Aarea(4),Catmos(3),caflux,albedoS,
     &       albedoC

c variables:

      real*8 x,ran,FW1,FW4,Kwa2

      real*8 fi(4,3),Plw,CO2a,Teta(4,3),Hout(4),Hin(4),Hao(4),Hlh(4),
     &       FT1122,FT2233,FT3344,FW1122,FW2233,FW3344,FT21,
     &       FT12,FT32,FT23,FT43,FT34,FW21,FW12,FW32,FW23,FW43,
     &       FW34,albedo(4),qsat(4),q(4),RH(4),FWAO(4),Ipnew(4)

      common /Aparameters1/Gama0,albedoL,albedoI,albedoSI,Solar,
     &       KT,KWa,Plw0,Ac1,Kapa,Aarea,Catmos,caflux

      common /Avar/fi,Plw,CO2a,Teta,Hout,Hin,Hao,Hlh,albedoC,
     &       FT1122,FT2233,FT3344,FW1122,FW2233,FW3344,FT21,FT12,
     &       FT32,FT23,FT43,FT34,FW21,FW12,FW32,FW23,FW43,FW34,
     &       albedo,qsat,q, RH,FWAO,Ipnew,FW1,FW4,Kwa2,albedoS

c      common /RANDOM/x,ran

c Plw0  emissivity of the atmosphere for current (320 ppm) CO2 concentration
c Ac1,Acoef1   some constant
c Kapa  - radiative forcing of CO2
c TOref(4)- ocean SST for atmosphere mode 
c KWa,KT  coefficient for meridional transport of humidity and heat
c Sref=35.   for transform fresh water flux into virtual salinity flux
c epsilon=R/Cp

c fi relative area cover by ice out of land area
c fsi relative area cover by sea-ice out of ocean area
c fo ocean fraction out of total area
c fl land fraction
c Plw emissivity of the atmosphere include CO2 effect
c CO2a CO2 concentration in the atmosphere in ppm
c Hout outgoing long-wave radiation from top of the atmosphere
c Hin  coming SW radiation on top of the atmosphere
c Solar solar constant for each box.( depends on latitude)
