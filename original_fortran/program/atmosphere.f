      subroutine atmosphere

      implicit none

      include 'atmosphere.h'
      include 'couple.h'
      include 'ocean.h'
      include 'seaice.h'
      include 'landice.h'
      include 'biochemical.h'

c  subroutine for the atmosphere part.the atmosphere consist of 4 boxs as 
c  in the following diagram
c
c   
c   -----|-----------|------------------|--------------|----
c   |        |                |                    |       |
c   |        |                |                    |       |
c   |   1  <-|-      2       -|->      3          -|-> 4   | Atmosphere
c   |        |                |                    |       | Height
c   ----------------|--------------------|------------------
c   |        |      |         |          |         |       |

c     SP                       EQ                         NP
c

c********************************************************************
c calculate terms which involve data from more than one boxes
c ********************************************************************

c calculate meridional heat  tranfer between the boxes in J/s using cubic 
c power of the temperature gradient
c [KT]=J*m^2/s*k^2.
c here we calculate the values without the coefficient which depends 
c on the boxes area. 
      
         FT1122=abs(KT*( Teta(2,k)-Teta(1,k) ) /
     &    ( ( L(1)+L(2) )  ))

         FT2233=abs(KT*( Teta(3,k)-Teta(2,k) ) /
     &    ( ( L(3)+L(2) )   ))

         FT3344=abs(KT*( Teta(3,k)-Teta(4,k) )/
     &    ( ( L(3)+L(4) )   ))

c to convert it from J/s into K/s, and normalized by the area of each box,
c we multiply by ac1


      if ( Teta(1,k).le.Teta(2,k) ) then
         FT21=    FT1122*ac1/Aarea(1)
         FT12=   -FT1122*ac1/Aarea(2)
      else
         FT21=   -FT1122*ac1/Aarea(1)
         FT12=    FT1122*ac1/Aarea(2)
      endif

      if ( Teta(2,k).le.Teta(3,k) ) then
         FT32=   FT2233*ac1/Aarea(2)
         FT23=  -FT2233*ac1/Aarea(3)
      else
         FT32=  -FT2233*ac1/Aarea(2)
         FT23=   FT2233*ac1/Aarea(3)
      endif

      if ( Teta(3,k).le.Teta(4,k) ) then
         FT43=    FT3344*ac1/Aarea(3)
         FT34=   -FT3344*ac1/Aarea(4)
      else
         FT43=   -FT3344*ac1/Aarea(3)
         FT34=    FT3344*ac1/Aarea(4)
      endif


c =====================================================================
c calculate new humidities for calculating meridional fresh water flux
c =====================================================================

      call humidity

c calculate meridional moisture fluxes in Sv
c hezi 20.2.98 



      if ( Teta(1,k).le.Teta(2,k) ) then
         FW1122=KWa*( Teta(2,k)-Teta(1,k) )*q(1) /
     &           ( 0.25*( L(1)+L(2) )  )
         FW21=    FW1122
         FW12=   -FW1122
      else
         FW1122=KWa*( Teta(1,k)-Teta(2,k) )*q(2) /
     &           ( 0.25*( L(1)+L(2) )  )
         FW21=   -FW1122
         FW12=    FW1122
      endif

      if ( Teta(2,k).le.Teta(3,k) ) then
         FW2233=KWa*( Teta(3,k)-Teta(2,k) )*q(2) /
     &      ( 0.25*( L(3)+L(2) )  )
         FW32=   FW2233
         FW23=  -FW2233
      else
         FW2233=KWa*( Teta(2,k)-Teta(3,k) )*q(3) /
     &      ( 0.25*( L(3)+L(2) )  )
         FW32=  -FW2233
         FW23=   FW2233
      endif

      if ( Teta(3,k).le.Teta(4,k) ) then
         FW3344=KWa*( Teta(4,k)-Teta(3,k) )*q(3) /
     &      ( 0.25*( L(3)+L(4) )  )
         FW43=    FW3344
         FW34=   -FW3344
      else
         FW3344=KWa*( Teta(3,k)-Teta(4,k) )*q(4) /
     &      ( 0.25*( L(3)+L(4) )  )
         FW43=   -FW3344
         FW34=    FW3344
      endif


c hezi 15.2 
      FW21=0.5*FW21
      FW12=0.5*FW12
      FW43=0.5*FW43
      FW34=0.5*FW34


      FW4=KWa2*fo(4)*(1.-fsi(4,k))*q(4)
      FW1=KWa2*fo(1)*(1.-fsi(1,k))*q(1)




c the water flux from atmosphere INTO ocean in Sv. these will be
c converted into virtual salt flux when run in coupled mode

c hezi 8/2 correct for precipitation falling on sea-ice
c second part of the correction is in seaice.f
      if(seaicectr) then
      Ipnew(1)=DMAX1((FW21+FW1)*fsi(1,k)*(1.-fl(1))*
     &         (rhoW/rhoice)*tdt,0.d0)
      Ipnew(4)=DMAX1((FW34+FW4)*fsi(4,k)*(1.-fl(4))*
     &         (rhoW/rhoice)*tdt,0.d0)
      endif

      FWAO(1)=FW21
      FWAO(2)=FW12+FW32
      FWAO(3)=FW23+FW43
      FWAO(4)=FW34

c calculate heat flux cause by latent heat. i evaporation region, the atmosphere 
c loose heat and in region with precipitation, the atmosphere gain heat
c to convert it from J/s into K/s, normalized by the area of each box,
c we multiply by ac1


      Hlh(1)=0.
      Hlh(2)=0.
      Hlh(3)=0.
      Hlh(4)=0.


c calculate land-ice

      if(landicectr) then
        call landice
      endif



      


c ***********************************************************************
c calculate terms which need data form one box only 
c i.e the next loop can be distributed between processors
c ************************************************************************

c======================================================================        
c  calculating boxes albedo
c  fl=amount of land, fi=amount of ice out of the amount of land, 
c  fo=amount of ocean,fsi=amount of sea-ice out of the amount of ocean
c  (albedoL,albedoI,albedoSI)= (land , ice, sea-ice) albedo
c
c in the first stage, we assume the f's are fixed. later we will use 
c ice and sea-ice model to calculate them.
c======================================================================



      do i=1,4

         albedo(i)=fl(i)*(1-fi(i,k))*albedoL+
     &             fi(i,k)*fl(i)*albedoI+
     &             fsi(i,k)*fo(i)*albedoSI+
     &             (1-fsi(i,k))*fo(i)*albedoS




c calculate the income radiation on top of the atmosphere
c the 0.7 account for cloud albedo of 0.3

         Hin(i)=Solar(i)*(1-albedo(i) )*(1.-albedoC)

c calculate outgoing radiation from the atmosphere, using Teta 
c hezi: chack the diffrence when using temeprature at 100mb

c calculate the emissivity. its depends on the [CO2]. we assume [CO2]
c is homogeneous in the atmosphere i.e does not depends on the box
c concentration of CO2 is in ppm. CO2ref=320 ppm

c later we will calculate the actual concentration, CO2t
c      CO2a=320. ! we start with constant concentration

c hezi comment the following line till use biology
c         Plw=Plw0(i)-Kapa*alog(real(CO2a/CO2ref))
         Plw=Plw0(i)

         Hout(i)=sigma*Plw*Teta(i,k)*Teta(i,k)*Teta(i,k)*Teta(i,k)
c      enddo




c====================================================================
c SURFACE B. C: heat fluxes in/out of the ocean
c====================================================================
c when calculating the fluxes, there are 2 option:
c (a) flux derived from ocean temperature calculated by the ocean model
c (b) flux derived from given ocean temperature
c 
c accordig to the variables coupledflux, true for (a).
c the refference values are specified in the job_cntrl

c calculating fluxes into the ocean for this time step. positive 
c for heat into the ocean

      if (coupleflux) then
c       do i=1,4
       Hao(i)=CwRhoW*D(i)*fo(i)*(1.-fsi(i,k))*Gama0*(Teta(i,k)-t(i,k))
     & +(CwRhoW*D(i)*fo(i)*fsi(i,k)*0.05*Gama0*(Teta(i,k)-t(i,k))
     &   /(Idepth(i) ))
c       end do
      else
c       do i=1,4
       Hao(i)=CwRhoW*D(i)*fo(i)*(1.-fsi(i,k))*Gama0*
     &  (Teta(i,k)-SSTdata(i))
     & +(CwRhoW*D(i)*fo(i)*fsi(i,k)*0.05*Gama0*(Teta(i,k)-SSTdata(i))
     &   /(Idepth(i) ))

c       end do

      endif

      enddo



c======================================================
c  step in time
c======================================================

c calculate the new Potencial Temperature
c
c   Acoef1=2^(R/Cp)*g/(Ps*Cp)
c

      Teta(1,kp)=tdt*(FT21+Hlh(1)+Ac1*(-Hao(1)-Hout(1)+Hin(1) ) )
     &           +Teta(1,km)

      Teta(2,kp)=tdt*(FT12+Hlh(2)+FT32+Ac1*(-Hao(2)-Hout(2)+Hin(2) ) )
     &           +Teta(2,km)

      Teta(3,kp)=tdt*(FT43+Hlh(3)+FT23+Ac1*(-Hao(3)-Hout(3)+Hin(3) ) )
     &           +Teta(3,km)

      Teta(4,kp)=tdt*(FT34+Hlh(4)+Ac1*(-Hao(4)-Hout(4)+Hin(4) ) )
     &           +Teta(4,km)

      return
      end 





