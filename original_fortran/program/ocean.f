      subroutine ocean(rtime )

      implicit none

      include 'ocean.h'
      include 'couple.h'
      include 'atmosphere.h'
      include 'seaice.h'
      include 'landice.h'
      include 'biochemical.h'
      real*8 t_state
c  subroutine for the ocean part.the ocean consist of 8 boxs as 
c  in the following diagram
c
c
c   --------------------------------------------------------
c   |        |                |                    |       |
c   |        |                |                    |       |
c   |   1  <-|-      2       -|->      3          -|-> 4   |deltaD            
c   |        |      ^         |          ^         |       |
c   ----------------|--------------------|------------------
c   |        |      |         |          |         |       |
c   |        |                |                    |       |
c   |        |                |                    |       |
c   |   5   -|->     6      <-|-       7         <-|-  8   |  D
c   |        |                |                    |       |
c   |        |                |                    |       | 
c   |        |                |                    |       |
c   -------------------------------------------------------
c     SP                       EQ                         NP
c
c     deltaL        L                    L             deltaL
c
      real*8 state, states,statet,unesco

c calculate the transport between boxes

c hezi 23.2 seperate the q to temperature driven and salinity driven parts
      
      q21=c21*( D(1)*(rho(1)-rho(2)) + D(5)*(rho(5)-rho(6)) )
      q23=c23*( D(2)*(rho(3)-rho(2)) + D(6)*(rho(7)-rho(6)) )
      q34=c34*( D(4)*(rho(4)-rho(3)) + D(8)*(rho(8)-rho(7)) )

      q21t=0.0 ! c21*( D(1)*(rhot(1)-rhot(2)) + D(5)*(rhot(5)-rhot(6)) )
      q23t=0.0 ! c23*( D(2)*(rhot(3)-rhot(2)) + D(6)*(rhot(7)-rhot(6)) )
      q34t=0.0 ! c34*( D(4)*(rhot(4)-rhot(3)) + D(8)*(rhot(8)-rhot(7)) )

      q21s=0.0 ! c21*( D(1)*(rhos(1)-rhos(2)) + D(5)*(rhos(5)-rhos(6)) )
      q23s=0.0 ! c23*( D(2)*(rhos(3)-rhos(2)) + D(6)*(rhos(7)-rhos(6)) )
      q34s=0.0 ! c34*( D(4)*(rhos(4)-rhos(3)) + D(8)*(rhos(8)-rhos(7)) )




      q56=q21
      q76=q23
      q87=q34
      q15=q21
      q48=q34

      q62=q56+q76
      q73=q87-q76


c ADVECTION

      if (q21.ge.0.0) then
         sf21=q21*s(2,k)
         sf15=q15*s(1,k)
         sf56=q56*s(5,k)

         tf21=q21*t(2,k)
         tf15=q15*t(1,k)
         tf56=q56*t(5,k)
      else
         sf21=q21*s(1,k)
         sf15=q15*s(5,k)
         sf56=q56*s(6,k)

         tf21=q21*t(1,k)
         tf15=q15*t(5,k)
         tf56=q56*t(6,k)
      endif

      if (q34.ge.0.0) then 
         sf34=q34*s(3,k)
         sf48=q48*s(4,k)
         sf87=q87*s(8,k)

         tf34=q34*t(3,k)
         tf48=q48*t(4,k)
         tf87=q87*t(8,k)
      else
         sf34=q34*s(4,k)
         sf48=q48*s(8,k)
         sf87=q87*s(7,k)

         tf34=q34*t(4,k)
         tf48=q48*t(8,k)
         tf87=q87*t(7,k)
      endif

      if (q62.ge.0.0) then
         sf62=q62*s(6,k)
         tf62=q62*t(6,k)
      else
         sf62=q62*s(2,k)
         tf62=q62*t(2,k)
      endif

      if (q23.ge.0.0) then
         sf23=q23*s(2,k)
         tf23=q23*t(2,k)
      else
         sf23=q23*s(3,k)
         tf23=q23*t(3,k)
      endif

      if (q73.ge.0.0) then
         sf73=q73*s(7,k)
         tf73=q73*t(7,k)
      else
         sf73=q73*s(3,k)
         tf73=q73*t(3,k)
      endif

      if (q76.ge.0.0) then
         sf76=q76*s(7,k)
         tf76=q76*t(7,k)
      else
         sf76=q76*s(6,k)
         tf76=q76*t(6,k)
      endif


c DIFFUSION

c vertical diffusion
c hezi vertical diffusion is twice bigger at polar boxes

c check for stable stratification:


         if (rho(1).gt.rho(5))then
       sdf15=Kv(1)*300*Oarea(1)*( s(1,k)-s(5,k) )
       tdf15=Kv(1)*300*Oarea(1)*( t(1,k)-t(5,k) )
         else
       sdf15=Kv(1)*Oarea(1)*( s(1,k)-s(5,k) )
       tdf15=Kv(1)*Oarea(1)*( t(1,k)-t(5,k) )
        endif

         if (rho(2).gt.rho(6))then
       sdf26=Kv(2)*300*Oarea(2)*( s(2,k)-s(6,k) )
       tdf26=Kv(2)*300*Oarea(2)*( t(2,k)-t(6,k) )
         else
       sdf26=Kv(2)*Oarea(2)*( s(2,k)-s(6,k) )
       tdf26=Kv(2)*Oarea(2)*( t(2,k)-t(6,k) )
        endif

         if (rho(3).gt.rho(7))then
       sdf37=Kv(3)*300*Oarea(3)*( s(3,k)-s(7,k) )
       tdf37=Kv(3)*300*Oarea(3)*( t(3,k)-t(7,k) )
         else
       sdf37=Kv(3)*Oarea(3)*( s(3,k)-s(7,k) )
       tdf37=Kv(3)*Oarea(3)*( t(3,k)-t(7,k) )
        endif

         if (rho(4).gt.rho(8))then
       sdf48=Kv(4)*300*Oarea(4)*( s(4,k)-s(8,k) )
       tdf48=Kv(4)*300*Oarea(4)*( t(4,k)-t(8,k) )
         else
       sdf48=Kv(4)*Oarea(4)*( s(4,k)-s(8,k) )
       tdf48=Kv(4)*Oarea(4)*( t(4,k)-t(8,k) )
        endif



c horizontal diffusion

      sdf12=Kh(1)*uppersection*( s(1,k)-s(2,k) )

      sdf23=Kh(2)*uppersection*( s(2,k)-s(3,k) )

      sdf34=Kh(3)*uppersection*( s(3,k)-s(4,k) )

      sdf56=Kh(4)*lowersection*( s(5,k)-s(6,k) )

      sdf76=Kh(5)*lowersection*( s(7,k)-s(6,k) )

      sdf78=Kh(6)*lowersection*( s(7,k)-s(8,k) )

      tdf12=Kh(1)*uppersection*( t(1,k)-t(2,k) )

      tdf23=Kh(2)*uppersection*( t(2,k)-t(3,k) )

      tdf34=Kh(3)*uppersection*( t(3,k)-t(4,k) )

      tdf56=Kh(4)*lowersection*( t(5,k)-t(6,k) )

      tdf76=Kh(5)*lowersection*( t(7,k)-t(6,k) )

      tdf78=Kh(6)*lowersection*( t(7,k)-t(8,k) )

c calculate seaice

      if(seaicectr) then
         call seaice
      endif



      do i=1,4

c SURFACE B. C:
c when calculating the fluxes, there are 2 option:
c restoring to specified value, or flux from atmosphere model
c accordig to the variables dataflux and coupleflux.
c only one of them can be true.
c the refference values are specified elsehwere

c update FWAO
c      do i=1,4
         FWAO(i)=FWAO(i)
c     &       (-Imelt(i)-Inewf(i))*(rhoIw/tdt)
     &       -LSIdiff(i)*(rhoIw/tdt)
     &       -LIdiff(i)*(rhoLIw/tdt)
c      end do

c calculating heat and salinity fluxes for this time step. 
c positive for heat into the ocean

      if (coupleflux) then  

c using heat and salinity fluxes from atmospheric model
c fresh-water flux from atmospheric model is converted into salinity flux
c heat flux is converted from W/m2 into m*K/s

c         do i=1,4
         tflux(i)=Hao(i)/CwrhoW
         sflux(i)=-1.*FWAO(i)*Sref
c hezi 8/2 make the precipitation falling on sea-ice to sea-ice by 
c changing FSI and not ny increasing thickness
         

         
c hezi 28/12 insert fresh-water flux perturbation
c      if ((KDAY.gt.(NDBEG+10)).and.(KDAY.lt.(NDBEG+36500))) then
c         FWAO(3)=FWAO(3)+0.75d+6
c      endif
c         end do 

      else                       ! restoring b.c.

c         do i=1,4
         tflux(i)=D(i)*gammaT*(SSTdata(i)-t(i,k))
         sflux(i)=D(i)*gammaS*(Sdata(i)-s(i,k))*Oarea(i)
c         end do
      endif


      end do

c=====================================================
c new salinity
c ====================================================

      s(1,kp)=tdt*recvol(1)*(   sf21-sf15
     &                         -sdf15-sdf12
     &              + sflux(1)     ) 
     &                                   + s(1,km)

      
      s(2,kp)=tdt*recvol(2)*(   sf62-sf21-sf23
     &                         +sdf12-sdf23-sdf26
     &              + sflux(2)      )
     &                                   +s(2,km)      

      s(3,kp)=tdt*recvol(3)*(   sf23+sf73-sf34
     &                         +sdf23-sdf34-sdf37
     &                +sflux(3)     )
     &                                   +s(3,km)

      s(4,kp)=tdt*recvol(4)*(   sf34-sf48
     &                         +sdf34-sdf48
     &              +  sflux(4)      ) 
     &                                   + s(4,km)

      s(5,kp)=tdt*recvol(5)*(   sf15-sf56   
     &                         +sdf15-sdf56         )
     &                                   +s(5,km)

      s(6,kp)=tdt*recvol(6)*(  sf56+sf76-sf62
     &                        +sdf56+sdf26+sdf76    )
     &                                   +s(6,km)

      s(7,kp)=tdt*recvol(7)*(  sf87-sf73-sf76 
     &                        -sdf76+sdf37-sdf78    )
     &                                   +s(7,km)

      s(8,kp)=tdt*recvol(8)*(   sf48-sf87 
     &                         +sdf78+sdf48         )
     &                                   +s(8,km)


c============================================
c new temperature
c===========================================
c add loosing heat because of evaporation
c      Hlh(2)=DMAX1(LH*rhoW*FWAO(2)*ac1/Aarea(2),0.)

c      Hlho(2)=LH*rhoW*FWAO(2)/CwrhoW
c      Hlho(3)=LH*rhoW*FWAO(3)/CwrhoW


      Hlho(1)=0.0
      Hlho(2)=0.0
      Hlho(3)=0.0
      Hlho(4)=0.0


cHavail(i)=CwrhoW

      t(1,kp)=tdt*recvol(1)*(   tf21-tf15
     &                         -tdf15-tdf12
     &              + tflux(1)*Oarea(1)+ Havail(1)/CwrhoW    ) 
     &                                   + t(1,km)



      t(2,kp)=tdt*recvol(2)*(   tf62-tf21-tf23
     &                         +tdf12-tdf23-tdf26+Hlho(2)
     &                +tflux(2)*Oarea(2)+ Havail(2)/CwrhoW    )
     &                                   +t(2,km)      

      t(3,kp)=tdt*recvol(3)*(   tf23+tf73-tf34
     &                         +tdf23-tdf34-tdf37+Hlho(3)
     &                +tflux(3)*Oarea(3)+ Havail(3)/CwrhoW    )
     &                                   +t(3,km)

      t(4,kp)=tdt*recvol(4)*(   tf34-tf48
     &                         +tdf34-tdf48
     &              +  tflux(4)*Oarea(4)+ Havail(4)/CwrhoW    ) 
     &                                   + t(4,km)

      t(5,kp)=tdt*recvol(5)*(   tf15-tf56 
     &                         +tdf15-tdf56    )
     &                                   +t(5,km)

      t(6,kp)=tdt*recvol(6)*(  tf56+tf76-tf62
     &                        +tdf56+tdf26+tdf76  )
     &                                   +t(6,km)

      t(7,kp)=tdt*recvol(7)*(  tf87-tf73-tf76
     &                        -tdf76+tdf37-tdf78  )
     &                                   +t(7,km)

      t(8,kp)=tdt*recvol(8)*(   tf48-tf87  
     &                         +tdf78+tdf48        )
     &                                   +t(8,km)



c ===============================================================
c calculate new density 
c ==============================================================

      do i=1,8
c         rhos(i)=states(t(i,kp),s(i,kp))
c         rhot(i)=statet(t(i,kp),s(i,kp))
         t_state=t(i,kp)-273.0
c        rho(i)=state(t_state,s(i,kp))
         rho(i)=unesco(s(i,kp),t_state)
      end do


c check for stable stratification:

c      do i=1,4
c         if (rho(i).gt.rho(i+4))then

c         t(i,kp)=(vol(i)*t(i,kp)+vol(i+4)*t(i+4,kp))/
c     &                   (vol(i)+vol(i+4))
c         s(i,kp)=(vol(i)*s(i,kp)+vol(i+4)*s(i+4,kp))/
c     &                   (vol(i)+vol(i+4))

c         AT(i,kp)=(vol(i)*AT(i,kp)+vol(i+4)*AT(i+4,kp))/
c     &                   (vol(i)+vol(i+4))

c         CT(i,kp)=(vol(i)*CT(i,kp)+vol(i+4)*CT(i+4,kp))/
c     &                   (vol(i)+vol(i+4))

c         PO4(i,kp)=(vol(i)*PO4(i,kp)+vol(i+4)*PO4(i+4,kp))/
c     &                   (vol(i)+vol(i+4))
c
c         rhos(i)=states(t(i,kp),s(i,kp))
c         rhos(i+4)=rhos(i)
c         rhot(i)=statet(t(i,kp),s(i,kp))
c         rhot(i+4)=rhot(i)
c         t_state=t(i,kp)-273.0
c         rho(i)=unesco(s(i,kp),t_state)
c         rho(i)=state(t_state,s(i,kp))
c         rho(i+4)=rho(i)
c         t(i+4,kp)=t(i,kp)
c         s(i+4,kp)=s(i,kp)
c         AT(i+4,kp)=AT(i,kp)
c         CT(i+4,kp)=CT(i,kp)
c         PO4(i+4,kp)=PO4(i,kp)

c         endif
c      enddo

      return
      end 





