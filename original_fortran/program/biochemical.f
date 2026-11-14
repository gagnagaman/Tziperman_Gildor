      subroutine biochemical

      implicit none  

      real*8 state

      include 'atmosphere.h'
      include 'ocean.h'
      include 'couple.h'
      include 'seaice.h'
      include 'landice.h'
      include 'biochemical.h'
      include 'biochem_coef.h'

c this subroutine calculate the pCO2 in water. 

c than calculate new Alk, PO4, O2, Total CO2



      do i=1,8

c calculate BT using the salinity, from Uppstorm (1974) .
c i took it from Yamanaka and there seems to be a problem

c      BT(i)=(0.000232/(10.81*1.80655))*s(i,k)*1.d+6

      BT(i)=11.87988687*s(i,k)

c solve for a
c the solver for a is in iteration method, as depicted in Yamanaka (1996)
c based on Peng et al. (1987) in Tellus
c we stop when the difference between two succesive values is reduce to
c less than PREC

      PREC=1.d-3
      a_temp1=a(i)   ! initial guess


 40   continue

      AC_temp=AT(i,k)-
     %     ( a_temp1*Kb(i)*BT(i) )/(1.+a_temp1*Kb(i) ) -
     %          ( a_temp1*Kw(i)-(1./a_temp1) )

c getting a better a

      gama_t=AC_temp/CT(i,k)
      a_temp2=( 
     $          -(1.-gama_t)*K1(i)+
     $ sqrt( 
     $       (1-gama_t)*(1-gama_t)*K1(i)*K1(i)
     $      +4*(2-gama_t)*K1(i)*K2(i)*gama_t     )   )/
     $  ( 2*(2-gama_t)*K1(i)*K2(i) )

      if ( abs((a_temp2-a_temp1)/a_temp2).le.PREC) then
         a(i)=a_temp2
         goto 50
      else
         a_temp1=a_temp2
         goto 40
      endif

 50   continue

c calculate concentration of the different species. for testing purpose,we
c do it for all species for the time being. later we can do it only for
c [CO2*] in order to reduce computing time.


      ct1=1+a(i)*K1(i)+a(i)*a(i)*K1(i)*K2(i)
      ct2=1+a(i)*Kb(i)

      CO2(i) =CT(i,k)/ct1
      HCO3(i)=a(i)*K1(i)*CT(i,k)/ct1
      CO3(i) =a(i)*a(i)*K1(i)*K2(i)*CT(i,k)/ct1

      BOH3(i)=BT(i)/ct2
      BOH4(i)=a(i)*Kb(i)*BT(i)/ct2


      enddo


c calculate flux of CO2 between the ocean and the atmosphere.
c later need to add also oxygen transfer
c 
c change atmospheric concentartion from ppm into mole/KgWater. than using
c K0, calculate the corresponding concentration in water. the flux is
c proportional to the difference between the last concentration and the 
c actual concentration of carbon in the water. positive flux into the ocean

      caflux=0.0 ! used for calculating new atmospheric concentration

      do i=1,4
c         CO2s(i)=(K0(i)*CO2a(k))/( R*teta(i,k)*rho(i) )
         CO2s(i)=K0(i)*CO2a
         CTflux(i)=(pv(i)*Oarea(i)*(1.-fsi(i,k)))*(CO2s(i)-CO2(i))
         caflux=caflux+CTflux(i)*1028.d-6 ! rhoW*1.d-6  1.d-6 from micromole into mole
      enddo

c calculate new atmospheric concentration in ppm

      Catmos(kp)=Catmos(km)-
     &           tdt*caflux

      CO2a=Catmos(kp)/1.77d+14  


c calculate CT (total carbon) , AT (total alkalinity), PO4 and O2
c use the transport between the boxes which is calculated in ocean.f

c advective fluxes

      if (q21.ge.0.0) then
         CTf21=q21*CT(2,k)
         CTf15=q15*CT(1,k)
         CTf56=q56*CT(5,k)

         PO4f21=q21*PO4(2,k)
         PO4f15=q15*PO4(1,k)
         PO4f56=q56*PO4(5,k)


         ATf21=q21*AT(2,k)
         ATf15=q15*AT(1,k)
         ATf56=q56*AT(5,k)
      else
         CTf21=q21*CT(1,k)
         CTf15=q15*CT(5,k)
         CTf56=q56*CT(6,k)

         PO4f21=q21*PO4(1,k)
         PO4f15=q15*PO4(5,k)
         PO4f56=q56*PO4(6,k)


         ATf21=q21*AT(1,k)
         ATf15=q15*AT(5,k)
         ATf56=q56*AT(6,k)
      endif

      if (q34.ge.0.0) then 
         CTf34=q34*CT(3,k)
         CTf48=q48*CT(4,k)
         CTf87=q87*CT(8,k)

         PO4f34=q34*PO4(3,k)
         PO4f48=q48*PO4(4,k)
         PO4f87=q87*PO4(8,k)

         ATf34=q34*AT(3,k)
         ATf48=q48*AT(4,k)
         ATf87=q87*AT(8,k)
      else
         CTf34=q34*CT(4,k)
         CTf48=q48*CT(8,k)
         CTf87=q87*CT(7,k)

         PO4f34=q34*PO4(4,k)
         PO4f48=q48*PO4(8,k)
         PO4f87=q87*PO4(7,k)

         ATf34=q34*AT(4,k)
         ATf48=q48*AT(8,k)
         ATf87=q87*AT(7,k)
      endif

      if (q62.ge.0.0) then
         CTf62=q62*CT(6,k)
         ATf62=q62*AT(6,k)
         PO4f62=q62*PO4(6,k)
      else
         CTf62=q62*CT(2,k)
         ATf62=q62*AT(2,k)
         PO4f62=q62*PO4(2,k)
      endif

      if (q23.ge.0.0) then
         CTf23=q23*CT(2,k)
         ATf23=q23*AT(2,k)
         PO4f23=q23*PO4(2,k)
      else
         CTf23=q23*CT(3,k)
         ATf23=q23*AT(3,k)
         PO4f23=q23*PO4(3,k)
      endif

      if (q73.ge.0.0) then
         CTf73=q73*CT(7,k)
         ATf73=q73*AT(7,k)
         PO4f73=q73*PO4(7,k)
      else
         CTf73=q73*CT(3,k)
         ATf73=q73*AT(3,k)
         PO4f73=q73*PO4(3,k)
      endif

      if (q76.ge.0.0) then
         CTf76=q76*CT(7,k)
         ATf76=q76*AT(7,k)
         PO4f76=q76*PO4(7,k)
      else
         CTf76=q76*CT(6,k)
         ATf76=q76*AT(6,k)
         PO4f76=q76*PO4(6,k)
      endif

c diffussive fluxes

c vertical diffusion

c hezi: vertical diffusion is twice bigger at polar boxes
c      if(Kv().gt.0.0) then ! check for vertical diffusion

c$$$         if (rho(1).gt.rho(5))then
c$$$       CTdf15=2.5d-6*Oarea(1)*( CT(1,k)-CT(5,k) )
c$$$       ATdf15=2.5d-6*Oarea(1)*( AT(1,k)-AT(5,k) )
c$$$       PO4df15=2.5d-6*Oarea(1)*( PO4(1,k)-PO4(5,k) )
c$$$         else
       CTdf15=dmax1(dmin1( Kv(1)/sqrt
     &  (dmax1( rho(5)-rho(1)+dsigma ,1.d-5)),
     &   7.5d-6),2.7d-8)*Oarea(1)*( CT(1,k)-CT(5,k) )
       ATdf15=dmax1(dmin1( Kv(1)/sqrt
     &  (dmax1( rho(5)-rho(1)+dsigma ,1.d-5)),
     &   7.5d-6),2.7d-8)*Oarea(1)*( AT(1,k)-AT(5,k) )
       PO4df15=dmax1(dmin1( Kv(1)/sqrt
     &  (dmax1( rho(5)-rho(1)+dsigma ,1.d-5)),
     &   7.5d-6),2.7d-8)*Oarea(1)*( PO4(1,k)-PO4(5,k) )
      
c$$$        endif
c$$$
c$$$         if (rho(2).gt.rho(6))then
c$$$       CTdf26=2.5d-6*Oarea(2)*( CT(2,k)-CT(6,k) )
c$$$       ATdf26=2.5d-6*Oarea(2)*( AT(2,k)-AT(6,k) )
c$$$       PO4df26=2.5d-6*Oarea(2)*( PO4(2,k)-PO4(6,k) )       
c$$$         else
c$$$       CTdf26=dmax1(dmin1( Kv(2)/sqrt( rho(6)-rho(2) )  ,
c$$$     &   7.5d-6),0.7d-8)*Oarea(2)*( CT(2,k)-CT(6,k) )
c$$$       ATdf26=dmax1(dmin1( Kv(2)/sqrt( rho(6)-rho(2) ) ,
c$$$     &   7.5d-6),0.7d-8)*Oarea(2)*( AT(2,k)-AT(6,k) )
c$$$       PO4df26=dmax1(dmin1( Kv(2)/sqrt( rho(6)-rho(2) ) ,
c$$$     &   7.5d-6),0.7d-8)*Oarea(2)*( PO4(2,k)-PO4(6,k) )
       CTdf26=1.d+6*( CT(2,k)-CT(6,k) )
       ATdf26=1.d+6*( AT(2,k)-AT(6,k) )
       PO4df26=1.d+6*( PO4(2,k)-PO4(6,k) )      
      
c$$$        endif
c$$$
c$$$
c$$$         if (rho(3).gt.rho(7))then
c$$$       ctdf37=2.5d-6*Oarea(3)*( ct(3,k)-ct(7,k) )
c$$$       atdf37=2.5d-6*Oarea(3)*( at(3,k)-at(7,k) )
c$$$       po4df37=2.5d-6*Oarea(3)*( po4(3,k)-po4(7,k) )
c$$$         else
c$$$       ctdf37=dmax1(dmin1( Kv(3)/sqrt( rho(7)-rho(3) )  ,
c$$$     &   7.5d-6),0.7d-8)*Oarea(3)*( ct(3,k)-ct(7,k) )
c$$$       atdf37=dmax1(dmin1( Kv(3)/sqrt( rho(7)-rho(3) )  ,
c$$$     &   7.5d-6),0.7d-8)*Oarea(3)*( at(3,k)-at(7,k) )
c$$$       po4df37=dmax1(dmin1( Kv(3)/sqrt( rho(7)-rho(3) ) , 
c$$$     &   7.5d-6),0.7d-8)*Oarea(3)*( po4(3,k)-po4(7,k) )
       ctdf37=1.d+6*( ct(3,k)-ct(7,k) )
       atdf37=1.d+6*( at(3,k)-at(7,k) )
       po4df37=1.d+6*( po4(3,k)-po4(7,k) )
       
c$$$        endif
       

c$$$         if (rho(4).gt.rho(8))then
c$$$       CTdf48=2.5d-6*Oarea(4)*( CT(4,k)-CT(8,k) )
c$$$       ATdf48=2.5d-6*Oarea(4)*( AT(4,k)-AT(8,k) )
c$$$       PO4df48=2.5d-6*Oarea(4)*( PO4(4,k)-PO4(8,k) )
c$$$         else
c$$$       CTdf48=dmax1(dmin1( Kv(4)/sqrt( rho(8)-rho(4)+dsigma )  ,
c$$$     &   7.5d-6),0.7d-8)*Oarea(4)*( CT(4,k)-CT(8,k) )
c$$$       ATdf48=dmax1(dmin1( Kv(4)/sqrt( rho(8)-rho(4)+dsigma ) ,
c$$$     &   7.5d-6),0.7d-8)*Oarea(4)*( AT(4,k)-AT(8,k) )
c$$$       PO4df48=dmax1(dmin1( Kv(4)/sqrt( rho(8)-rho(4)+dsigma )  ,
c$$$     &   7.5d-6),0.7d-8)*Oarea(4)*( PO4(4,k)-PO4(8,k) )
c       CTdf48=1.d+6*( CT(4,k)-CT(8,k) )
c       ATdf48=1.d+6*( AT(4,k)-AT(8,k) )
c       PO4df48=1.d+6*( PO4(4,k)-PO4(8,k) )

       CTdf48=dmax1(dmin1( Kv(1)/sqrt
     &  (dmax1( rho(8)-rho(4)+dsigma2 ,1.d-5)),
     &   7.5d-6),2.7d-8)*Oarea(1)*( CT(4,k)-CT(8,k) )
       ATdf48=dmax1(dmin1( Kv(1)/sqrt
     &  (dmax1( rho(8)-rho(4)+dsigma2 ,1.d-5)),
     &   7.5d-6),2.7d-8)*Oarea(1)*( AT(4,k)-AT(8,k) )
       PO4df48=dmax1(dmin1( Kv(1)/sqrt
     &  (dmax1( rho(8)-rho(4)+dsigma2 ,1.d-5)),
     &   7.5d-6),2.7d-8)*Oarea(1)*( PO4(4,k)-PO4(8,k) )

       
c$$$        endif
c horizontal diffusion

      CTdf12=Kh(1)*uppersection*( CT(1,k)-CT(2,k) )

      CTdf23=Kh(2)*uppersection*( CT(2,k)-CT(3,k) )

      CTdf34=Kh(3)*uppersection*( CT(3,k)-CT(4,k) )

      CTdf56=Kh(4)*lowersection*( CT(5,k)-CT(6,k) )

      CTdf76=Kh(5)*lowersection*( CT(7,k)-CT(6,k) )

      CTdf78=Kh(6)*lowersection*( CT(7,k)-CT(8,k) )

      ATdf12=Kh(1)*uppersection*( AT(1,k)-AT(2,k) )

      ATdf23=Kh(2)*uppersection*( AT(2,k)-AT(3,k) )

      ATdf34=Kh(3)*uppersection*( AT(3,k)-AT(4,k) )

      ATdf56=Kh(4)*lowersection*( AT(5,k)-AT(6,k) )

      ATdf76=Kh(5)*lowersection*( AT(7,k)-AT(6,k) )

      ATdf78=Kh(6)*lowersection*( AT(7,k)-AT(8,k) )

      PO4df12=Kh(1)*uppersection*( PO4(1,k)-PO4(2,k) )

      PO4df23=Kh(2)*uppersection*( PO4(2,k)-PO4(3,k) )

      PO4df34=Kh(3)*uppersection*( PO4(3,k)-PO4(4,k) )

      PO4df56=Kh(4)*lowersection*( PO4(5,k)-PO4(6,k) )

      PO4df76=Kh(5)*lowersection*( PO4(7,k)-PO4(6,k) )

      PO4df78=Kh(6)*lowersection*( PO4(7,k)-PO4(8,k) )

c ==================================================
c  export production+calcification
c==================================================

      totProd=0.0   ! total EP
      sProd  =0.0
      hProd  =0.0

c hezi: chabge Redfield ratio to be larger with higher stratification
c      RedC(2)=120+1000*(rho(6)-rho(2)-3.47)
c      RedC(3)=RedC(2)

       do i=1,4
c pc1 include normalization to the light-factor and 0.1
c for the fact that about 10% of the production sink below 
c 400 m (Suess 1980)

         PROD(i)=pc1(i)*1.d-6*D(i)*Oarea(i)*(1.-fsi(i,k))*solar(i)*
     &    PO4(i,k)*PO4(i,k)/(h+PO4(i,k))
        
         PRODC(i)=-RedC(i)*PROD(i)
         ATprod(i)=RedN*PROD(i)
c calculate temperature-dependent rain ratio as in Maier-Reimer 1993
c i change in the exponent from 0.1 to 0.05
         Rain(i)=61*exp(0.1*(t(i,k)-283.))/
     &            (1.+exp(0.1*(t(i,k)-283.)))
         CaCo3C(i)=Rain(i)*PRODC(i)/RedC(i)
         ATcaco3(i)=2.*CaCo3C(i)

c       sProd=sProd-PRODC(i)*rho(i)*1.d-6 
c       hProd=hProd-CaCo3C(i)*rho(i)*1.d-6

       sProd=sProd-PRODC(i)
       hProd=hProd-CaCo3C(i)        

       enddo
       totProd=totProd+sProd+hProd
      
       do i=5,8
            PRODC(i)=RedC(i-4)*PROD(i-4)
            PROD(i)=PROD(i-4)
            ATprod(i)=-RedN*PROD(i)
            CaCo3C(i)=-CaCo3C(i-4)
            ATcaco3(i)=2.*CaCo3C(i)
        enddo



c=====================================================
c new Total Carbon
c ====================================================

      CT(1,kp)=tdt*recvol(1)*(   CTf21-CTf15
     &                         -CTdf15-CTdf12
     &              + CTflux(1)+
     &    PRODC(1)+CaCo3C(1)      ) 
     &                                   + CT(1,km)

      
      CT(2,kp)=tdt*recvol(2)*(   CTf62-CTf21-CTf23
     &                         +CTdf12-CTdf23-CTdf26
     &              + CTflux(2)+
     &    PRODC(2)+CaCo3C(2)     )
     &                                   +CT(2,km)      

      CT(3,kp)=tdt*recvol(3)*(   CTf23+CTf73-CTf34
     &                         +CTdf23-CTdf34-CTdf37
     &                +CTflux(3)+
     &    PRODC(3)+CaCo3C(3)      )
     &                                   +CT(3,km)

      CT(4,kp)=tdt*recvol(4)*(   CTf34-CTf48
     &                         +CTdf34-CTdf48
     &              +  CTflux(4)+
     &    PRODC(4)+CaCo3C(4)       ) 
     &                                   + CT(4,km)

      CT(5,kp)=tdt*recvol(5)*(   CTf15-CTf56   
     &                         +CTdf15-CTdf56+
     &    PRODC(5)+CaCo3C(5)        )
     &                                   +CT(5,km)

      CT(6,kp)=tdt*recvol(6)*(  CTf56+CTf76-CTf62
     &                        +CTdf56+CTdf26+CTdf76+
     &    PRODC(6)+CaCo3C(6)        )
     &                                   +CT(6,km)

      CT(7,kp)=tdt*recvol(7)*(  CTf87-CTf73-CTf76 
     &                        -CTdf76+CTdf37-CTdf78+
     &    PRODC(7)+CaCo3C(7)       )
     &                                   +CT(7,km)

      CT(8,kp)=tdt*recvol(8)*(   CTf48-CTf87 
     &                         +CTdf78+CTdf48+
     &    PRODC(8)+CaCo3C(8)         )
     &                                   +CT(8,km)


c============================================
c new Total Alkalinity
c===========================================

      AT(1,kp)=tdt*recvol(1)*(   ATf21-ATf15
     &                         -ATdf15-ATdf12
     &        +ATprod(1)+ATcaco3(1)                 ) 
     &                                   + AT(1,km)

      AT(2,kp)=tdt*recvol(2)*(   ATf62-ATf21-ATf23
     &                         +ATdf12-ATdf23-ATdf26
     &            + ATprod(2)+ATcaco3(2)             )
     &                                   +AT(2,km)      

      AT(3,kp)=tdt*recvol(3)*(   ATf23+ATf73-ATf34
     &                         +ATdf23-ATdf34-ATdf37
     &       +ATprod(3)+ATcaco3(3)                    )
     &                                   +AT(3,km)

      AT(4,kp)=tdt*recvol(4)*(   ATf34-ATf48
     &                         +ATdf34-ATdf48
     &        +ATprod(4)+ATcaco3(4)                ) 
     &                                   + AT(4,km)

      AT(5,kp)=tdt*recvol(5)*(   ATf15-ATf56 
     &                         +ATdf15-ATdf56   
     &       +ATprod(5)+ATcaco3(5) )
     &                                   +AT(5,km)

      AT(6,kp)=tdt*recvol(6)*(  ATf56+ATf76-ATf62
     &                        +ATdf56+ATdf26+ATdf76
     &       +ATprod(6)+ATcaco3(6)  )
     &                                   +AT(6,km)

      AT(7,kp)=tdt*recvol(7)*(  ATf87-ATf73-ATf76
     &                        -ATdf76+ATdf37-ATdf78 
     &       +ATprod(7)+ATcaco3(7) )
     &                                   +AT(7,km)

      AT(8,kp)=tdt*recvol(8)*(   ATf48-ATf87  
     &                         +ATdf78+ATdf48      
     &       +ATprod(8)+ATcaco3(8)  )
     &                                   +AT(8,km)

c =========================================================
c    new phosphate
c =========================================================
      PO4(1,kp)=tdt*recvol(1)*(   PO4f21-PO4f15
     &                         -PO4df15-PO4df12
     &        -prod(1)                 ) 
     &                                   + PO4(1,km)

      PO4(2,kp)=tdt*recvol(2)*(   PO4f62-PO4f21-PO4f23
     &                         +PO4df12-PO4df23-PO4df26
     &        - prod(2)                )
     &                                   +PO4(2,km)      

      PO4(3,kp)=tdt*recvol(3)*(   PO4f23+PO4f73-PO4f34
     &                         +PO4df23-PO4df34-PO4df37
     &       - prod(3)                   )
     &                                   +PO4(3,km)

      PO4(4,kp)=tdt*recvol(4)*(   PO4f34-PO4f48
     &                         +PO4df34-PO4df48
     &        - prod(4)                ) 
     &                                   + PO4(4,km)

      PO4(5,kp)=tdt*recvol(5)*(   PO4f15-PO4f56 
     &                         +PO4df15-PO4df56   
     &       + prod(5)    )
     &                                   +PO4(5,km)

      PO4(6,kp)=tdt*recvol(6)*(  PO4f56+PO4f76-PO4f62
     &                        +PO4df56+PO4df26+PO4df76
     &       + prod(6)    )
     &                                   +PO4(6,km)

      PO4(7,kp)=tdt*recvol(7)*(  PO4f87-PO4f73-PO4f76
     &                        -PO4df76+PO4df37-PO4df78 
     &       + prod(7)   )
     &                                   +PO4(7,km)

      PO4(8,kp)=tdt*recvol(8)*(   PO4f48-PO4f87  
     &                         +PO4df78+PO4df48      
     &       + prod(8)    )
     &                                   +PO4(8,km)



c mixed the two mid-latitude boxes to get the effect of one box

      CT(2,kp)=(vol(2)*CT(2,kp)+vol(3)*CT(3,kp))/(vol(2)+vol(3))
      CT(3,kp)=CT(2,kp)

      CT(6,kp)=(vol(6)*CT(6,kp)+vol(7)*CT(7,kp))/(vol(6)+vol(7))
      CT(7,kp)=CT(6,kp)



      AT(2,kp)=(vol(2)*AT(2,kp)+vol(3)*AT(3,kp))/(vol(2)+vol(3))
      AT(3,kp)=AT(2,kp)

      AT(6,kp)=(vol(6)*AT(6,kp)+vol(7)*AT(7,kp))/(vol(6)+vol(7))
      AT(7,kp)=AT(6,kp)



      PO4(2,kp)=(vol(2)*PO4(2,kp)+vol(3)*PO4(3,kp))/(vol(2)+vol(3))
      PO4(3,kp)=PO4(2,kp)

      PO4(6,kp)=(vol(6)*PO4(6,kp)+vol(7)*PO4(7,kp))/(vol(6)+vol(7))
      PO4(7,kp)=PO4(6,kp)     

      

      return
      end
