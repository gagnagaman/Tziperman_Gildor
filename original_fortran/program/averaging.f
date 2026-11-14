      subroutine averaging
c this subroutine calculate avereged means

      implicit none
      include 'averaging.h'
      include 'insolation1.h'
      include 'atmosphere.h'
      include 'ocean.h'
      include 'couple.h'
      include 'seaice.h'
      include 'landice.h'
      include 'biochemical.h'
      include 'biochem_coef.h'

c  assume we already zeroed all terms
      
      j=13
      
      if (NMSTEP.le.NJM1S(1)) then
         j=1
       endif  ! of month 1
      if (NMSTEP.le.NJM1S(2).and.NMSTEP.gt.NJM1S(1)) then
         j=2
       endif  ! of month 2
      if (NMSTEP.le.NJM1S(3).and.NMSTEP.gt.NJM1S(2)) then
         j=3
       endif  ! of month 3
      if (NMSTEP.le.NJM1S(4).and.NMSTEP.gt.NJM1S(3)) then
         j=4
       endif  ! of month 4
      if (NMSTEP.le.NJM1S(5).and.NMSTEP.gt.NJM1S(4)) then
         j=5
       endif  ! of month 5
      if (NMSTEP.le.NJM1S(6).and.NMSTEP.gt.NJM1S(5)) then
         j=6
       endif  ! of month 6
      if (NMSTEP.le.NJM1S(7).and.NMSTEP.gt.NJM1S(6)) then
         j=7
       endif  ! of month 7
      if (NMSTEP.le.NJM1S(8).and.NMSTEP.gt.NJM1S(7)) then
         j=8
       endif  ! of month 8
      if (NMSTEP.le.NJM1S(9).and.NMSTEP.gt.NJM1S(8)) then
         j=9
       endif  ! of month 9
      if (NMSTEP.le.NJM1S(10).and.NMSTEP.gt.NJM1S(9)) then
         j=10 
       endif  ! of month 10
      if (NMSTEP.le.NJM1S(11).and.NMSTEP.gt.NJM1S(10)) then
         j=11
       endif  ! of month 11
      if (NMSTEP.le.NJM1S(12).and.NMSTEP.gt.NJM1S(11)) then
         j=12
       endif  ! of month 12

       if(j.lt.12) then
          q21m(j)=q21m(j)+q21
          q62m(j)=q62m(j)+q62
          q23m(j)=q23m(j)+q23
          q73m(j)=q73m(j)+q73
          q87m(j)=q87m(j)+q87
          q34m(j)=q34m(j)+q34
          q76m(j)=q76m(j)+q76
          q15m(j)=q15m(j)+q15
          q56m(j)=q56m(j)+q56
          q48m(j)=q48m(j)+q48

          q21sm(j)=q21sm(j)+q21s
          q23sm(j)=q23sm(j)+q23s
          q34sm(j)=q34sm(j)+q34s 

          q21tm(j)=q21tm(j)+q21t
          q23tm(j)=q23tm(j)+q23t
          q34tm(j)=q34tm(j)+q34t

          Savm(j)=Savm(j)+Sav
          Tavm(j)=Tavm(j)+Tav
          tetaavm(j)=tetaavm(j)+tetaav

          co2am(j)=co2am(j)+co2a
          CTavm(j)=CTavm(j)+CTav
          ATavm(j)=ATavm(j)+ATav
          PO4avm(j)=PO4avm(j)+PO4av
          sPRODm(j)=sPRODm(j)+sPROD
          hPRODm(j)=hPRODm(j)+hPROD
          totPRODm(j)=totPRODm(j)+totPROD

          Imm(1,j)=Imm(1,j)+Im(1)
          Iflowm(1,j)=Iflowm(1,j)+Iflow(1)
          ILnewm(1,j)=ILnewm(1,j)+InewTemp(1)
          IceHsm(j)=IceHsm(j)+IceHs
          IvolSm(j)=IvolSm(j)+IcevolS(k)
          Imm(4,j)=Imm(4,j)+Im(4)
          Iflowm(4,j)=Iflowm(4,j)+Iflow(4)
          ILnewm(4,j)=ILnewm(4,j)+InewTemp(4)
          IceHNm(j)=IceHNm(j)+IceHN
          IvolNm(j)=IvolNm(j)+IcevolN(k)

          FT21m(j)=FT21m(j)+FT21
          FT23m(j)=FT23m(j)+FT23
          FT34m(j)=FT34m(j)+FT34
          FW21m(j)=FW21m(j)+FW21
          FW23m(j)=FW23m(j)+FW23
          FW34m(j)=FW34m(j)+FW34

          FW1m(j)=FW1m(j)+FW1
          FW4m(j)=FW4m(j)+FW4

c hezi 16.9.1999 add the following
          tf21m(j)=tf21m(j)+tf21
          tf15m(j)=tf15m(j)+tf15
          tdf12m(j)=tdf12m(j)+tdf12
          tdf15m(j)=tdf15m(j)+tdf15

          ctf21m(j)=ctf21m(j)+ctf21
          ctf15m(j)=ctf15m(j)+ctf15
          ctf56m(j)=ctf56m(j)+ctf56
          ctf34m(j)=ctf34m(j)+ctf34
          ctf48m(j)=ctf48m(j)+ctf48
          ctf87m(j)=ctf87m(j)+ctf87
          ctf62m(j)=ctf62m(j)+ctf62
          ctf73m(j)=ctf73m(j)+ctf73

          ctdf15m(j)=ctdf15m(j)+ctdf15
          ctdf26m(j)=ctdf26m(j)+ctdf26
          ctdf48m(j)=ctdf48m(j)+ctdf48
          ctdf12m(j)=ctdf12m(j)+ctdf12
          ctdf34m(j)=ctdf34m(j)+ctdf34
          ctdf56m(j)=ctdf56m(j)+ctdf56
          ctdf78m(j)=ctdf78m(j)+ctdf78

          po4f21m(j)=po4f21m(j)+po4f21
          po4f15m(j)=po4f15m(j)+po4f15
          po4f56m(j)=po4f56m(j)+po4f56
          po4f34m(j)=po4f34m(j)+po4f34
          po4f48m(j)=po4f48m(j)+po4f48
          po4f87m(j)=po4f87m(j)+po4f87
          po4f62m(j)=po4f62m(j)+po4f62
          po4f73m(j)=po4f73m(j)+po4f73

          po4df15m(j)=po4df15m(j)+po4df15
          po4df26m(j)=po4df26m(j)+po4df26
          po4df48m(j)=po4df48m(j)+po4df48
          po4df12m(j)=po4df12m(j)+po4df12
          po4df34m(j)=po4df34m(j)+po4df34
          po4df56m(j)=po4df56m(j)+po4df56
          po4df78m(j)=po4df78m(j)+po4df78
          

          vmixm(1,j)=vmixm(1,j)+(tdf15/( t(1,k)-t(5,k) ))
          vmixm(2,j)=vmixm(2,j)+(tdf26/( t(2,k)-t(6,k) ))
          vmixm(3,j)=vmixm(3,j)+(tdf37/( t(3,k)-t(7,k) ))
          vmixm(4,j)=vmixm(4,j)+(tdf48/( t(4,k)-t(8,k) ))
          
          do i=1,8
             tm(i,j)=tm(i,j)+t(i,k)
             sm(i,j)=sm(i,j)+s(i,k)
             tcm(i,j)=tcm(i,j)+tc(i)
             scm(i,j)=scm(i,j)+sc(i)
             rhom(i,j)=rhom(i,j)+rho(i)
             rhosm(i,j)=rhosm(i,j)+rhos(i)
             rhotm(i,j)=rhotm(i,j)+rhot(i)

             CTm(i,j)=CTm(i,j)+CT(i,k)
             ATm(i,j)=ATm(i,j)+AT(i,k)
             HCO3m(i,j)=HCO3m(i,j)+HCO3(i)
             CO3m(i,j)=CO3m(i,j)+CO3(i)
             CO2m(i,j)=CO2m(i,j)+CO2(i)
             CO2aam(i,j)=CO2aam(i,j)+CO2(i)/K0(i)
             am(i,j)=am(i,j)+a(i)
             PO4m(i,j)=PO4m(i,j)+PO4(i,k)
             PRODm(i,j)=PRODm(i,j)+PROD(i)
             RAINm(i,j)=RAINm(i,j)+RAIN(i)

          enddo

          do i=1,4
             tetam(i,j)=tetam(i,j)+teta(i,k)
             Haom(i,j)=Haom(i,j)+Hao(i)
             Hlhm(i,j)=Hlhm(i,j)+Hlh(i)
             Houtm(i,j)=Houtm(i,j)+Hout(i)
             Hinm(i,j)=Hinm(i,j)+Hin(i)
             tfluxm(i,j)=tfluxm(i,j)+tflux(i)
             sfluxm(i,j)=sfluxm(i,j)+sflux(i)
             FWAOm(i,j)=FWAOm(i,j)+FWAO(i)
             qm(i,j)=qm(i,j)+q(i)
             albedom(i,j)=albedom(i,j)+albedo(i)
             solarm(i,j)=solarm(i,j)+solar(i)
             Havailm(i,j)=Havailm(i,j)+Havail(i)
             Hneedm(i,j)=Hneedm(i,j)+Hneed(i)
             sivolm(i,j)=sivolm(i,j)+sivol(i,k)

             CTfluxm(i,j)=CTfluxm(i,j)+CTflux(i)
             RedCm(i,j)=RedCm(i,j)+RedC(i)
          enddo

          do i=1,4,3
             fsim(i,j)=fsim(i,j)+fsi(i,k)
             Imeltm(i,j)=Imeltm(i,j)+Imelt(i)
             LSIdiffm(i,j)=LSIdiffm(i,j)+LSIdiff(i)
             Inewfm(i,j)=Inewfm(i,j)+Inewf(i)
             Idepthm(i,j)=Idepthm(i,j)+Idepth(i)
             Ipnewm(i,j)=Ipnewm(i,j)+Ipnew(i)
             fim(i,j)=fim(i,j)+fi(i,k)
             Inewm(i,j)=Inewm(i,j)+inew(i)
             LIdiffm(i,j)=LIdiffm(i,j)+LIdiff(i)
          enddo

            goto 100

          endif

       if(j.eq.12) then
          q21m(j)=q21m(j)+q21
          q62m(j)=q62m(j)+q62
          q23m(j)=q23m(j)+q23
          q73m(j)=q73m(j)+q73
          q87m(j)=q87m(j)+q87
          q34m(j)=q34m(j)+q34
          q76m(j)=q76m(j)+q76
          q15m(j)=q15m(j)+q15
          q56m(j)=q56m(j)+q56
          q48m(j)=q48m(j)+q48

          q21sm(j)=q21sm(j)+q21s
          q23sm(j)=q23sm(j)+q23s
          q34sm(j)=q34sm(j)+q34s 

          q21tm(j)=q21tm(j)+q21t
          q23tm(j)=q23tm(j)+q23t
          q34tm(j)=q34tm(j)+q34t

          Savm(j)=Savm(j)+Sav
          Tavm(j)=Tavm(j)+Tav
          tetaavm(j)=tetaavm(j)+tetaav

          co2am(j)=co2am(j)+co2a
          CTavm(j)=CTavm(j)+CTav
          ATavm(j)=ATavm(j)+ATav
          PO4avm(j)=PO4avm(j)+PO4av
          sPRODm(j)=sPRODm(j)+sPROD
          hPRODm(j)=hPRODm(j)+hPROD
          totPRODm(j)=totPRODm(j)+totPROD

          Imm(1,j)=Imm(1,j)+Im(1)
          Iflowm(1,j)=Iflowm(1,j)+Iflow(1)
          ILnewm(1,j)=ILnewm(1,j)+InewTemp(1)
          IceHsm(j)=IceHsm(j)+IceHs
          IvolSm(j)=IvolSm(j)+IcevolS(k)
          Imm(4,j)=Imm(4,j)+Im(4)
          Iflowm(4,j)=Iflowm(4,j)+Iflow(4)
          ILnewm(4,j)=ILnewm(4,j)+InewTemp(4)
          IceHNm(j)=IceHNm(j)+IceHN
          IvolNm(j)=IvolNm(j)+IcevolN(k) 

          FT21m(j)=FT21m(j)+FT21
          FT23m(j)=FT23m(j)+FT23
          FT34m(j)=FT34m(j)+FT34
          FW21m(j)=FW21m(j)+FW21
          FW23m(j)=FW23m(j)+FW23
          FW34m(j)=FW34m(j)+FW34

          FW1m(j)=FW1m(j)+FW1
          FW4m(j)=FW4m(j)+FW4

          tf21m(j)=tf21m(j)+tf21
          tf15m(j)=tf15m(j)+tf15
          tdf12m(j)=tdf12m(j)+tdf12
          tdf15m(j)=tdf15m(j)+tdf15


          ctf21m(j)=ctf21m(j)+ctf21
          ctf15m(j)=ctf15m(j)+ctf15
          ctf56m(j)=ctf56m(j)+ctf56
          ctf34m(j)=ctf34m(j)+ctf34
          ctf48m(j)=ctf48m(j)+ctf48
          ctf87m(j)=ctf87m(j)+ctf87
          ctf62m(j)=ctf62m(j)+ctf62
          ctf73m(j)=ctf73m(j)+ctf73

          ctdf15m(j)=ctdf15m(j)+ctdf15
          ctdf26m(j)=ctdf26m(j)+ctdf26
          ctdf48m(j)=ctdf48m(j)+ctdf48
          ctdf12m(j)=ctdf12m(j)+ctdf12
          ctdf34m(j)=ctdf34m(j)+ctdf34
          ctdf56m(j)=ctdf56m(j)+ctdf56
          ctdf78m(j)=ctdf78m(j)+ctdf78

          po4f21m(j)=po4f21m(j)+po4f21
          po4f15m(j)=po4f15m(j)+po4f15
          po4f56m(j)=po4f56m(j)+po4f56
          po4f34m(j)=po4f34m(j)+po4f34
          po4f48m(j)=po4f48m(j)+po4f48
          po4f87m(j)=po4f87m(j)+po4f87
          po4f62m(j)=po4f62m(j)+po4f62
          po4f73m(j)=po4f73m(j)+po4f73

          po4df15m(j)=po4df15m(j)+po4df15
          po4df26m(j)=po4df26m(j)+po4df26
          po4df48m(j)=po4df48m(j)+po4df48
          po4df12m(j)=po4df12m(j)+po4df12
          po4df34m(j)=po4df34m(j)+po4df34
          po4df56m(j)=po4df56m(j)+po4df56
          po4df78m(j)=po4df78m(j)+po4df78

          vmixm(1,j)=vmixm(1,j)+(tdf15/( t(1,k)-t(5,k) ))
          vmixm(2,j)=vmixm(2,j)+(tdf26/( t(2,k)-t(6,k) ))
          vmixm(3,j)=vmixm(3,j)+(tdf37/( t(3,k)-t(7,k) ))
          vmixm(4,j)=vmixm(4,j)+(tdf48/( t(4,k)-t(8,k) ))          

          do i=1,8
             tm(i,j)=tm(i,j)+t(i,k)
             sm(i,j)=sm(i,j)+s(i,k)
             tcm(i,j)=tcm(i,j)+tc(i)
             scm(i,j)=scm(i,j)+sc(i)
             rhom(i,j)=rhom(i,j)+rho(i)
             rhosm(i,j)=rhosm(i,j)+rhos(i)
             rhotm(i,j)=rhotm(i,j)+rhot(i)

             CTm(i,j)=CTm(i,j)+CT(i,k)
             ATm(i,j)=ATm(i,j)+AT(i,k)
             HCO3m(i,j)=HCO3m(i,j)+HCO3(i)
             CO3m(i,j)=CO3m(i,j)+CO3(i)
             CO2m(i,j)=CO2m(i,j)+CO2(i)
             CO2aam(i,j)=CO2aam(i,j)+CO2(i)/K0(i)
             am(i,j)=am(i,j)+a(i)
             PO4m(i,j)=PO4m(i,j)+PO4(i,k)
             PRODm(i,j)=PRODm(i,j)+PROD(i)
             RAINm(i,j)=RAINm(i,j)+RAIN(i)

          enddo

          do i=1,4
             tetam(i,j)=tetam(i,j)+teta(i,k)
             Haom(i,j)=Haom(i,j)+Hao(i)
             Hlhm(i,j)=Hlhm(i,j)+Hlh(i)
             Houtm(i,j)=Houtm(i,j)+Hout(i)
             Hinm(i,j)=Hinm(i,j)+Hin(i)
             tfluxm(i,j)=tfluxm(i,j)+tflux(i)
             sfluxm(i,j)=sfluxm(i,j)+sflux(i)
             FWAOm(i,j)=FWAOm(i,j)+FWAO(i)
             qm(i,j)=qm(i,j)+q(i)
             albedom(i,j)=albedom(i,j)+albedo(i)
             solarm(i,j)=solarm(i,j)+solar(i)
             Havailm(i,j)=Havailm(i,j)+Havail(i)
             Hneedm(i,j)=Hneedm(i,j)+Hneed(i)
             sivolm(i,j)=sivolm(i,j)+sivol(i,k)

             CTfluxm(i,j)=CTfluxm(i,j)+CTflux(i)
             RedCm(i,j)=RedCm(i,j)+RedC(i)
          enddo

          do i=1,4,3
             fsim(i,j)=fsim(i,j)+fsi(i,k)
             Imeltm(i,j)=Imeltm(i,j)+Imelt(i)
             LSIdiffm(i,j)=LSIdiffm(i,j)+LSIdiff(i)
             Inewfm(i,j)=Inewfm(i,j)+Inewf(i)
             Idepthm(i,j)=Idepthm(i,j)+Idepth(i)
             Ipnewm(i,j)=Ipnewm(i,j)+Ipnew(i)
             fim(i,j)=fim(i,j)+fi(i,k)
             Inewm(i,j)=Inewm(i,j)+inew(i)
             LIdiffm(i,j)=LIdiffm(i,j)+LIdiff(i)
          enddo


      if (NMSTEP.eq.NJM1S(12)) then
       do m=1,12

          q21m(m)=q21m(m)/NJMS(m)
          q62m(m)=q62m(m)/NJMS(m)
          q23m(m)=q23m(m)/NJMS(m) 
          q73m(m)=q73m(m)/NJMS(m) 
          q87m(m)=q87m(m)/NJMS(m) 
          q34m(m)=q34m(m)/NJMS(m) 
          q76m(m)=q76m(m)/NJMS(m) 
          q15m(m)=q15m(m)/NJMS(m) 
          q56m(m)=q56m(m)/NJMS(m) 
          q48m(m)=q48m(m)/NJMS(m)

          q21sm(m)=q21sm(m)/NJMS(m) 
          q23sm(m)=q23sm(m)/NJMS(m) 
          q34sm(m)=q34sm(m)/NJMS(m) 

          q21tm(m)=q21tm(m)/NJMS(m) 
          q23tm(m)=q23tm(m)/NJMS(m) 
          q34tm(m)=q34tm(m)/NJMS(m)

          Savm(m)=Savm(m)/NJMS(m) 
          Tavm(m)=Tavm(m)/NJMS(m) 
          tetaavm(m)=tetaavm(m)/NJMS(m)

          co2am(m)=co2am(m)/NJMS(m) 
          CTavm(m)=CTavm(m)/NJMS(m) 
          ATavm(m)=ATavm(m)/NJMS(m) 
          PO4avm(m)=PO4avm(m)/NJMS(m) 
          sPRODm(m)=sPRODm(m)/NJMS(m)
          hPRODm(m)=hPRODm(m)/NJMS(m)
          totPRODm(m)=totPRODm(m)/NJMS(m) 

          Imm(1,m)=Imm(1,m)/NJMS(m)
          Iflowm(1,m)=Iflowm(1,m)/NJMS(m)
          ILnewm(1,m)=ILnewm(1,m)/NJMS(m) 
          IceHsm(m)=IceHsm(m)/NJMS(m) 
          IvolSm(m)=IvolSm(m)/NJMS(m) 
          Imm(4,m)=Imm(4,m)/NJMS(m) 
          Iflowm(4,m)=Iflowm(4,m)/NJMS(m) 
          ILnewm(4,m)=ILnewm(4,m)/NJMS(m) 
          IceHNm(m)=IceHNm(m)/NJMS(m) 
          IvolNm(m)=IvolNm(m)/NJMS(m) 

          FT21m(m)=FT21m(m)/NJMS(m) 
          FT23m(m)=FT23m(m)/NJMS(m) 
          FT34m(m)=FT34m(m)/NJMS(m) 
          FW21m(m)=FW21m(m)/NJMS(m) 
          FW23m(m)=FW23m(m)/NJMS(m) 
          FW34m(m)=FW34m(m)/NJMS(m) 

          FW1m(m)=FW1m(m)/NJMS(m)
          FW4m(m)=FW4m(m)/NJMS(m)

          tf21m(m)=tf21m(m)/NJMS(m)
          tf15m(m)=tf15m(m)/NJMS(m)
          tdf12m(m)=tdf12m(m)/NJMS(m)
          tdf15m(m)=tdf15m(m)/NJMS(m)


          ctf21m(m)=ctf21m(m)/NJMS(m)
          ctf15m(m)=ctf15m(m)/NJMS(m)
          ctf56m(m)=ctf56m(m)/NJMS(m)
          ctf34m(m)=ctf34m(m)/NJMS(m)
          ctf48m(m)=ctf48m(m)/NJMS(m)
          ctf87m(m)=ctf87m(m)/NJMS(m)
          ctf62m(m)=ctf62m(m)/NJMS(m)
          ctf73m(m)=ctf73m(m)/NJMS(m)

          ctdf15m(m)=ctdf15m(m)/NJMS(m)
          ctdf26m(m)=ctdf26m(m)/NJMS(m)
          ctdf48m(m)=ctdf48m(m)/NJMS(m)
          ctdf12m(m)=ctdf12m(m)/NJMS(m)
          ctdf34m(m)=ctdf34m(m)/NJMS(m)
          ctdf56m(m)=ctdf56m(m)/NJMS(m)
          ctdf78m(m)=ctdf78m(m)/NJMS(m)

          po4f21m(m)=po4f21m(m)/NJMS(m)
          po4f15m(m)=po4f15m(m)/NJMS(m)
          po4f56m(m)=po4f56m(m)/NJMS(m)
          po4f34m(m)=po4f34m(m)/NJMS(m)
          po4f48m(m)=po4f48m(m)/NJMS(m)
          po4f87m(m)=po4f87m(m)/NJMS(m)
          po4f62m(m)=po4f62m(m)/NJMS(m)
          po4f73m(m)=po4f73m(m)/NJMS(m)

          po4df15m(m)=po4df15m(m)/NJMS(m)
          po4df26m(m)=po4df26m(m)/NJMS(m)
          po4df48m(m)=po4df48m(m)/NJMS(m)
          po4df12m(m)=po4df12m(m)/NJMS(m)
          po4df34m(m)=po4df34m(m)/NJMS(m)
          po4df56m(m)=po4df56m(m)/NJMS(m)
          po4df78m(m)=po4df78m(m)/NJMS(m)        

          vmixm(1,m)=vmixm(1,m)/NJMS(m)
          vmixm(2,m)=vmixm(2,m)/NJMS(m)
          vmixm(3,m)=vmixm(3,m)/NJMS(m)
          vmixm(4,m)=vmixm(4,m)/NJMS(m)

          do i=1,8
             tm(i,m)=tm(i,m)/NJMS(m) 
             sm(i,m)=sm(i,m)/NJMS(m)
             tcm(i,m)=tcm(i,m)/NJMS(m) 
             scm(i,m)=scm(i,m)/NJMS(m)
             rhom(i,m)=rhom(i,m)/NJMS(m) 
             rhosm(i,m)=rhosm(i,m)/NJMS(m)
             rhotm(i,m)=rhotm(i,m)/NJMS(m)

             CTm(i,m)=CTm(i,m)/NJMS(m)
             ATm(i,m)=ATm(i,m)/NJMS(m)
             HCO3m(i,m)=HCO3m(i,m)/NJMS(m)
             CO3m(i,m)=CO3m(i,m)/NJMS(m)
             CO2m(i,m)=CO2m(i,m)/NJMS(m)
             CO2aam(i,m)=CO2aam(i,m)/NJMS(m)
             am(i,m)=am(i,m)/NJMS(m)
             PO4m(i,m)=PO4m(i,m)/NJMS(m)
             PRODm(i,m)=PRODm(i,m)/NJMS(m)


          enddo

          do i=1,4
             tetam(i,m)=tetam(i,m)/NJMS(m) 
             Haom(i,m)=Haom(i,m)/NJMS(m) 
             Hlhm(i,m)=Hlhm(i,m)/NJMS(m)
             Houtm(i,m)=Houtm(i,m)/NJMS(m) 
             Hinm(i,m)=Hinm(i,m)/NJMS(m) 
             tfluxm(i,m)=tfluxm(i,m)/NJMS(m) 
             sfluxm(i,m)=sfluxm(i,m)/NJMS(m) 
             FWAOm(i,m)=FWAOm(i,m)/NJMS(m) 
             qm(i,m)=qm(i,m)/NJMS(m) 
             albedom(i,m)=albedom(i,m)/NJMS(m) 
             solarm(i,m)=solarm(i,m)/NJMS(m)
             Havailm(i,m)=Havailm(i,m)/NJMS(m)
             Hneedm(i,m)=Hneedm(i,m)/NJMS(m)
             sivolm(i,m)=sivolm(i,m)/NJMS(m)

c             CTfluxm(i,m)=CTfluxm(i,m)/(NJMS(m)*Oarea(i))
             CTfluxm(i,m)=CTfluxm(i,m)/NJMS(m)
             RedCm(i,m)=RedCm(i,m)/NJMS(m)
             RAINm(i,m)=RAINm(i,m)/(NJMS(m)*RedC(i))
          enddo

          do i=1,4,3
             fsim(i,m)=fsim(i,m)/NJMS(m) 
             Imeltm(i,m)=Imeltm(i,m)/NJMS(m)
             LSIdiffm(i,m)=LSIdiffm(i,m)/NJMS(m)
             Inewfm(i,m)=Inewfm(i,m)/NJMS(m) 
             Idepthm(i,m)=Idepthm(i,m)/NJMS(m) 
             Ipnewm(i,m)=Ipnewm(i,m)/NJMS(m) 
             fim(i,m)=fim(i,m)/NJMS(m)
             Inewm(i,m)=Inewm(i,m)/NJMS(m)
             LIdiffm(i,m)=LIdiffm(i,m)/NJMS(m)
          enddo
          inewn2=0.0
          do i=1,12
             inewn2=inewn2+ILnewm(4,i)
          enddo
          inewn2=inewn2/12.
      
       enddo
       ndat=ndat+1
       call output
       MFLAG=.false.
       NPSTP=0
       NMSTEP=0

c       ndat=ndat+1
c       call output
c       MFLAG=.false.
c       NPSTP=0
c       NMSTEP=0
c       do m=1,12
c       print *,'m=',m,'hlh4=',hlhm(4,m),'fw34=',fw34/1.d+6
c       enddo
c hezi 23.8.98
c       print *,'kday=',kday
       endif

c          goto 100
             
       endif  ! of month 12








 100   continue
       return
       end
