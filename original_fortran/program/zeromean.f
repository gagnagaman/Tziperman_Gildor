      subroutine zeromean
      implicit none
      include 'averaging.h'
      integer i,m

       do m=1,12

          q21m(m)=0.
          q62m(m)=0.
          q23m(m)=0.
          q73m(m)=0.
          q87m(m)=0.
          q34m(m)=0.
          q76m(m)=0.
          q15m(m)=0.
          q56m(m)=0.
          q48m(m)=0.

          q21sm(m)=0.
          q23sm(m)=0.
          q34sm(m)=0.

          q21tm(m)=0.
          q23tm(m)=0.
          q34tm(m)=0.

          Savm(m)=0.
          Tavm(m)=0.
          tetaavm(m)=0.


          co2am(m)=0.
          CTavm(m)=0.
          ATavm(m)=0.
          PO4avm(m)=0.
          sPRODm(m)=0.
          hPRODm(m)=0.
          totPRODm(m)=0.

          Imm(1,m)=0.
          Iflowm(1,m)=0.
          ILnewm(1,m)=0.
          IceHsm(m)=0.
          IvolSm(m)=0.
          Imm(4,m)=0.
          Iflowm(4,m)=0.
          ILnewm(4,m)=0.
          IceHNm(m)=0.
          IvolNm(m)=0.

          FT21m(m)=0.
          FT23m(m)=0. 
          FT34m(m)=0. 
          FW21m(m)=0.
          FW23m(m)=0.
          FW34m(m)=0. 
*dbgz:
          FW4m(m)=0.
          FW1m(m)=0.

          tf21m(m)=0.
          tf15m(m)=0.
          tdf15m(m)=0.
          tdf12m(m)=0.

          ctf21m(m)=0.
          ctf15m(m)=0.
          ctf56m(m)=0.
          ctf34m(m)=0.
          ctf48m(m)=0.
          ctf87m(m)=0.
          ctf62m(m)=0.
          ctf73m(m)=0.

          ctdf15m(m)=0.
          ctdf26m(m)=0.
          ctdf48m(m)=0.
          ctdf12m(m)=0.
          ctdf34m(m)=0.
          ctdf56m(m)=0.
          ctdf78m(m)=0.

          po4f21m(m)=0.
          po4f15m(m)=0.
          po4f56m(m)=0.
          po4f34m(m)=0.
          po4f48m(m)=0.
          po4f87m(m)=0.
          po4f62m(m)=0.
          po4f73m(m)=0.

          po4df15m(m)=0.
          po4df26m(m)=0.
          po4df48m(m)=0.
          po4df12m(m)=0.
          po4df34m(m)=0.
          po4df56m(m)=0.
          po4df78m(m)=0.

          do i=1,8
             tm(i,m)=0.
             sm(i,m)=0.
             rhom(i,m)=0.
             CTm(i,m)=0.
             ATm(i,m)=0.
             HCO3m(i,m)=0.
             CO3m(i,m)=0.
             CO2m(i,m)=0.
             CO2aam(i,m)=0.
             am(i,m)=0.
             PO4m(i,m)=0.
             PRODm(i,m)=0.
             RAINm(i,m)=0.
          enddo

          do i=1,4
             tetam(i,m)=0. 
             Haom(i,m)=0.
             Hlhm(i,m)=0.
             Houtm(i,m)=0.
             Hinm(i,m)=0. 
             tfluxm(i,m)=0. 
             sfluxm(i,m)=0. 
             FWAOm(i,m)=0.
             qm(i,m)=0.
             albedom(i,m)=0.
             solarm(i,m)=0.
             Havailm(i,m)=0.
             Hneedm(i,m)=0.
             sivolm(i,m)=0.

             CTfluxm(i,m)=0.
             vmixm(i,m)=0.
             RedCm(i,m)=0.
          enddo

          do i=1,4,3
             fsim(i,m)=0. 
             Imeltm(i,m)=0. 
             Inewfm(i,m)=0. 
             Idepthm(i,m)=0.
             Ipnewm(i,m)=0. 
             fim(i,m)=0.
             Inewm(i,m)=0.
             lidiffm(i,m)=0.
             lsidiffm(i,m)=0.
          enddo
      
       enddo

       return
       end
