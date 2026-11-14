      subroutine output

c store the snapshot data

      implicit none
      include 'couple.h'
      include 'atmosphere.h'
      include 'ocean.h'
      include 'seaice.h'
      include 'landice.h'
      include 'biochemical.h'
      include 'averaging.h'
      real*8 aa_sm

        do i=1,4
           do j=1,12
              write(12,rec=irec) sngl(tetam(i,j)-273.)
              irec=irec+1
           enddo
        enddo
c hezi: 16.9.1999 change from Haom(i,j)*Aarea(i) to Haom(i,j)*Oarea(i)
        do i=1,4
           do j=1,12
              write(12,rec=irec) sngl(Haom(i,j)*Oarea(i)/1.d+15)
              irec=irec+1
           enddo
        enddo

        do i=1,4
           do j=1,12
              write(12,rec=irec) sngl(Houtm(i,j)*Aarea(i)/1.d+15)
              irec=irec+1
           enddo
        enddo

        do i=1,4
           do j=1,12
*dbgz:
              write(12,rec=irec) sngl(Hinm(i,j)*Aarea(i)/(ac1*1.d15))
              irec=irec+1
           enddo
        enddo

        do i=1,4
           do j=1,12
              write(12,rec=irec) sngl(tfluxm(i,j)*CwrhoW)
              irec=irec+1
          enddo
        enddo

        do i=1,4
           do j=1,12
              write(12,rec=irec) sngl(-1.*sfluxm(i,j)*3.15d+10/Sref)
              irec=irec+1
           enddo
        enddo

        do i=1,4
           do j=1,12
              write(12,rec=irec) sngl(FWAOm(i,j)/1.d+6)
              irec=irec+1
           enddo
        enddo

           do j=1,12
              write(12,rec=irec) sngl(FT21m(j)*Aarea(1)/1.d+15/ac1)
              irec=irec+1
           enddo

           do j=1,12
              write(12,rec=irec) sngl(FT23m(j)*Aarea(3)/1.d+15/ac1)
              irec=irec+1
           enddo

           do j=1,12
              write(12,rec=irec) sngl(FT34m(j)*Aarea(4)/1.d+15/ac1)
              irec=irec+1
           enddo

           do j=1,12
              write(12,rec=irec) sngl(FW21m(j)/1.d+6)
              irec=irec+1
           enddo

           do j=1,12
              write(12,rec=irec) sngl(FW23m(j)/1.d+6)
              irec=irec+1
           enddo

           do j=1,12
              write(12,rec=irec) sngl(FW34m(j)/1.d+6)
              irec=irec+1
           enddo

           do j=1,12
              write(12,rec=irec) sngl(FW1m(j)/1.d+6)
              irec=irec+1
           enddo

           do j=1,12
              write(12,rec=irec) sngl(FW4m(j)/1.d+6)
              irec=irec+1
           enddo

        do i=1,8
           do j=1,12
              write(12,rec=irec) sngl(tm(i,j)-273.)
              irec=irec+1
           enddo
        enddo

        do i=1,8
           do j=1,12
              write(12,rec=irec) sngl(sm(i,j)-1.d-20)
              irec=irec+1
           enddo
        enddo

c        do i=1,8
c           do j=1,12
c              write(12,rec=irec) sngl(rhom(i,j)-1000.0)
c              irec=irec+1
c           enddo
c        enddo

        do i=1,8
           do j=1,12
              write(12,rec=irec) sngl(rhom(i,j)-1.d-20)
              irec=irec+1
           enddo
        enddo

        do i=1,8
           do j=1,12
              write(12,rec=irec) sngl(rhosm(i,j)-1000.0)
              irec=irec+1
           enddo
        enddo

        do i=1,8
           do j=1,12
              write(12,rec=irec) sngl(rhotm(i,j)-1000.0)
              irec=irec+1
           enddo
        enddo

           do j=1,12
              write(12,rec=irec) sngl(q21m(j)/1.d+6)
              irec=irec+1
           enddo

           do j=1,12
              write(12,rec=irec) sngl(q34m(j)/1.d+6)
              irec=irec+1
           enddo

           do j=1,12
              write(12,rec=irec) sngl(q23m(j)/1.d+6)
              irec=irec+1
           enddo

           do j=1,12
              write(12,rec=irec) sngl(q62m(j)/1.d+6)
              irec=irec+1
           enddo

           do j=1,12
              write(12,rec=irec) sngl(q73m(j)/1.d+6)
              irec=irec+1
           enddo

           do j=1,12
              write(12,rec=irec) sngl(q76m(j)/1.d+6)
              irec=irec+1
           enddo

           do j=1,12
              write(12,rec=irec) sngl(q56m(j)/1.d+6)
              irec=irec+1
           enddo

           do j=1,12
              write(12,rec=irec) sngl(q87m(j)/1.d+6)
              irec=irec+1
           enddo


 
           do j=1,12
              write(12,rec=irec) sngl(fsim(1,j)-1.d-20)
              irec=irec+1
           enddo


           do j=1,12
              write(12,rec=irec) sngl(fsim(4,j)-1.d-20)
              irec=irec+1
           enddo

           do j=1,12
              write(12,rec=irec) sngl(fim(1,j)-1.d-20)
              irec=irec+1
           enddo

           do j=1,12
              write(12,rec=irec) sngl(fim(4,j)-1.d-20)
              irec=irec+1
           enddo

        do i=1,4
           do j=1,12
              write(12,rec=irec) sngl(qm(i,j)-1.d-20)
              irec=irec+1
           enddo
        enddo
 
           do j=1,12
              write(12,rec=irec) sngl(Tavm(j)-1.d-20)
              irec=irec+1
           enddo

           do j=1,12
              write(12,rec=irec) sngl(tetaavm(j)-1.d-20)
              irec=irec+1
           enddo

           do j=1,12
              write(12,rec=irec) sngl(Savm(j)-1.d-20)
              irec=irec+1
           enddo

        do j=1,12
           write(12,rec=irec) sngl((Imm(1,j))*(rhoLice/rhoW)/1.d+6)
           irec=irec+1
        enddo

        do j=1,12
           write(12,rec=irec) sngl((Imm(4,j))*(rhoLice/rhoW)/1.d+6)
           irec=irec+1
        enddo

        do j=1,12
         write(12,rec=irec) sngl((Iflowm(1,j)/tdt)*(rhoLice/rhoW)/1.d+6)
         irec=irec+1
        enddo

        do j=1,12
         write(12,rec=irec) sngl((Iflowm(4,j)/tdt)*(rhoLice/rhoW)/1.d+6)
         irec=irec+1
        enddo

        do j=1,12
         write(12,rec=irec) sngl((ILnewm(1,j))*(rhoLice/rhoW)/1.d+6)
          irec=irec+1
        enddo

        do j=1,12
         write(12,rec=irec) sngl((ILnewm(4,j))*(rhoLice/rhoW)/1.d+6)
          irec=irec+1
        enddo

           do j=1,12
              write(12,rec=irec) sngl(IceHsm(j)-1.d-20)
              irec=irec+1
           enddo

           do j=1,12
              write(12,rec=irec) sngl(IceHnm(j)-1.d-20)
              irec=irec+1
           enddo

           do j=1,12
              write(12,rec=irec) sngl(IvolSm(j)*(rhoLice/rhoW))
              irec=irec+1
           enddo

           do j=1,12
              write(12,rec=irec) sngl(IVolNm(j)*(rhoLice/rhoW))
              irec=irec+1
           enddo

        do i=1,4
           do j=1,12
              write(12,rec=irec) sngl(albedom(i,j)-1.d-20)
              irec=irec+1
           enddo
        enddo

        do i=1,4,3
           do j=1,12
              write(12,rec=irec) sngl(Idepthm(i,j)-1.d-20)
              irec=irec+1 
           enddo
        enddo


        do j=1,12
          write(12,rec=irec) sngl((Imeltm(1,j)/tdt)*(rhoIce/rhoW)/1.d+6)
          irec=irec+1
        enddo

        do j=1,12
         write(12,rec=irec) sngl((Imeltm(4,j)/tdt)*(rhoLice/rhoW)/1.d+6)
         irec=irec+1
        enddo

        do j=1,12
          write(12,rec=irec) sngl((Inewfm(1,j)/tdt)*(rhoIce/rhoW)/1.d+6)
          irec=irec+1
        enddo

        do j=1,12
          write(12,rec=irec) sngl((Inewfm(4,j)/tdt)*(rhoIce/rhoW)/1.d+6)
          irec=irec+1
        enddo

        do j=1,12
          write(12,rec=irec) sngl((Ipnewm(1,j)/tdt)*(rhoIce/rhoW)/1.d+6)
          irec=irec+1
        enddo

        do j=1,12
          write(12,rec=irec) sngl((Ipnewm(4,j)/tdt)*(rhoIce/rhoW)/1.d+6)
          irec=irec+1
        enddo
        
 

        do i=1,4
           do j=1,12
              write(12,rec=irec) sngl(solarm(i,j)-1.d-20)
              irec=irec+1
           enddo
        enddo

           do j=1,12
              write(12,rec=irec) sngl(q21sm(j)/1.d+6)
              irec=irec+1
           enddo

           do j=1,12
              write(12,rec=irec) sngl(q21tm(j)/1.d+6)
              irec=irec+1
           enddo

           do j=1,12
              write(12,rec=irec) sngl(q23sm(j)/1.d+6)
              irec=irec+1
           enddo

           do j=1,12
              write(12,rec=irec) sngl(q23tm(j)/1.d+6)
              irec=irec+1
           enddo

           do j=1,12
              write(12,rec=irec) sngl(q34sm(j)/1.d+6)
              irec=irec+1
           enddo

           do j=1,12
              write(12,rec=irec) sngl(q34tm(j)/1.d+6)
              irec=irec+1
           enddo

           do j=1,12
              write(12,rec=irec) sngl(Havailm(1,j)-1.d-20)
              irec=irec+1
           enddo

           do j=1,12
              write(12,rec=irec) sngl(Havailm(4,j)-1.d-20)
              irec=irec+1
           enddo

           do j=1,12
              write(12,rec=irec) sngl(Hneedm(1,j)-1.d-20)
              irec=irec+1
           enddo

           do j=1,12
              write(12,rec=irec) sngl(Hneedm(4,j)-1.d-20)
              irec=irec+1
           enddo

           do j=1,12
              write(12,rec=irec) sngl(sivolm(1,j)-1.d-20)
              irec=irec+1
           enddo

           do j=1,12
              write(12,rec=irec) sngl(sivolm(4,j)-1.d-20)
              irec=irec+1
           enddo

           do j=1,12
              write(12,rec=irec) sngl(LIdiffm(1,j)-1.d-20)
              irec=irec+1
           enddo

           do j=1,12
              write(12,rec=irec) sngl(LIdiffm(4,j)-1.d-20)
              irec=irec+1
           enddo

           do j=1,12
              write(12,rec=irec) sngl(LSIdiffm(1,j)-1.d-20)
              irec=irec+1
           enddo

           do j=1,12
              write(12,rec=irec) sngl(LSIdiffm(4,j)-1.d-20)
              irec=irec+1
           enddo

           do j=1,12
              write(12,rec=irec) sngl(inewn2*(rhoLice/rhoW)/1.d+6)
              irec=irec+1
           enddo

        do i=1,4
           do j=1,12
              write(12,rec=irec) sngl(hlhm(i,j)*Aarea(i)/ac1/1.d+15)
              irec=irec+1
           enddo
        enddo
c ***
 
           do j=1,12
              write(12,rec=irec) sngl(co2am(j)-1.d-20)
              irec=irec+1
           enddo

           do j=1,12
              write(12,rec=irec) sngl(ctavm(j)-1.d-20)
              irec=irec+1
           enddo

           do j=1,12
              write(12,rec=irec) sngl(PO4avm(j)-1.d-20)
              irec=irec+1
           enddo

           do j=1,12
              write(12,rec=irec) sngl(ATavm(j)-1.d-20)
              irec=irec+1
           enddo

           do j=1,12
              write(12,rec=irec) sngl(sprodm(j)*389028/1.d+15)
              irec=irec+1
           enddo

           do j=1,12
              write(12,rec=irec) sngl(hprodm(j)*389028/1.d+15)
              irec=irec+1
           enddo

           do j=1,12
              write(12,rec=irec) sngl(totprodm(j)*389028/1.d+15)
              irec=irec+1
           enddo

        do i=1,8
           do j=1,12
              write(12,rec=irec) sngl(ctm(i,j))
              irec=irec+1
           enddo
        enddo

        do i=1,8
           do j=1,12
              write(12,rec=irec) sngl(atm(i,j))
              irec=irec+1
           enddo
        enddo

        do i=1,8
           do j=1,12
              write(12,rec=irec) sngl(hco3m(i,j))
              irec=irec+1
           enddo
        enddo


        do i=1,8
           do j=1,12
              write(12,rec=irec) sngl(co3m(i,j))
              irec=irec+1
           enddo
        enddo

        do i=1,8
           do j=1,12
              write(12,rec=irec) sngl(co2m(i,j))
              irec=irec+1
           enddo
        enddo

        do i=1,8
           do j=1,12
              write(12,rec=irec) sngl(co2aam(i,j))
              irec=irec+1
           enddo
        enddo

        do i=1,4
           do j=1,12
              write(12,rec=irec) sngl(ctfluxm(i,j)*389028/1.d+15)
              irec=irec+1
           enddo
        enddo

        do i=1,8
           do j=1,12
              write(12,rec=irec) sngl(am(i,j))
              irec=irec+1
           enddo
        enddo

        do i=1,8
           do j=1,12
              write(12,rec=irec) sngl(po4m(i,j))
              irec=irec+1
           enddo
        enddo

        do i=1,8
           do j=1,12
              write(12,rec=irec) sngl(rainm(i,j))
              irec=irec+1
           enddo
        enddo

        do i=1,8
           do j=1,12
         write(12,rec=irec) sngl(prodm(i,j)*389028*RedC(i)/1.d+15)
         irec=irec+1
           enddo
        enddo

ccc hezi 16.9.1999 add the following
 
         do j=1,12
           write(12,rec=irec) sngl(tf21m(j)*CwrhoW/1.d+15)
           irec=irec+1
         enddo

         do j=1,12
           write(12,rec=irec) sngl(tf15m(j)*CwrhoW/1.d+15)
           irec=irec+1
         enddo

         do j=1,12
         write(12,rec=irec) sngl(-1.*tdf12m(j)*CwrhoW/1.d+15)
         irec=irec+1
         enddo

         do j=1,12
           write(12,rec=irec) sngl(tdf15m(j)*CwrhoW/1.d+15)
           irec=irec+1
         enddo

         do i=1,4
           do j=1,12
           write(12,rec=irec) sngl(vmixm(i,j)/1.d+6)
           irec=irec+1
           enddo
         enddo

         do i=1,4
           do j=1,12
           write(12,rec=irec) sngl(RedCm(i,j))
           irec=irec+1
           enddo
         enddo

           do j=1,12
              write(12,rec=irec) sngl(ctf21m(j)*389028/1.d+15)
              irec=irec+1
           enddo

           do j=1,12
              write(12,rec=irec) sngl(ctf15m(j)*389028/1.d+15)
              irec=irec+1
           enddo
           
           do j=1,12
              write(12,rec=irec) sngl(ctdf56m(j)*389028/1.d+15)
              irec=irec+1
           enddo
           do j=1,12
              write(12,rec=irec) sngl(ctf34m(j)*389028/1.d+15)
              irec=irec+1
           enddo
           do j=1,12
              write(12,rec=irec) sngl(ctf48m(j)*389028/1.d+15)
              irec=irec+1
           enddo
           do j=1,12
              write(12,rec=irec) sngl(ctf87m(j)*389028/1.d+15)
              irec=irec+1
           enddo
           do j=1,12
              write(12,rec=irec) sngl(ctf62m(j)*389028/1.d+15)
              irec=irec+1
           enddo
           do j=1,12
              write(12,rec=irec) sngl(ctf73m(j)*389028/1.d+15)
              irec=irec+1
           enddo
           do j=1,12
              write(12,rec=irec) sngl(ctdf15m(j)*389028/1.d+15)
              irec=irec+1
           enddo
           do j=1,12
              write(12,rec=irec) sngl(ctdf26m(j)*389028/1.d+15)
              irec=irec+1
           enddo
           do j=1,12
              write(12,rec=irec) sngl(ctdf48m(j)*389028/1.d+15)
              irec=irec+1
           enddo
           do j=1,12
              write(12,rec=irec) sngl(ctdf12m(j)*389028/1.d+15)
              irec=irec+1
           enddo
           do j=1,12
              write(12,rec=irec) sngl(ctdf34m(j)*389028/1.d+15)
              irec=irec+1
           enddo
           do j=1,12
              write(12,rec=irec) sngl(ctdf56m(j)*389028/1.d+15)
              irec=irec+1
           enddo
           do j=1,12
              write(12,rec=irec) sngl(ctdf78m(j)*389028/1.d+15)
              irec=irec+1
           enddo


           do j=1,12
              write(12,rec=irec) sngl(po4f21m(j)*389028/1.d+15)
              irec=irec+1
           enddo

           do j=1,12
              write(12,rec=irec) sngl(po4f15m(j)*389028/1.d+15)
              irec=irec+1
           enddo
           
           do j=1,12
              write(12,rec=irec) sngl(po4df56m(j)*389028/1.d+15)
              irec=irec+1
           enddo
           do j=1,12
              write(12,rec=irec) sngl(po4f34m(j)*389028/1.d+15)
              irec=irec+1
           enddo
           do j=1,12
              write(12,rec=irec) sngl(po4f48m(j)*389028/1.d+15)
              irec=irec+1
           enddo
           do j=1,12
              write(12,rec=irec) sngl(po4f87m(j)*389028/1.d+15)
              irec=irec+1
           enddo
           do j=1,12
              write(12,rec=irec) sngl(po4f62m(j)*389028/1.d+15)
              irec=irec+1
           enddo
           do j=1,12
              write(12,rec=irec) sngl(po4f73m(j)*389028/1.d+15)
              irec=irec+1
           enddo
           do j=1,12
              write(12,rec=irec) sngl(po4df15m(j)*389028/1.d+15)
              irec=irec+1
           enddo
           do j=1,12
              write(12,rec=irec) sngl(po4df26m(j)*389028/1.d+15)
              irec=irec+1
           enddo
           do j=1,12
              write(12,rec=irec) sngl(po4df48m(j)*389028/1.d+15)
              irec=irec+1
           enddo
           do j=1,12
              write(12,rec=irec) sngl(po4df12m(j)*389028/1.d+15)
              irec=irec+1
           enddo
           do j=1,12
              write(12,rec=irec) sngl(po4df34m(j)*389028/1.d+15)
              irec=irec+1
           enddo
           do j=1,12
              write(12,rec=irec) sngl(po4df56m(j)*389028/1.d+15)
              irec=irec+1
           enddo
           do j=1,12
              write(12,rec=irec) sngl(po4df78m(j)*389028/1.d+15)
              irec=irec+1
           enddo           
        return
        end
