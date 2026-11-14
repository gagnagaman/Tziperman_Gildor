C  RESTART IS TRUE IF RESTARTING FROM PREVIOUS RUN, OTHERWISE FALSE
C  NSTART IS THE DAY NUMBER OF RESTART (IT IS ONLY MEANINGFUL WITH START=TRUE)
C  NDBEG AND NDEND ARE START AND END TIMES IN DAYS.
C  NBEG AND NEND ARE THE CORRESPONDING TIME STEP NUMBERS.
C  NDINC IS THE NUMBER OF DAYS TO BE RUN.
C  NDPRIN IS THE PRINT INTERVAL IN DAYS.
C  NDSTOR IS THE STORING INTERVAL IN DAYS.
C  ISMTH IS THE NUMBER OF TIME STEPS BETWEEN TIME SMOOTHING 
c  dtcouple is the time step 

      LOGICAL RESTART,STORIT,PRNTIT,coupleflux,seaicectr,
     &        landicectr,geochemctr,bioctr,MFLAG

      real*8 L(4),width(8),SSTdata(4),Sdata(4),Sref,Tref,fl(4),
     &       fi0(4),fsi0(4),CO2ref,CO2in

      real*8 rhoW,alpha,beta,G,sigma,RR,CpA,CpW,Ps,epsilon,CwRhoW,
     &       RhoIce, Lice,LiSp,LiNp,Tper(4),CTp(4),LH

      real*8 rtime,tdt,dt,dtcouple,fo(4),Tav,CTav,totalArea,
     &       tetaav,Sav

      character*32 FILRESI, FILRESO, FILPLT, FILCNT

      integer  NDBEG, NDINC, NSTART, NDPRIN, NDSTOR,ndend,kday,
     1         ISMTH, k,km,kp,ks,NSMTH1,NDSAVE,i,nstep,nmstep,
     1         ndstep,nbeg,nend,nprint,nstore,nsstp,nsmth,npstp,
     1         irec,ndat,j,m,NSMILAN,NSMILAN1,NSMILAN2,NMSTEP2,
     1         MID_JUNE

       common /cntrl1/RESTART,STORIT,PRNTIT,coupleflux,NSMILAN,NSMILAN2,
     &        NDBEG, NDINC, NSTART, NDPRIN, NDSTOR,ndend,kday,
     &        ISMTH, k,km,kp,ks,NSMTH1,NDSAVE,nstep,MFLAG,nmstep,
     &        FILPLT, FILCNT,seaicectr,landicectr,geochemctr,bioctr



       common /cntrl2/ndstep,nbeg,nend,nprint,nstore,nsstp,nsmth,npstp,
     &        irec,ndat,tdt,dt,dtcouple,NMSTEP2

       common /const/rhoW,alpha,beta,G,sigma,RR,CpA,CpW,Ps,epsilon,
     &               CwRhoW,RhoIce, Lice,LiSp,LiNp,Tper,CTp,LH
 
       common /parameters1/L,width,SSTdata,Sdata,Sref,Tref,fl,
     &       fi0,fsi0,CO2ref,fo,CO2in,m,MID_JUNE


       common /globalvar/Tav,CTav,totalArea,tetaav,Sav
