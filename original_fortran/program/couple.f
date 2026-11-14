c     =====================================
c     driver program for the coupled box model
c     =====================================
c
c     MKS units

      implicit none
      real*8 state,unesco
      include 'insolation1.h'
      include 'atmosphere.h'
      include 'ocean.h'
      include 'couple.h'
      include 'seaice.h'
      include 'landice.h'
      include 'biochemical.h'
      include 'averaging.h'
      include 'biochem_coef.h'
      real*8 t_state

C  NSTART IS THE DAY NUMBER OF RESTART (IT IS ONLY MEANINGFUL WITH START=TRUE)
C  NDBEG AND NDEND ARE START AND END TIMES IN DAYS.
C  NDINC IS THE NUMBER OF DAYS TO BE RUN.
C  NDPRIN IS THE PRINT INTERVAL IN DAYS.
C  NDSTOR IS THE STORING INTERVAL IN DAYS.
C  ISMTH IS THE NUMBER OF TIME STEPS BETWEEN TIME SMOOTHING  

       NAMELIST /CONTRL/ RESTART,NDBEG,NDINC,NSTART,NDPRIN,NDSTOR,
     &               ISMTH,dtcouple

       NAMELIST /CONFIGURE/ coupleflux,seaicectr,landicectr,
     &                      geochemctr,bioctr 

C  INPUT FILE IF RESTARTING,OUTPUT RESTART FILE, OUTPUT PLOTTING FILE
       NAMELIST /FILES/  FILRESI, FILRESO, FILPLT, FILCNT

c SOME PARAMETERS:
c L is the length of atmospheric box
c width is the width of atmospheric box
c SSTdata - sst input for restoring b.c
c Sdata   - salinity input for restoring b.c
c Sref & Tref are for the linear equation of state
c fl-width fraction which is land
c fi0-meridional fraction which is land ice at time zero
c fsi0-meridional fraction which is sea ice at time zero
c CO2ref - used in calculating radiative forcing of CO2

       namelist /parameters/L,width,SSTdata,Sdata,Sref,Tref,fl,
     &                      fi0,fsi0,CO2ref,CO2in

c PHYSICAL CONSTANTS:
c rhoW - density of water
c alpha & beta are used in the linear equation of state. i took them from Eli.
c G is gravity and sigma Stefan-Boltzman constants (W/m^2*K^4)
c RR - Gas constant for dry air (J/Kg*K)
c CpA - atmospheric specific heat at constant pressure (J/Kg*K). From Ilya.
c Ps- surface pressure taken as 1000mb
c epsilon - RR/CpA
c Lice fusion heat of ice. 
c Rhoice- density of ice
c LIf
c rhoLice

       namelist /constants/rhoW,alpha,beta,G,sigma,RR,CpA,CpW,Ps,
     &             epsilon,CwRhoW,Lice,Rhoice,LH,rhoLice

C OCEAN PARAMETERS
c D are the depths of the ocean boxes
c r is friction coefficient
c gammaS & gammaT are restoring coefficient for salinity and temperature
c **section are coefficient to convert circulation into Sv
c Kv - vertical diffusion coefficient, Kh for horizontal diffusion

       namelist /Oparameters/D,r, gammaS, gammaT,uppersection, 
     &           lowersection,horisection,Kv,Kh,dsigma,dsigma2

C ATMOSPHERIC PARAMETERS
c Gama0 - coefficient for calculating heat flux between ocean and atmosphere (Hany type) 
c albedoL,albedoI,albedoSI- albedos of Land, land-ice and sea-ice
c Solar - area weighted average solar radiation
c KT & KW are for calculating meridional transport of heat and humidity
c Kapa - for calculating radiative forcing of CO2
c Plw0 - background value of emisivity
c RH - relative humidity

       namelist /Aparameters/Gama0,albedoL,albedoI,albedoSI,albedoS,
     &           albedoC,Solar,KT,KWa,Kwa2,Kapa,Plw0,RH

c ICE PARAMETERS
c Di - thickness of sea-ice created 
c Tice - temperature in which the model create ice
c Igama - gama in the land-ice model. e.g Ghil. include the yield stress
c TLice1 & TLice2 are for use in calculating precipitation distribution function
c Icoef1 - restoring coefficient for calculating heat flux between atmosphere 
c          and land-ice
c IceHc - critical height for calculating precipitation distribution function

       namelist /SeaIcepar/Di,Tice

       namelist /LandIcepar/Igama,TLice1,TLice2,Icoef1,IceHc,T12preci,
     &           T32preci,T34preci,alphaP,T1melt,T4melt


c CT0 & AT0 - total C02 and total alkalinity at time zero
c pv - coefficient for gas exchange between atmosphere and ocean
c Rc 

       namelist /BioChem/CT0,AT0,pv0,pc10,PO40,O20,RedC0,RedN,RedO2,h

       namelist /Perturbation/LiSp,LiNp,Tper,CTp

C  READ NAME LIST PARAMETER
      OPEN(UNIT=8,FILE='coupledata.list',status='old')
      READ (8,CONTRL)
      READ (8,CONFIGURE)
      READ (8,FILES)
      READ (8,parameters)
      READ (8,constants)
      READ (8,Oparameters)
      READ (8,Aparameters)
      READ (8,SeaIcepar)
      READ (8,LandIcepar)
      READ (8,BioChem)
      READ (8,Perturbation)
      CLOSE(UNIT=8)

      print *,'dsigma=',dsigma
c OPEN OUTPUT FILE FOR GRAD
             OPEN(UNIT=12,FILE=FILPLT
     1,            FORM='UNFORMATTED'
     1,            ACCESS='DIRECT'
     1,            RECL=1*1)
c     1,            RECL=1*4)

             irec=1
             ndat=0

c =============================================
C  RESTART OR INITIALIZE ARRAYS
c==============================================
        IF(RESTART) THEN
           TDT=2.*dtcouple
           
c          IceVolS & IceVolN are land-ice volume in the south and north boxes
c          ILenS & ILenN are meridional extension of land-ice sheet
c          fi - sea-ice fraction

           OPEN(10,FILE=FILRESI,FORM='FORMATTED')
           READ(10,*,END=9876) NDSAVE,NSMTH1,KP,K,KM,Teta,t,s,fsi,
     &      rho,fi,ILenS,ILenN,IceVolS,IceVolN,CT,AT,a,Catmos,rhos,rhot,
     &      CO2a,PO4,Idepth,Tav,Sav,tetaav,q21,q34,q23,fw21,fw34,
     &      ft21,ft34,hlh(1),hlh(4),hao(1),hao(4),q21s,q21t,q34s,q34t,
     &      IceSflux,Ipnew,SIvol,Inewf,Imelt,Im,Iflow,InewTemp,Havail,
     &      FW1,FW4



c          hezi 23.12
c          CO2a=CO2in
           
           IF(NDSAVE.EQ.NSTART) GOTO 7000
 9876      continue
           PRINT *,' ** ERROR ** COULD NOT FIND CORRECT RESTART DAY'
           CLOSE(10)
           STOP

 7000      continue
           print *,' '
           PRINT *,' ===  SUCCESS TO READ IN RESTART DATA'
           print *,' '
           print *, 'FILRESI=',FILRESI
           print *,' '
           print *,' NDSAVE ,      NSMTH1 ,      KP ,       K ,      KM'
           print *,NDSAVE,NSMTH1,KP,K,KM
           CLOSE(10)
           NDBEG=NDSAVE

c          swap k values
           ks=km
           km=k
           k=kp
           kp=ks



        ENDIF

c hezi: tmporary
c        do i=1,8
c           do j=1,3
c              CT(i,j)=CT0
c              AT(i,j)=AT0
c           end do
c        end do

c hezi: tmporary
        do i=1,4
              RedC(i)=RedC0
              pc1(i)=pc10(i)
              pv(i)=pv0
        end do
        do i=1,8
          k0(i)=1.
        enddo


        
        IF(.NOT.RESTART) THEN 
           TDT=dtcouple

c          initialized arrays
c          atmosphere and ice variable
           do i=1,4
              IceSflux(i)=0.0
              Ipnew(i)=0.0
              do k=1,3
                 Teta(1,k)=273.
                 Teta(2,k)=283.
                 Teta(3,k)=283.
                 Teta(4,k)=273.
                 Catmos(k)=5.7d+16
                 SIvol(i,k)=fsi0(i)*L(i)*Owidth(i)*Di(i)
                 fsi(i,k)=fsi0(i)
                 fi(i,k)=fi0(i)
                 ILenS(k)=fi(1,k)*L(1)
                 ILenN(k)=0.5*fi(4,k)*L(4)
         IceVolS(k)=fl(1)*width(1)*sqrt(Igama)*(ILenS(k)**(3./2.))
         IceVolN(k)=2*fl(4)*width(4)*sqrt(Igama)*(ILenN(k)**(3./2.))
              end do
              Idepth(i)=Di(i)
              Inewf(i)=0.0
              Imelt(i)=0.0
              Im(i)=0.0
              Iflow(i)=0.0
              InewTemp(i)=0.0
              Havail(i)=0.0
              Inew(i)=0.0
           enddo

           CO2a=CO2in

c          ocean variables
           do i=1,8
              do k=1,3
                 s(i,k)=35.
                 t(i,k)=278.
                 a(i)=1.d+8
                 CT(i,k)=CT0
                 AT(i,k)=AT0
                 PO4(i,k)=PO40
                 O2(i,k)=O20
              end do
              t_state=t(i,3)-273.0
c              rho(i)=state(t_state,s(i,3))
              rho(i)=unesco(s(i,3),t_state)

           end do

           KP=3
           K= 2
           KM=1
           NSMTH1=0
        ENDIF


C  DEFINE VARIABLES DEPENDING ON RUN PARAMETERS
      NDEND=NDBEG+NDINC
      NDSTEP=INT(86400./DTcouple +.5) ! number of time steps each day
      NBEG=NDBEG
      NEND=NDEND*NDSTEP
c hezi start to calculate monthely mean fron ndprint -365
      NPRINT=NDPRIN*NDSTEP
      NSTORE=NDSTOR*NDSTEP

c calculate some constant for use in atmosphere.f & ocean.f. 

      ac1=(2.**(RR/CpA))*G/(Ps*CpA)     ! 2^(R/Cp)*g/(Ps*Cp)

      totvolume=0.0
      totalArea=0.0
      do i=1,4
         Aarea(i)=L(i)*width(i)
         totalArea=totalArea+Aarea(i)
         fo(i)=1.-fl(i)
         Owidth(i)=fo(i)*width(i)
         recL(i)=1./L(i)
         recD(i)=1./D(i)
         vol(i)=L(i)*D(i)*Owidth(i)
         recvol(i)=1./vol(i)
         totvolume=totvolume+vol(i)
         Oarea(i)=fo(i)*Aarea(i)
         OWarea(i)=(1.-fsi(i,k))*Oarea(i)
      enddo

      do i=5,8
         Owidth(i)=Owidth(i-4)
         Oarea(i)=Oarea(i-4)
         recL(i)=1./L(i-4)
         recD(i)=1./D(i)
         vol(i)=L(i-4)*D(i)*Owidth(i)
         recvol(i)=1./vol(i)
         totvolume=totvolume+vol(i)
      enddo

c calculate coefficients needed when calculating the meridional transport 
c between ocean boxes

      c21=D(5)*G/(r*rhoW*(D(1)+D(5))*(L(1)+L(2)))
      c23=D(6)*G/(r*rhoW*(D(2)+D(6))*(L(2)+L(3)))
      c34=D(8)*G/(r*rhoW*(D(4)+D(8))*(L(3)+L(4)))

      c56=c21*D(1)/D(5)
      c76=c23*D(2)/D(6)
      c87=c34*D(4)/D(8)

      c21=c21*uppersection
      c23=c23*uppersection
      c34=c34*uppersection

      c56=c56*lowersection
      c76=c76*lowersection
      c87=c87*lowersection

      rhoIw=RhoIce/rhoW
      rhoLIw=RhoLice/rhoW
 
c for seaice

      rhoILice=rhoice*Lice
      do i=1,4
         DiOarea(i)=Di(i)*Oarea(i)
      enddo


c insert perturbation

c          fi(1)=fi(1)+LiSp
c          fi(4)=fi(4)+LiNp
c      do i=1,4
c         T(i,k)=T(i,k)+Tper(i)
c         CT(i,k)=CT(i,k)+CTp(i)
c      enddo

c some variables relating to Land ice

c      ISmax=3.04d+16*sqrt(Igama)   ! maximux land-ice volume
c      INmax=5.38d+16*sqrt(Igama)

      rhowLice=rhoW/rhoLice

      LIc1=fl(1)*width(1)*sqrt(Igama)
      LIc4=fl(4)*width(4)*sqrt(Igama)

      ISmax=(L(1)**(3./2.))*LIc1
      INmax=((0.5*L(4))**(3./2.))*2*LIc4

c variables and calculation related to insolation.f
      pi=3.14159265358979D0
      radius=6.368d+6
      pi2rr=2*pi*radius*radius
      do i=1,36
         phi=-87.5+(i-1)*5
         belt(i)=pi2rr*
     &   abs((sin((phi-2.5)*pi/180.)-sin((phi+2.5)*pi/180.)))
      enddo
      do i=1,4
         boxarea(i)=0.0
         do j=1,9
           boxarea(i)=boxarea(i)+belt((i-1)*9+j)

         enddo
           print *,'boxarea no ',i,' = ',boxarea(i)
      enddo

      DATA NJM/31,28,31,30,31,30,31,31,30,31,30,31/
      DATA NJM1/31,59,90,120,151,181,212,243,273,304,334,365/

      open(9,file="insolation.input",form='formatted')

      PIR=PI/180.0D0                                                    
      PIRR=PIR/3600.0D0                                                 
      STEP=360.0D0/365.25D0                                             
      TEST=0.0001D0


C                                                                       
C   1.EARTH ORBITAL ELEMENTS : ECCENTRICITY           ECC   TABLE 1     
C***************************   PRECESSIONAL PARAMETER PRE               
C                              OBLIQUITY              XOB   TABLE 2     
C                              GENERAL PRECESSION     PRG               
C                              LONGITUDE PERIHELION   PERH  TABLE 3     
C
C         NEF  NOB  NOP  MAY BE REDUCED TO  19  18  9                   
C         BUT THE INPUT DATA MUST BE CHANGED ACCORDINGLY 
C                                                                       
      READ(9,5001)                                                      
5001  FORMAT(/////)                                                     
C         READ AMPLITUDE A  MEAN RATE B  PHASE C                        
C         THEY ARE IMMEDIATELY CONVERTED IN RADIANS                
C                                                                       
C   ECCENTRICITY table 4 from JAS 1978
C                                                                       
      NEF=19                                                            
      DO I=1,NEF                                                      
      READ(9,5000) AE(I),Y,Z                                            
 5000 FORMAT (13X,F11.8,F20.7,F20.6)                                    
      BE(I)=Y*PIRR                     ! convert into radians                                  
      CE(I)=Z*PIR                                                       
      enddo
C                                                                       
C   OBLIQUITY      table 1 in JAS                                                     
C                                                                       
      XOD=23.320556D0                                                   
      NOB=47                                                            
      DO I=1,NOB                                                      
      READ(9,5002) AOB(I),Y,Z                                           
 5002 FORMAT(7X,F13.7,2X,F10.6,2X,F10.4)                                
      BOB(I)=Y*PIRR                                                     
      COB(I)=Z*PIR                                                      
      enddo
C                                                                       
C   GENERAL PRECESSION IN LONGITUDE table 5 in JAS                                 
C                                                                       
      XOP=3.392506D0                                                    
      PRM=50.439273D0                                                   
      NOP=78                                                            
      DO I=1,NOP                                                     
      READ(9,5002) AOP(I),Y,Z                                           
      BOP(I)=Y*PIRR                                                     
      COP(I)=Z*PIR                                                      
      enddo
C                                                                       
C    2.INPUT PARAMETERS :  TIME T ; when investigation seasonal cycle alone,  T=0

      T1=0.0
C                                                                        
C      NEFF  NOBB  NOPP  NUMBER OF TERMS KEPT FOR COMPUTATION OF EARTH ORBITAL ELEMENTS 
                       
C      THEY CAN BE REDUCED TO 19,18,9 RESPECTIVELY   <- this is what we do                     
C          FOR A MINIMUM ACCURACY                                       
C
      NEFF=19
      NOBB=18
      NOPP=9

C                                                                       
C   3.NUMERICAL VALUE FOR ECC PRE XOB                                   
C     T IS NEGATIVE FOR THE PAST ;  T IS IN 1000 YEARS ; T=0 in seasonal cycles mode                
C                                                                       
      T1=T1*1000.0D0                                                      
      XES=0.0D0                                                         
      XEC=0.0D0                                                         
      DO I=1,NEFF                                                     
      ARG=BE(I)*T1+CE(I)                                                 
      XES=XES+AE(I)*DSIN(ARG)                                           
      XEC=XEC+AE(I)*DCOS(ARG)                                           
      enddo
      ECC=DSQRT(XES*XES+XEC*XEC)                                        
      TRA=DABS(XEC)                                                     
      IF(TRA.LE.1.0D-08) GO TO 10                                       
      RP=DATAN(XES/XEC)                                                 
      IF(XEC) 11,10,12                                                  
   11 RP=RP+PI                                                          
      GO TO 13                                                          
   12 IF(XES) 14,13,13                                                  
   14 RP=RP+2.0D0*PI                                                    
      GO TO 13                                                          
   10 IF(XES) 15,16,17                                                  
   15 RP=1.5D0*PI                                                       
      GO TO 13                                                          
   16 RP=0.0D0                                                          
      GO TO 13                                                          
   17 RP=PI/2.0D0                                                       
   13 PERH=RP/PIR                                                       
C                                                                       
      PRG=PRM*T1                                                         
      DO I=1,NOP                                                      
      ARG=BOP(I)*T1+COP(I)                                               
      PRG=PRG+AOP(I)*DSIN(ARG)                                          
      enddo
      PRG=PRG/3600.0D0+XOP                                              
      PERH=PERH+PRG                                                     
   54 IF(PERH) 51,55,53                                                 
   51 PERH=PERH+360.0D0                                                 
      GO TO 54                                                          
   53 IF(PERH.LT.360.0D0) GO TO 55                                      
      PERH=PERH-360.0D0                                                 
      GO TO 53                                                          
   55 CONTINUE                                                          
C                                                                       
      PRE=ECC*DSIN(PERH*PIR)                                            
C                                                                       
      XOB=XOD                                                           
      DO I=1,NOBB                                                     
      ARG=BOB(I)*T1+COB(I)                                               
      XOB=XOB+AOB(I)/3600.0D0*DCOS(ARG)                                 
      enddo

C                                                                       
C   4.DAILY INSOLATION                                                  
C*********************                                                  
C                                                                       
C HEZI: use CALENDAR DATE                                   
C       DAILY INSOLATION IN LY DAY(-1)    OR    KJ M(-2) DAY(-1)        
C                  IF S0 IN LY MIN(-1)    OR    W M(-2)                 
C                     TAU = 24*60 MIN     OR    24*60*60 SEC / 1000     
C                                                                       
C                        IN W M(-2)                                     
C                  IF S0 IN W M(-2)      AND    TAU=1.0                 
C                                                                       
      SS=1353.0D0      ! solar constant                                                 
      TAU=86.4D0                                                        
      SF=TAU*SS/PI                                                      
      SO=DSIN(XOB*PIR)                                                  
      XL=PERH+180.0D0                                                   

c end of insolation initialization


c  DEFINE VARIABLES DEPENDING ON RUN PARAMETERS NEEDED FOR 
c  CALCULATING MONTHLY MEANS
      do i=1,12
         NJMS(i)=NJM(i)*NDSTEP
         NJM1S(i)=NJM1(i)*NDSTEP
         NJM2S(i)=NJM1(i)-INT(0.5*NJM(i))*NDSTEP
         write(6,*) NJM(i),' ',njms(i),'  ',NJM1(i),' ',
     &              njm1s(i),' ',njm2s(i)
      enddo

      MFLAG=.false.

C  PRINT OUT CONFIGURATION PARAMETERS
      write(*,'(/,1x,"CONFIGURATION PARAMETERS:"/,
     1         1X,"coupleflux=",L8/,
     1         1X,"seaicectr= ",L8/,
     1         1X,"landicectr=",L8/,
     1         1X,"geochemctr=",L8/,
     1         1X,"bioctr=    ",L8/)')
     &      coupleflux,seaicectr,landicectr,
     1      geochemctr,bioctr 



C  BEGIN TIME-STEPPING 
      NPSTP=0
      NSSTP=0
      NSMTH=NSMTH1
      NSTEP=0
      KDAY=NDBEG     
      NDSAVE=KDAY
      RTIME=KDAY
      NMSTEP=0  

      call insolation1


c     ------------
c     step in time:
c     -------------

 3000 continue

      call atmosphere

      call ocean

      if (geochemctr) then
         call biochem_coef
         call biochemical
      endif

C  INCREMENT TIME PARAMETERS
      NSTEP=NSTEP+1
      NPSTP=NPSTP+1
      NSSTP=NSSTP+1
      NSMTH=NSMTH+1
c hezi, add for monthly means
      NMSTEP=NMSTEP+1

      RTIME=RTIME+DT/86400.
      IF(MOD(NSTEP,NDSTEP).EQ.0) THEN ! if another day elapsed
        NSTEP=0
        KDAY=KDAY+1
        NDSAVE=KDAY
        RTIME=KDAY
      ENDIF

c hezi 22.11.98: test every 4 days
c calculate insolation every two days
      IF(MOD(NSTEP,NDSTEP).EQ.0) THEN ! if another day elapsed
         call insolation1
      ENDIF

C  SMOOTH FIELDS, SM2DT FOR 2DT NOISE
c      IF(MOD(NSMTH,ISMTH).EQ.0) THEN
c        CALL ASM2DT
c        CALL OSM2DT
c        NSMTH=0
c      ENDIF
c      NSMTH1=NSMTH 
c hezi 9 september
        CALL ASM2DT1
        CALL OSM2DT1

        

C  CHECK FOR PRINTING FIELDS
      PRNTIT=.FALSE.
      IF((MOD(KDAY,NDPRIN).EQ.0).and.(nstep.eq.1)) THEN  ! ****  PRINT FIELDS
c        NPSTP=0        ! hezi: zero it in averaging.f as well
        PRNTIT=.TRUE.
        MFLAG=.TRUE.
c        print *,'save at ',kday
        NMSTEP=1

        call zeromean


      ENDIF  

       IF(MFLAG) THEN  ! STORE FOR PRINTING
c hezi 25/9 : calculate some global mean values
c 1) average total co2 in the ocean
         ATav=0.0
         PO4av=0.0
         CTav=0.0

         do i=1,8
           CTav=CTav+vol(i)*CT(i,kp)
           PO4av=PO4av+vol(i)*PO4(i,kp)
           ATav=ATav+vol(i)*AT(i,kp)
         enddo
         CTav=CTav/totvolume
         PO4av=PO4av/totvolume
         ATav=ATav/totvolume
c 2) average atmospheric temperature
         Tetaav=0.0
         do i=1,4
            Tetaav=Tetaav+teta(i,kp)*Aarea(i)
         enddo
         Tetaav=(Tetaav/totalArea)-273.

c 3) average ocean temperature and salinity
         Tav=0.0
         Sav=0.0
         do i=1,8
            Tav=Tav+T(i,kp)*vol(i)
            Sav=Sav+S(i,kp)*vol(i)
         enddo
         Tav=(Tav/totvolume)-273.
         Sav=(Sav/totvolume)

c        ndat=ndat+1
        call averaging

      ENDIF

c hezi 22.6
C  CHECK FOR   STORING FIELDS
      STORIT=.FALSE.
      IF(NSTORE.NE.0.AND.MOD(NSSTP,NSTORE).EQ.0) THEN      
C ! ****  STORE FIELDS
        NSSTP=0
        STORIT=.TRUE.
        OPEN(UNIT=11,FILE=FILRESO,FORM='FORMATTED')
        WRITE(11,*) NDSAVE,NSMTH1,KP,K,KM,Teta,t,s,fsi,
     &      rho,fi,ILenS,ILenN,IceVolS,IceVolN,CT,AT,a,Catmos,rhos,rhot,
     &      CO2a,PO4,Idepth,Tav,Sav,tetaav,q21,q34,q23,fw21,fw34,
     &      ft21,ft34,hlh(1),hlh(4),hao(1),hao(4),q21s,q21t,q34s,q34t,
     &      IceSflux,Ipnew,SIvol,Inewf,Imelt,Im,Iflow,InewTemp,Havail,
     &      FW1,FW4

        print *,'store at ',NDSAVE,NSMTH1,KP,K,KM

        close (11)

      ENDIF

      IF(.NOT.(KDAY.EQ.NDEND.OR.NSTEP.EQ.NEND))THEN
        TDT=2.*DTcouple
        KS=KM
        KM=K
        K=KP
        KP=KS

        GOTO 3000  ! <======================== DO ANOTHER TIME STEP
      ELSE
        GOTO 9000  ! <======================== DONE
      ENDIF

9000  CONTINUE   ! ****   PROGRAM COMPLETED 

c prepare control file for GRADS

       call p4grads
C

c store restart data


        OPEN(UNIT=11,FILE=FILRESO,FORM='FORMATTED')
        WRITE(11,*) NDSAVE,NSMTH1,KP,K,KM,Teta,t,s,fsi,
     &      rho,fi,ILenS,ILenN,IceVolS,IceVolN,CT,AT,a,Catmos,rhos,rhot,
     &      CO2a,PO4,Idepth,Tav,Sav,tetaav,q21,q34,q23,fw21,fw34,
     &      ft21,ft34,hlh(1),hlh(4),hao(1),hao(4),q21s,q21t,q34s,q34t,
     &      IceSflux,Ipnew,SIvol,Inewf,Imelt,Im,Iflow,InewTemp,Havail,
     &      FW1,FW4

        print *,NDSAVE,NSMTH1,KP,K,KM

        stop
        end

