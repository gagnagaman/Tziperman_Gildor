      subroutine insolation1

      implicit none
      include 'insolation1.h'
      include 'atmosphere.h'
      include 'couple.h'
      include 'ocean.h'

C*******************************************                            
C   DAILY INSOLATION - LONG TERM VARIATION *                            
C*******************************************                            
C                                                                       

c subroutine for calculation of daily insolation. based on Berger 1978 
c but modified by me.
c we calculate the insolation in latitudes belt of 5 degree and averaged
c for each box.


                                             

C   4.2 CALENDAR DATE  MA-JA                                            
C--------------------                                                   
C      ND  NUMBER OF THIS DAY IN A YEAR OF 365 DAYS                     
C      XLAM = MEAN LONG. SUN FOR TRUE LONG. = 0                         
C      DLAMM = MEAN LONG. SUN FOR MA-JA                                 
C                                                                       
c determaine month and day
c      WRITE(18,6100)
c      WRITE(18,6003) 
c 6003 FORMAT(1X,'   Day   ',' Latitude   ',' Berger     ') 

c      WRITE(18,6100)                                                     
c 6100 FORMAT(/)

      ND=mod(kday,365)
                                             
      XLLP=XL*PIR                                                       
      XEE=ECC*ECC                                                       
      XSE=DSQRT(1.0D0-XEE)                                              
      XLAM=(ECC/2.0D0+ECC*XEE/8.0D0)*(1.0D0+XSE)*DSIN(XLLP)-XEE/4.0D0*(0
     & .5D0+XSE)*DSIN(2.0D0*XLLP)+ECC*XEE/8.0D0*(1.0D0/3.0D0+XSE)*DSIN(3.
     &  0D0*XLLP)                                                         
      XLAM=2.0D0*XLAM/PIR                                               
      DLAMM=XLAM+(ND-80)*STEP                                           
      ANM=DLAMM-XL                                                      
      RANM=ANM*PIR                                                      
      XEC=XEE*ECC                                                       
      RANV=RANM+(2.0D0*ECC-XEC/4.0D0)*DSIN(RANM)+5.0D0/4.0D0*ECC*ECC*   
     &     DSIN(2.0D0*RANM)+13.0D0/12.0D0*XEC*DSIN(3.0D0*RANM)               
      ANV=RANV/PIR                                                      
      TLS=ANV+XL 

c calculate for 10 degrees latitude belts
      ny=36
      Pi=3.141593d0
      do i=1,ny
       phi=(-87.5+(i-1)*5)
       RPHI=PHI*PIR                                                      
       RANV=(TLS-XL)*PIR                                                
       RAU=(1.0D0-ECC*ECC)/(1.0D0+ECC*DCOS(RANV))                        
       S1=SF/RAU/RAU                                                      
       RLAM=TLS*PIR                                                     
       SD=SO*DSIN(RLAM)                                                  
       CD=DSQRT(1.0D0-SD*SD)                                             
       RDELTA=DATAN(SD/CD)                                               
       DELTA=RDELTA/PIR                                                  
       SP=SD*DSIN(RPHI)                                                  
       CP=CD*DCOS(RPHI)                                                  
       APHI=DABS(PHI)                                                    
       ADELTA=DABS(DELTA)                                                
C                                                                       
C   SINGULARITY FOR APHI=90 AND DELTA=0                                 
C   PARTICULAR CASES FOR PHI=0  OR  DELTA=0                             
C                                                                       
      TT=DABS(APHI-90.0D0)                                              
      IF ((TT.LE.TEST).AND.(ADELTA.LE.TEST)) GO TO 22                    
      IF(ADELTA.LE.TEST) GO TO 26                                        
      IF(APHI.LE.TEST) GO TO 27                                          
C                                                                       
C   LABEL 2 : POLAR CONTINUAL NIGHT OR W=0  DAYL=0                      
C   LABEL 4 : POLAR CONTINUAL DAY                                       
C   LABEL 3 : DAILY SUNRISE AND SUNSET                                  
C   LABEL 6 : EQUINOXES                                                 
C   LABEL 7 : EQUATOR                                                   
C                                                                       
      ATT=90.0D0-ADELTA                                                  
      SPD=PHI*DELTA                                                     
      IF (APHI.LE.ATT) GO TO 23                                           
      IF (SPD) 22,23,24                                                    
 22   DAYL=0.00D0                                                       
      WW=0.00D0                                                         
      GO TO 35                                                           
 24   DAYL=24.00D0                                                      
      WW=S1*SP*PI                                                        
      GO TO 35                                                           
 23   TP=-SP/CP                                                         
      STP=DSQRT(1.0D0-TP*TP)                                            
      RDAYL=dACOS(TP)                                                  
      DAYL=24.0D0*RDAYL/PI                                              
      WW=S1*(RDAYL*SP+CP*STP)                                            
      GO TO 35                                                           
 26   DAYL=12.0D0                                                       
      WW=S1*DCOS(RPHI)                                                   
      GO TO 35                                                           
 27   DAYL=12.0D0                                                       
      WW=S1*DCOS(RDELTA)

 35   continue

      shine(i)=ww

      enddo  ! on latitudes

c start aevraging
      do i=1,4
         solar(i)=0.0
         do j=1,9
            solar(i)=solar(i)+shine((i-1)*9+j)*belt((i-1)*9+j)
         enddo
         solar(i)=solar(i)*1000/84600/boxarea(i)
      enddo


      return

      END    ! of program

                                                           
                                                            
