      real*8 radius,belt(36),pi2rr,boxarea(4)
      real*8 PI,PIR,PIRR,STEP,TEST,y,z,xod,xes,xec,xob,xee,
     &       xse,xop, prm,t1,ss,ecc,arg,tra,rp,phi,perh,prg,
     &       sf,so,xl,xllp,xlam,lp1,spd,tp,stp,ranv

      real*8 rau,rlam,sd,cd,sp,cp,tt,att,aphi,rphi,tls,dayl,
     &       time_fa1,pre,lp2,dlam,tau,anm,ranm,anv,
     &       DLAMM,s1,rdelta,delta,adelta,ww,rdayl,shine(36)

      real*8 AE(19),BE(19),CE(19),AOB(47),BOB(47),COB(47),AOP(78)
     &       ,BOP(78),COP(78)

      integer NJM(12),NJM1(12),ny,nef,nob,nop,NEFF,NOBB,NOPP,nd



	
      common /inso1/radius,belt,pi2rr,boxarea
      common /inso2/PI,PIR,PIRR,STEP,TEST,y,z,xod,xes,xec,xob,xee,
     &       xse,xop, prm,t1,ss,ecc,arg,tra,rp,phi,perh,prg,
     &       sf,so,xl,xllp,xlam,lp1,spd,tp,stp,ranv
      common /inso3/rau,rlam,sd,cd,sp,cp,tt,att,aphi,rphi,tls,dayl,
     &       time_fa1,pre,lp2,dlam,tau,anm,ranm,anv,
     &       DLAMM,s1,rdelta,delta,adelta,ww,rdayl,shine

      common /inso4/NJM,NJM1,ny,nef,nob,nop,NEFF,NOBB,NOPP,nd 
      common /inso5/AE,BE,CE,AOB,BOB,COB,AOP,BOP,COP
