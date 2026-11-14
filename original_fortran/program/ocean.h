      real*8 D(8), Owidth(8), r, gammaS,gammaT,vol(8),recvol(8),
     &       recL(8),recD(8),totheat,Oarea(8),totsalt,totvolume,
     & lowersection,uppersection,horisection,Kv(4),Kh(6),Hlho(4)

      real*8 q21,q62,q23,q73,q87,q34,q76,q15,q48,q56,
     $       t(8,3),s(8,3),rho(8),sflux(4),tflux(4),
     $       c21,c23,c34,c56,c76,c87,fsi(4,3)

      real*8 q21s,q62s,q23s,q73s,q87s,q34s,q76s,q15s,q56s,
     $       q21t,q62t,q23t,q73t,q87t,q34t,q76t,q15t,q56t,
     $       rhos(8),rhot(8),tc(8),sc(8)

      real*8 sf21,sf62,sf23,sf73,sf87,   
     $       sf34,sf76,sf15,sf48,sf56,
     $       tf21,tf62,tf23,tf73,tf87,   
     $       tf34,tf76,tf15,tf48,tf56

      real*8 sdf15,sdf51,sdf12,sdf21,sdf56,sdf65,
     $       sdf62,sdf26,sdf23,sdf32,sdf67,sdf76,
     $       sdf73,sdf37,sdf34,sdf43,sdf48,sdf84,
     $       sdf87,sdf78,
     $       tdf15,tdf51,tdf12,tdf21,tdf56,tdf65,
     $       tdf62,tdf26,tdf23,tdf32,tdf67,tdf76,
     $       tdf73,tdf37,tdf34,tdf43,tdf48,tdf84,
     $       tdf87,tdf78

      real*8 rhoIw,rhoLIw,dsigma,dsigma2

      COMMON/OVAR/q21,q62,q23,q73,q87,q34,q76,q15,q48,q56,
     $    sf21,sf62,sf23,sf73,sf87,sf34,sf76,sf15,sf48,sf56,
     $    tf21,tf62,tf23,tf73,tf87,tf34,tf76,tf15,tf48,tf56,
     $    t,s,rho,sflux,tflux,tc,sc

      COMMON/OVAR2/sdf15,sdf51,sdf12,sdf21,sdf56,sdf65,
     $       sdf62,sdf26,sdf23,sdf32,sdf67,sdf76,
     $       sdf73,sdf37,sdf34,sdf43,sdf48,sdf84,
     $       sdf87,sdf78,
     $       tdf15,tdf51,tdf12,tdf21,tdf56,tdf65,
     $       tdf62,tdf26,tdf23,tdf32,tdf67,tdf76,
     $       tdf73,tdf37,tdf34,tdf43,tdf48,tdf84,
     $       tdf87,tdf78,fsi


      common/Oparameters1/D, Owidth, r, gammaS,gammaT,
     &       vol,recvol,c21,c23,c34,c56,c76,c87,
     &       recL,recD,totheat,Oarea,totsalt,totvolume,
     &       lowersection,uppersection,horisection,Kv,Kh

      common/Ovar3/q21s,q62s,q23s,q73s,q87s,q34s,q76s,q15s,q56s,
     $       q21t,q62t,q23t,q73t,q87t,q34t,q76t,q15t,q56t,
     $       rhos,rhot,Hlho,rhoIw,rhoLIw,dsigma,dsigma2
