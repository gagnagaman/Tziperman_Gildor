
      real*8 tetam(4,12),Haom(4,12),Houtm(4,12),Hinm(4,12),
     $       tfluxm(4,12),sfluxm(4,12),FWAOm(4,12),qm(4,12),
     $   albedom(4,12),solarm(4,12),Havailm(4,12),ILnewm(4,12)

      real*8 tm(8,12),sm(8,12),rhom(8,12),LIdiffm(4,12),
     $     LSIdiffm(4,12),rhosm(8,12),rhotm(8,12),vmixm(4,12),
     $     ctf21m(12),ctf15m(12),ctf56m(12),ctf34m(12),ctf48m(12),
     $     ctf87m(12),ctf62m(12),ctf73m(12),po4f21m(12),po4f15m(12),
     $     po4f56m(12),po4f34m(12),po4f48m(12),po4f87m(12),
     $     po4f62m(12),po4f73m(12)

      real*8 ctdf15m(12),ctdf26m(12),ctdf48m(12),ctdf12m(12),
     $     ctdf34m(12),po4df26m(12),po4df48m(12),po4df12m(12),
     $     po4df34m(12),caco3m(4,12),ctdf56m(12),ctdf78m(12),
     $ po4df78m(12),po4df56m(12),po4df15m(12)

      
      real*8 CTm(8,12),ATm(8,12),HCO3m(8,12),CO3m(8,12),
     $  CTfluxm(4,12),am(8,12),co2am(12),CTavm(12),PO4m(8,12),
     $  Rainm(8,12),Prodm(8,12),PO4avm(12),ATavm(12),sprodm(12),
     $  hprodm(12),totprodm(12),CO2m(8,12),CO2aam(8,12)



      real*8 q21m(12),q62m(12),q23m(12),q73m(12),q87m(12),q34m(12),
     $       q76m(12),q15m(12),q48m(12),q56m(12),Savm(12),Tavm(12),
     $ tetaavm(12),Imm(4,12),Iflowm(4,12),Hneedm(4,12),SIvolm(4,12),
     $       IceHsm(12),IceHnm(12),IvolSm(12),IvolNm(12),FT21m(12),
     $       FT23m(12),FT34m(12),FW21m(12),FW23m(12),FW34m(12),
     $   q21sm(12),q21tm(12),q23sm(12),q23tm(12),q34sm(12),q34tm(12)


      real*8 fsim(4,12),Idepthm(4,12),fim(4,12),Imeltm(4,12),
     $    Ipnewm(4,12),Inewfm(4,12),Inewm(4,12),tcm(8,12),scm(8,12),
     $    inewn2,FW1m(12),FW4m(12),hlhm(4,12),tf21m(12),tf15m(12),
     $    tdf12m(12), tdf15m(12), RedCm(4,12)

c      real*8 NJMS(12),NJM1S(12),NJM2S(12)
       integer NJMS(12),NJM1S(12),NJM2S(12)
      
      common/aver1/tetam,Haom,Houtm,Hinm,tfluxm,sfluxm,FWAOm,qm,
     $             albedom,solarm,tm,sm,rhom,fsim,Idepthm,fim,
     $             Imeltm,Ipnewm,Inewfm,hlhm,NJMS,NJM1S,NJM2S

      common/aver2/q21m,q62m,q23m,q73m,q87m,q34m,q76m,q15m,q48m,
     $     q56m,Savm,Tavm,tetaavm,Imm,Inewm,Iflowm,Havailm,Lidiffm,
     $     IceHsm,IceHnm,IvolSm,IvolNm,FT21m,FT23m,FT34m,FW21m,ILnewm,
     $     FW23m,FW34m,q21sm,q21tm,q23sm,q23tm,q34sm,q34tm,Hneedm,
     $     SIvolm,LSIdiffm,rhosm,rhotm,tcm,scm,inewn2,FW1m,FW4m

      common/aver3/CTm,ATm,HCO3m,CO3m,CO2m,CTfluxm,am,co2am,CTavm,
     $       PO4m,Rainm,Prodm,PO4avm,ATavm,sprodm,hprodm,totprodm,
     $       tf21m,tf15m,tdf12m,tdf15m,CO2aam,vmixm,RedCm

      common/aver4/ctdf15m,ctdf26m,ctdf48m ,ctdf12m ,
     $     ctdf34m ,po4df26m ,po4df48m ,po4df12m ,
     $     po4df34m ,caco3m,ctf21m ,ctf15m ,ctf56m ,
     $     ctf34m ,ctf48m ,ctf87m ,ctf62m ,ctf73m ,
     $     po4f21m ,po4f15m ,po4f56m ,po4f34m ,po4f48m ,
     $     po4f87m ,po4f62m ,po4f73m ,ctdf56m,ctdf78m,
     $     po4df78m,po4df56m,po4df15m
