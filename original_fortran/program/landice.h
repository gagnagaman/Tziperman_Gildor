
      real*8 Igama,TLice1,Icoef1,LIf,IceHc,TLice2,rhoLice,gamaLI

      real*8 IceHs,ILenS(3),DhS,DtS,Im(4),HaiS,IceVolS(3),ISmax,
     &       InewTemp(4),Iflow(4),T1melt,T4melt,T12preci,T32preci,
     &       T34preci,alphaP,InewMaxS, LIdiff(4)

      real*8 IceHn,ILenN(3),DhN,DtN,HaiN,IceVolN(3),INmax,InewMaxN

      real*8 rhowLice,LIc1,LIc4,SINK_PER

      common/LiceS/IceHs,ILenS,DhS,DtS,Im,HaiS,IceVolS,ISmax,InewMaxS,
     &       InewTemp,Iflow,T1melt,T4melt,T12preci,T32preci,
     &       T34preci,alphaP,LIdiff,gamaLI

      common/LiceN/IceHn,ILenN,DhN,DtN,HaiN,IceVolN,INmax,InewMaxN

      common/Lice/Igama,TLice1,Icoef1,LIf,IceHc,TLice2,
     &       rhoLice,rhowLice,LIc1,LIc4,SINK_PER




