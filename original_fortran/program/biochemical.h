      real*8  BT(8),CT(8,3),AT(8,3),CO2(8),HCO3(8),CO3(8),BOH3(8),
     $        BOH4(8),CTflux(4),CO2s(4),pv(4),CT0,AT0,a(8),PO4(8,3),
     $        O2(8,3),PROD(8),PRODC(8),ATprod(8),Rain(8),CaCo3C(8),
     $        ATcaco3(8),O2flux(8),PO4flux(8),pv0

      real*8  CTf21,CTf15,CTf56,CTf34,CTf48,CTf87,CTf62,CTf23,
     $        CTf73,CTf76,CTdf15,CTdf26,CTdf37,CTdf48,CTdf12,
     $        CTdf23,CTdf34,CTdf56,CTdf76,CTdf78

      real*8  ATf21,ATf15,ATf56,ATf34,ATf48,ATf87,ATf62,ATf23,
     $        ATf73,ATf76,ATdf15,ATdf26,ATdf37,ATdf48,ATdf12,
     $        ATdf23,ATdf34,ATdf56,ATdf76,ATdf78

      real*8  PO4f21,PO4f15,PO4f56,PO4f34,PO4f48,PO4f87,PO4f62,PO4f23,
     $        PO4f73,PO4f76,PO4df15,PO4df26,PO4df37,PO4df48,PO4df12,
     $        PO4df23,PO4df34,PO4df56,PO4df76,PO4df78

      real*8  O215,O256, O248,O287,O262,O273,O276,O2df15,O2df26,
     $        O2df37,O2df48,O2df56,O2df76,O2df78

      real*8 ct1,ct2,a_temp1,a_temp2,PREC,AC_temp,gama_t,pc1(4),
     $       RedC(8),RedN,RedO2,h,PO40,O20,totPROD,PO4av,ATav,
     $       sProd,hProd, RedC0, pc10(4)


      common /biochem1/BT,CT,AT,CO2,HCO3,CO3,BOH3,BOH4,CTflux,
     $                 CO2s,pv,CT0,AT0,a,PO4,O2,PROD,PRODC,ATprod,
     $                 Rain,CaCo3C,ATcaco3,O2flux,PO4flux,pv0 

      common /biochem2/CTf21,CTf15,CTf56,CTf34,CTf48,CTf87,CTf62,
     $        CTf23,CTf73,CTf76,CTdf15,CTdf26,CTdf37,CTdf48,CTdf12,
     $        CTdf23,CTdf34,CTdf56,CTdf76,CTdf78

      common /biochem3/ATf21,ATf15,ATf56,ATf34,ATf48,ATf87,ATf62,
     $        ATf23,ATf73,ATf76,ATdf15,ATdf26,ATdf37,ATdf48,ATdf12,
     $        ATdf23,ATdf34,ATdf56,ATdf76,ATdf78

      common /biochem4/PO4f21,PO4f15,PO4f56,PO4f34,PO4f48,PO4f87,PO4f62,
     $        PO4f23,PO4f73,PO4f76,PO4df15,PO4df26,PO4df37,PO4df48,
     $        PO4df12,PO4df23,PO4df34,PO4df56,PO4df76,PO4df78

      common /biochem5/O215,O256, O248,O287,O262,O273,O276,O2df15,
     $                 O2df26,O2df37,O2df48,O2df56,O2df76,O2df78


      common /biochem4/ct1,ct2,a_temp1,a_temp2,PREC,AC_temp,PO4av,
     $        gama_t,pc1,RedC,RedN,RedO2,h,PO40,O20,totPROD,ATav,
     $       sProd,hProd,redc0,pc10
