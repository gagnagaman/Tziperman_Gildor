c=======================================================================        
                real*8 function UNESCO(sin, tin)                                  
c=======================================================================        
c     this subroutine calculates the density of seawater using the              
c     standard equation of state recommended by unesco(1981).                   
c                                                                               
c     input [units]:                                                            
c       in-situ temperature (t): [degrees centigrade]                           
c       salinity (s): [practical salinity units]                                
c       pressure (pin): [decibars, approx. as meters of depth]                  
c     output [units]:                                                           
c       density(rho): kilograms per cubic meter                                 
c         or                                                                    
c       sigma(sigma):sigma units                                                
c                                                                               
c     references:                                                               
c        Gill, A., Atmosphere-Ocean Dynamics: International Geophysical         
c         Series No. 30. Academic Press, London, 1982, pp 599-600.              
c        UNESCO, 10th report of the joint panel on oceanographic tables         
c          and standards. UNESCO Tech. Papers in Marine Sci. No. 36,            
c          Paris, 1981.                                                         
c                                                                               
c-----------------------------------------------------------------------        
c                                                                               
      real*8 sin, tin, pin
      doubleprecision t,s,p,rho,c1p5                                            
     $,               rw,rsto,xkw,xksto,xkstp                                   
      c1p5 = 1.5d0                                                              
c                                                                               
c-----------------------------------------------------------------------        
c    


c hezi 14.3.99

      pin=0.0
c  temperature [degrees centigrade]    
      if (tin.eq.-100) then
         unesco=-100.0
         else
      t = dble(tin)  
c  salinity [psu] 
      s = dble(sin) 
c  convert from depth [m] (decibars) to bars                                    
      p = dble(pin*0.1)                                                         
c        
      rw =     9.99842594d2 + 6.793952d-2*t - 9.095290d-3*t**2                  
     $        + 1.001685d-4*t**3 - 1.120083d-6*t**4 + 6.536332d-9*t**5          
c                                                                               
      rsto =   rw + (8.24493d-1 - 4.0899d-3*t + 7.6438d-5*t**2                  
     $        - 8.2467d-7*t**3 + 5.3875d-9*t**4) * s                            
     $       + (-5.72466d-3 + 1.0227d-4*t - 1.6546d-6*t**2) * s**c1p5           
     $       + 4.8314d-4 * s**2                                                 
c                                                                               
      xkw =    1.965221d4 + 1.484206d2*t - 2.327105d0*t**2 +                    
     $         1.360477d-2*t**3 - 5.155288d-5*t**4                              
c                                                                               
      xksto =  xkw + (5.46746d1 - 6.03459d-1*t + 1.09987d-2*t**2                
     $        - 6.1670d-5*t**3) * s                                             
     $       + (7.944d-2 + 1.6483d-2*t - 5.3009d-4*t**2) * s**c1p5              
c                                                                               
      xkstp =  xksto + (3.239908d0 + 1.43713d-3*t + 1.16092d-4*t**2             
     $        - 5.77905d-7*t**3) * p                                            
     $       + (2.2838d-3 - 1.0981d-5*t - 1.6078d-6*t**2) * p * s               
     $       + 1.91075d-4 * p * s**c1p5                                         
     $       + (8.50935d-5 - 6.12293d-6*t + 5.2787d-8*t**2) * p**2              
     $       + (-9.9348d-7 + 2.0816d-8*t + 9.1697d-10*t**2) * p**2 * s          
c                                                                               
      rho   =  rsto / (1.0d0 - p/xkstp)                                         
c     sigma =  rho                                                              
      unesco =  rho - 1.0d3 
      endif
      return                                                                    
      end       
