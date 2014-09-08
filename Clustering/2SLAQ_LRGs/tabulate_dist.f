c------------------------------------------------------------------------------
c Tabulate redshift, z, luminosity distance, d_L, and the combination
c d_L^2 (1+z)^-beta against comoving distance, given Omega_0 and Lambda_0
c     
      subroutine  tabulate_dist(omega0,lambda0,beta)
*************************************variables*********************************
      implicit none
      integer ntable,NT,it,NNT,jt
      real*8 EPS,beta
c      parameter(NT=1000,EPS=1.0e-04,NNT=10000)
      parameter(NT=10000,EPS=1.0e-04,NNT=10000)
      real*8 rcomov(NT),z(NT),dl(NT),comb(NT),zplus1,x,om,atanhom,rchar
      real*8 drc,zplus1t(NNT),rct(NNT),dz,intp,intt,h
      real*8 omega0,lambda0,rcuml
      common /dist_table/ ntable,drc,rcomov,z,dl,comb
*******************************************************************************
      
      ntable=NT
      if (omega0.lt.0.0) then   !flags non-relativistic assumption
      else if (abs(omega0-1.0).le.EPS .and. (lambda0).le.EPS) then !omega=1
      else if (abs(omega0-1.0).gt.EPS .and. (lambda0).le.EPS) then !open
         om=sqrt(1.0-omega0)    !useful numbers in open model
         atanhom=0.5*log((1.0+om)/(1.0-om) )
      else if (abs(omega0+lambda0-1.0).lt.EPS ) then !flat with lambda
c     First tabulate rcomov as with redshift by integrating
         dz=10.0/real(NNT)
         rcuml=0.0
         rct(1)=0.0
         zplus1t(1)=1.0
         intp=0.0
         do it=2,NNT
            zplus1t(it)=1.0+real(it-1)*dz
            intt=1.0/sqrt(omega0*zplus1t(it)**3 + 1.0 - omega0)
            rcuml=rcuml+0.5*(intt+intp)*dz
            rct(it)=rcuml*3000.0
            intp=intt
         end do
      else
         stop 'Not programmed for this cosmology'
      end if
c
c     Loop over comoving distances and compute corresponding
c     redshift and luminosity distances...
      drc=6000.0/real(NT)
      do it=1,NT   !make table uniformly spaced in comoving distance.
         rcomov(it)=real(it-1)*drc
c
         if (omega0.lt.0.0) then
            zplus1=1.0 +rcomov(it)/3000.0
            z(it)= rcomov(it)/3000.0
            dl(it)=rcomov(it)
            comb(it)= dl(it)**2 * zplus1**(-beta)
         else if (abs(omega0-1.0).le.EPS .and. (lambda0).le.EPS) then !omega=1
            zplus1=1.0/(1.0-rcomov(it)/6000.0)**2
            z(it)=zplus1-1.0
            dl(it)=rcomov(it)*zplus1
            comb(it)=dl(it)**2 * zplus1**(-beta)
c
         else if (abs(omega0-1.0).gt.EPS .and. (lambda0).le.EPS) then !open
            x=tanh(atanhom-om*rcomov(it)/6000.0)
            zplus1=((om/x)**2-om**2)/omega0
            z(it)=zplus1-1.0
            rchar=(6000.0/(zplus1*omega0**2)) *
     &           (omega0*z(it)+(omega0-2.0)*(sqrt(1+omega0*z(it))-1.0))
            dl(it)=rchar*zplus1
            comb(it)=dl(it)**2 * zplus1**(-beta)
c
         else if (abs(omega0+lambda0-1.0).lt.EPS ) then !flat
c           look up redshift from temporary table
            jt=1
            do while (rcomov(it).ge.rct(jt)) 
               jt=jt+1
               if (jt.gt.NNT)   stop
     &         'tabulate_dist(): rct() not tabulated to sufficient z'
            end do
            h=(rcomov(it)-rct(jt-1))/(rct(jt)-rct(jt-1))
            zplus1=zplus1t(jt-1)*(1.0-h) + zplus1t(jt)*h
            z(it)=zplus1-1.0
            dl(it)=rcomov(it)*zplus1
            comb(it)=dl(it)**2 * zplus1**(-beta)
c     
         end if
      end do
      
      return
      end





