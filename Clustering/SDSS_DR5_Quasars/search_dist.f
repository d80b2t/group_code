c------------------------------------------------------------------------------
c Given either rcomov, z, dl or dl^2(1+z)^-beta look up the corresponding
c values of the other quantities. The integer isel specifies which argument
c is set on input.
c     
      subroutine  search_dist(rc,zplus1,dlum,dlzcomb,isel)
************************************variables*********************************
      implicit none
      integer ntable,NT,it,isel,itlo,ithi
      parameter(NT=10000)
      real*8 rcomov(NT),z(NT),dl(NT),comb(NT),rc,zplus1,dlum,dlzcomb
      real*8 drc,h,rit,zt
      save itlo,ithi
      common /dist_table/ ntable,drc,rcomov,z,dl,comb
      data itlo/1/
      data ithi/1/
******************************************************************************
      if (isel.eq.1) then !rc set so lookup zplus1,dlum and dlum^2(1+z)^-beta
c      This is the fast/easy case as the table is equally spaced in rc
         rit=1.0+rc/drc
         it = int(rit)
         h= rit-real(it) 
         if (it.ge.NT) stop 'search_dist() r beyond tabulated range'
         zplus1=1.0+ z(it)*(1.0-h)    +    z(it+1)*h
         dlum=       dl(it)*(1.0-h)  +    dl(it+1)*h
         dlzcomb=  comb(it)*(1.0-h)  +  comb(it+1)*h
c     
      else if (isel.eq.2) then  !zplus1 set 
c     Search for corresponding redshift
         zt=zplus1-1.0
         if (z(ithi).lt.zt .or. z(itlo).gt.zt ) then !short cut if zt close to
            itlo=1                                   !last call
            ithi=NT                                !otherwise do binary search
            do while (ithi-itlo.gt.1) 
               it=(ithi+itlo)/2
               if(z(it).gt.zt)then
                  ithi=it
               else
                  itlo=it
               endif
            end do      
         end if
         h=(zt-z(itlo))/(z(ithi)-z(itlo))
         rc=     rcomov(itlo)*(1.0-h)  +rcomov(ithi)*h
         dlum=       dl(itlo)*(1.0-h)  +    dl(ithi)*h
         dlzcomb=  comb(itlo)*(1.0-h)  +  comb(ithi)*h
c
      else if (isel.eq.3) then  !luminosity distance set 
c      Search for corresponding luminosity distance
c         write(*,*)'Gth #e'
c         write(*,*)'dl(ithi)',dl(ithi),'dl(itlo)',dl(itlo),'dlum', dlum
         if (dl(ithi).lt.dlum .or. dl(itlo).gt.dlum ) then
            itlo=1                                  
            ithi=NT                                 
            do while (ithi-itlo.gt.1) 
               it=(ithi+itlo)/2
               if(dl(it).gt.dlum)then
                  ithi=it
               else
                  itlo=it
               endif
            end do      
         end if
         h=(dlum-dl(itlo))/(dl(ithi)-dl(itlo))
         rc=     rcomov(itlo)*(1.0-h)  +rcomov(ithi)*h
         zplus1=1.0+ z(itlo)*(1.0-h)    +    z(ithi)*h
         dlzcomb=  comb(itlo)*(1.0-h)  +  comb(ithi)*h

c
      else if (isel.eq.4) then  !dl^2 zplus1^-beta set
c      Search for corresponding  dlum^2 zplus1^-beta set 
         if (comb(ithi).lt.dlzcomb .or. comb(itlo).gt.dlzcomb ) then
            itlo=1                                  
            ithi=NT                                 
            do while (ithi-itlo.gt.1) 
               it=(ithi+itlo)/2
               if(comb(it).gt.dlzcomb)then
                  ithi=it
               else
                  itlo=it
               endif
            end do      
         end if
         h=(dlzcomb-comb(itlo))/(comb(ithi)-comb(itlo))
         rc=     rcomov(itlo)*(1.0-h)  +rcomov(ithi)*h
         zplus1=1.0+ z(itlo)*(1.0-h)    +    z(ithi)*h
         dlum=       dl(itlo)*(1.0-h)  +    dl(ithi)*h
c
      end if
      return
      end
c-----------------------------------------------------------------------------

