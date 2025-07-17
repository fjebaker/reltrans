!-----------------------------------------------------------------------
      subroutine getcont(Cp, earx, nex, Gamma, Ecut_obs, logxi, logne, contx)
!!! Calculates continuum spectrum calling nthComp with the correct normalisation
!!!based on the xillver spectrum 
!!!  Arg:
        !  earx: energy grid
        !  nex: number of grid points
        !  Gamma: continuum spectrum inclination
        !  Ecut_obs: high energy cut-off or electron temperature
        !  logxi: ionisation parameter
        !  logne: density 
        !  (output) contx: continuum spectrum
!!! Last change: Gullo 2023 Nov; adapted to match the xillver tables call 

      use gr_continuum
      implicit none
      integer, intent(in)           :: nex, Cp
      real   , intent(in)           :: earx(0:nex), Ecut_obs, logxi, logne
      real   , intent(out)          :: contx(nex)
      double precision , intent(in) :: Gamma

      real   , parameter  :: pi = acos(-1.0),ergsev  = 1.602197e-9 ! Convert keV to ergs
      integer :: i, ifl
      real    :: nth_par(5), photer(nex), E, Icomp, inc_flux
      real    :: get_norm_cont_local

      Icomp = 0.0

      if (Cp .eq. 2) then
!So far this works only with kTe, so only with nthComp continuum model
         nth_par(1) = real(Gamma)
         nth_par(2) = Ecut_obs
         nth_par(3) = 0.05
         nth_par(4) = 1.0
         ! nth_par(5) = 0.0
         nth_par(5) = (1.0/ real(gso(1))) - 1.0
         Ifl=1

      ! write(*,*) 'continuum parameters', nth_par
         call donthcomp(earx, nex, nth_par, ifl, contx, photer)
!the continuum needs to be renormalised according to the illuminating flux that was considered in xillver 
! Plus we divide by a factor that depends on ionisation and density to agree with the first versions of reltrans

         do i = 1, nex
            E   = 0.5 * ( earx(i) + earx(i-1) )
            if (E .ge. 0.1 .and. E .le. 1e3) then
               Icomp = Icomp + ((earx(i) + earx(i-1)) * 0.5 * contx(i))
            endif
         enddo
         inc_flux = 10**(logne + logxi) / (4.0 * pi) / ergsev !calculate incident flux in units  [keV/cm^2/s]
         get_norm_cont_local = inc_flux/ Icomp / 1e20
         
         ! contx = contx * get_norm_cont(real(Gamma), Ecut_obs, logxi, logne)
         ! contx = contx  / 10**(logxi + logne - 15)
         contx = contx * get_norm_cont_local / (10**(logxi + logne - 15))
         ! write(*,*) 'continuum normalization parameters', get_norm_cont_local, logxi, logne

      ! endif
         
      else

         ! write(*,*) 'continuum parameters ', Cp, Gamma, Ecut_obs
         do i = 1, nex
            E   = 0.5 * ( earx(i) + earx(i-1) )
            contx(i) = E**(-1.0*real(Gamma)+1) * exp(-E/(Ecut_obs)) 
            if (E .ge. 0.1 .and. E .le. 1e3) then
               Icomp = Icomp + ((earx(i) + earx(i-1)) * 0.5 * contx(i))
            endif
         enddo
         inc_flux = 10**(logne + logxi) / (4.0 * pi) / ergsev !calculate incident flux in units  [keV/cm^2/s]
         get_norm_cont_local = inc_flux/ Icomp / 1e20
         contx = contx * get_norm_cont_local / (10**(logxi + logne - 15))
      endif
      ! do i = 1, nex
      !    write(20, *) (earx(i-1) + earx(i)) *0.5, contx(i)
      ! enddo
      ! write(20, *) 'no no'
         
      return
    end subroutine getcont
!-----------------------------------------------------------------------
