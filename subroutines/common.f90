module common_types
    implicit none

    ! this should be moved so that all of the wrappers pass the following
    ! instead of a params array
    type :: t_model_arguments
        ! These are hardcoded with a fixed size of 2, since nlp <= 2. by not
        ! having them be pointers or dynamically allocated, it makes reasoning
        ! about the code a little easier.
        double precision :: h(2)
        real :: DelAB(2), g(2)
        ! Black hole spin (dimensionless):
        double precision :: a
        ! Observer inclination in degrees:
        double precision :: inc
        ! Inner and outer radii of the accretion disc (rg):
        double precision :: rin, rout
        ! Scale height of the accretion disc
        double precision :: honr
        ! Cosmological redshift:
        double precision :: zcos
        ! Photon index
        double precision :: Gamma
        real :: logxi, Afe, lognep, Cutoff_obs, Cutoff_s, Dkpc, Anorm, beta_p
        real :: Nh, boost, Mass, floHz, fhiHz, DelA
        integer :: nlp, ReIm, resp_matr, Cp

        ! eta: Fourier frequency dependent normalisation ratio C1(vc) / C2(vc)
        ! See Lucchini et al. 2023 for details.
        double precision :: eta

        ! eta_0: Time-averaged normalisation ratio C1 / C2 between the continuum
        ! of the two lampposts. It sets the continuum cutoff and disc
        ! ionisation.
        ! See Lucchini et al. 2023 for details.
        double precision :: eta_0

        ! Asymmetry parameter of angular emissivity function:
        double precision :: qboost

        ! Linear coefficient of angular emissivity function (b1) and quadratic
        ! coefficient of angular emissivity function (b2):
        double precision :: b1, b2

        ! The below are computed from the above
        ! The cosine angle
        double precision :: muobs
    end type t_model_arguments

    type :: t_cached_parameters
        type(t_model_arguments) :: args
        ! Previous number of frequency bins, and low and high frequency:
        integer :: nf = -1
        double precision :: flo, fhi
    end type t_cached_parameters

    type :: t_config
        integer :: verbose = 0
        ! firstcall: is this the first time the model has been called?
        logical :: firstcall = .true., needtrans = .true., needconv = .true.,test = .false.
        ! Disable the parameter caching. If this is true, then `needs_check`
        ! will always configure everything to be recalculated regardless of the
        ! state of the cache.
        logical :: disable_cache = .false.

        ! me: Number of mue bins
        ! xe: Number of logr bins: bins 1:xe-1 are logarithmically spaced, bin
        ! xe is everything else
        integer :: me, xe, nex, m, ionvar, refvar

        ! TODO: are these really constants, or are they hidden variables?
        ! Constants

        ! Relativistic grid size, i.e. pixel count splitting up the image plane
        ! in polar coordinates.
        integer :: nphi = 200, nro = 200
        ! Non-relativistic grid size:
        integer :: nron = 100, nphin = 100

        ! Observer-to-black-hole distance, set in `do_first_call_setup`
        double precision :: distance

        ! The number of time bins used for calculating the impulse response function:
        integer :: nt = 2**9

        ! Emin, Emax: min and max of the internal energy grid
        ! (different from output grid)
        real :: Emin = 1e-2, Emax = 3e3, dyn = 0.0

        ! rnmax: max radius for which GR ray-tracing is used
        ! dlogf: a resolution parameter in base 10
        double precision :: rnmax = 300.d0, dlogf = 0.09

        real :: DeltaGamma = 0.01

        ! Use ring-like coronal model. Should future models get added, this
        ! could be promoted to an enumeration of some description.
        logical :: ring_like = .false.

        ! Toggle whether to calculate the impulse response or not. Since the
        ! impulse response is not directly used in calculations, and must be
        ! read out in another way, this defaults to false.
        logical :: calculate_impulse_response = .false.

        ! Internal frequency grid. These values are all set in the
        ! `config_frequency` subrtouine.
        ! Number of frequency bins:
        integer :: nf
        real :: f, fac
        ! Frequency center, low, and high:
        double precision :: fc, flo, fhi
        ! internal frequency grid, for when we do lag/frequency spectra
        integer :: fbinx
        real, allocatable :: fix(:)

        ! relativistic parameters and limit on rin and h
        double precision :: rmin, rh
        double precision, allocatable :: height(:)

        ! internal energy grid (nex) and output/xspec (ne) energy grid
        ! dloge: logarithmic resolution of the internal energy grid
        real :: E, dE, dloge

        ! Radial and angle profile
        integer :: mubin, rbin, ibin

        ! variable for non linear effects
        integer :: DC, ionvariation
        real :: dlogxi1, dlogxi2

        ! Parameters that are cached and saved between successive runs of the model.
        type(t_cached_parameters) :: cached

     end type t_config

    type :: t_arrays
        ! earx: internal energy grid array (0:nex)
        real, allocatable :: earx(:), ear(:), fix(:)
        real, allocatable :: ReGbar(:), ImGbar(:)
        real, allocatable :: contx(:,:)
        double precision, allocatable :: contx_int(:)
        ! TRANSFER FUNCTIONS and Cross spectrum dynamic allocation + variables
        complex, dimension(:,:,:,:,:), allocatable :: ker_W0, ker_W1, ker_W2, ker_W3
        ! ker_W0(nlp,ne,nf,me,xe) Transfer function W0 - linear transfer function
        ! ker_W1(nlp,ne,nf,me,xe) Transfer function W1 - one aspect of photon index variations
        ! ker_W2(nlp,ne,nf,me,xe) Transfer function W2 - other aspect of photon index variations
        ! ker_W3(nlp,ne,nf,me,xe) Transfer function W3 - ionization variations
        real, dimension(:,:,:), allocatable :: ReW0, ImW0, ReW1, ImW1
        real, dimension(:,:,:), allocatable :: ReW2, ImW2, ReW3, ImW3
        real, dimension(:,:), allocatable :: ReSraw, ImSraw, ReSrawa, ImSrawa, ReGrawa, ImGrawa, ReG, ImG

        ! Observed fraction, reflection fraction, lensing factors.
        double precision, dimension(:), allocatable :: frobs, frrel, lens
    end type t_arrays

    type(t_config), target, save :: global_config
contains

    subroutine reset_reltrans() bind(C, name="reset_reltrans")
      !The routine reset the check variables, in order to force
      ! a fresh start of the model
      global_config%firstcall = .true.
    end subroutine reset_reltrans


    ! Unwraps the arguments from a parameter array into `args`.
    subroutine unwrap_arguments(args, nlp, dset, params, cutoff_powerlaw)
        use rtconstants, only: parse_reim
        double precision, parameter :: pi = acos(-1.d0)
        integer, intent(in) :: nlp, dset, cutoff_powerlaw
        real, target, intent(in) :: params(32)
        type(t_model_arguments), intent(out) :: args
        integer :: i
        do i = 1,nlp
            args%DelAB(i) = params(27 + (i - 1) * nlp)
            args%g(i) = params(28 + (i - 1) * nlp)
        end do
        if (dset .eq. 1) then
           args%Dkpc = params(9)
           args%logxi = 0.0
        else
           args%logxi = params(9)
        end if
        args%h(1) = dble(params(1))
        args%h(2) = dble(params(2))
        args%nlp = nlp
        args%a = dble(params(3))
        args%inc = dble(params(4))
        args%muobs = cos(args%inc * pi / 180.d0)
        args%rin = dble(params(5))
        args%rout = dble(params(6))
        args%zcos = dble(params(7))
        args%Gamma = dble(params(8))
        args%Afe = params(10)
        args%lognep = params(11)
        args%Cutoff_s = params(12)
        args%Cutoff_obs = params(12)
        args%eta_0 = params(13)
        args%eta = params(14)
        args%beta_p = params(15)
        args%Nh = params(16)
        args%boost = params(17)
        args%qboost = dble(params(18))
        args%Mass = dble(params(19))
        args%honr = dble(params(20))
        args%b1 = dble(params(21))
        args%b2 = dble(params(22))
        args%floHz = params(23)
        args%fhiHz = params(24)
        args%ReIm = parse_reim(int(params(25)))
        args%DelA = params(26)
        args%Anorm = params(31)
        args%resp_matr = params(32)
        args%Cp = cutoff_powerlaw
    end subroutine unwrap_arguments

    ! Adjust the model parameters to sane values and set the derived values in
    ! `config`. It performs the following steps:
    ! - Checks `a`, `rin`, `h` are in bounds.
    ! - Sets the inner radius to the ISCO.
    subroutine arguments_check(config, model_args)
        use rtconstants, only: MODE_CROSS_SPEC_REAL_REF_FOLDED
        type(t_config), intent(inout) :: config
        type(t_model_arguments), intent(inout) :: model_args
        integer :: i
        double precision :: disco

        ! TODO: should this be printing if it modifies the input parameters?
        ! some kind of warning perhaps?
        if (abs(model_args%a) .gt. 0.999) then
            model_args%a = sign(model_args%a,1.d0) * 0.999
        end if
        config%rmin = disco(model_args%a)
        config%rh = 1.d0+sqrt(1.d0-model_args%a**2)
        if (model_args%rin .lt. 0.d0) then
            model_args%rin = abs(model_args%rin) * config%rmin
        end if
        if (model_args%rin .lt. config%rmin)then
            write(*,*)"Warning! rin<ISCO! Set to ISCO"
            model_args%rin = config%rmin
        end if
        do i=1,model_args%nlp
            if (model_args%h(i) .lt. 0.d0) then
                model_args%h(i) = abs(model_args%h(i)) * config%rh
            end if
            if (model_args%h(i) .lt. 1.5d0*config%rh)then
                write(*,*)"Warning! h<1.5*rh! Set to 1.5*rh"
                model_args%h(i) = 1.5d0 * config%rh
            end if
        end do

        ! Decide if this is the DC component/time averaged spectrum or not
        if (config%flo .lt. tiny(config%flo) .or. config%fhi .lt. tiny(config%fhi))then
            config%DC = 1
            model_args%g = 0.0
            model_args%DelAB = 0.0
            model_args%DelA = 0.0
            model_args%ReIm = MODE_CROSS_SPEC_REAL_REF_FOLDED
            model_args%eta = model_args%eta_0
            ! this is an ugly hack for the double LP model to calculate the time-
            ! averaged spectrum
            model_args%beta_p = 1.
        else
            config%DC = 0
            model_args%boost = abs(model_args%boost)
        end if

        ! Determine what needs to be recalculated.
        call need_check(config, model_args)
    end subroutine arguments_check

    subroutine need_check(config, model_args)
        !> Checks if reltrans needs to calculate the kernel. Updated fields in
        !> `config` with the outcome of the checks.
    ! Parameters that rtrans() is sensitive to:
    ! (1-9):   h1,h2,a,inc,rin,rout,zcos,Gamma,logxi/Dkpc
    ! (11):    lognep
    ! (13):    eta_0
    ! (18):    qboost
    ! (20-22): honr,b1,b2
    ! (31):    Anorm
    ! Also need to check if the frequency range changes
    !
    ! Parameters the restframe spec is sensitive to
    ! (10):     Afe
    ! (12):    Ecut/kTe

    !!! Arg:
      ! INPUTS
      !   Cp:        defines which model
      !   Cpsave:    saved Cp
      !   param:     parameter array
      !   paramsave: saved array
      !   fhi:       high frequency range
      !   flo:       low frequency range
      !   fhisave:   saved frequency
      !   flosave:   saved frequency
      !   nf:        number of frequency bins
      !   nfsave:    saved number
      ! OUTPUTS
      !   needtrans: if true, we must do the kernel calculation
    !> Checks if reltrans needs to calculate the kernel
    !> Parameters that rtrans() is sensitive to:
    !> (1-9):   h1,h2,a,inc,rin,rout,zcos,Gamma,logxi/Dkpc
    !> (11):    lognep
    !> (13):    eta_0
    !> (18):    qboost
    !> (20-22): honr,b1,b2
    !> (31):    Anorm
    !> Also need to check if the frequency range changes
    !>
    !> Parameters the restframe spec is sensitive to
    !> (10):     Afe
    !> (12):    Ecut/kTe
    !> Inputs:
    !>     Cp:        defines which model
    !>     Cpsave:    saved Cp
    !>     param:     parameter array
    !>     paramsave: saved array
    !>     fhi:       high frequency range
    !>     flo:       low frequency range
    !>     fhisave:   saved frequency
    !>     flosave:   saved frequency
    !>     nf:        number of frequency bins
    !>     nfsave:    saved number
    !> Outputs:
    !>     needtrans: if true, we must do the kernel calculation
    !>     neecconv:  if true, we must do the convolution
        implicit none
        type(t_config), intent(inout) :: config
        type(t_model_arguments), intent(in) :: model_args
        real            , parameter   :: tol = 1e-7
        double precision, parameter   :: dtol = 1e-5

        ! functions
        integer :: get_env_int

        if (config%disable_cache) then
            config%needtrans = .true.
            config%needconv = .true.
            return
        end if

        config%needtrans = .false.
        config%needconv = .false.

        ! First check the parameter entries:
        if (abs(model_args%h(1) - config%cached%args%h(1)) > tol) then
            config%needtrans = .true.
        end if
        if (abs(model_args%h(2) - config%cached%args%h(2)) > tol) then
            config%needtrans = .true.
        end if
        if (abs(model_args%a - config%cached%args%a) > tol) then
            config%needtrans = .true.
        end if
        if (abs(model_args%inc - config%cached%args%inc) > tol) then
            config%needtrans = .true.
        end if
        if (abs(model_args%rin - config%cached%args%rin) > tol) then
            config%needtrans = .true.
        end if
        if (abs(model_args%rout - config%cached%args%rout) > tol) then
            config%needtrans = .true.
        end if
        if (abs(model_args%zcos - config%cached%args%zcos) > tol) then
            config%needtrans = .true.
        end if
        if (abs(model_args%gamma - config%cached%args%gamma) > tol) then
            config%needtrans = .true.
        end if
        if (abs(model_args%logxi - config%cached%args%logxi) > tol) then
            config%needtrans = .true.
        end if
        if (abs(model_args%lognep - config%cached%args%lognep) > tol) then
            config%needtrans = .true.
        end if
        if (abs(model_args%eta_0 - config%cached%args%eta_0) > tol) then
            config%needtrans = .true.
        end if
        if (abs(model_args%qboost - config%cached%args%qboost) > tol) then
            config%needtrans = .true.
        end if
        if (abs(model_args%honr - config%cached%args%honr) > tol) then
            config%needtrans = .true.
        end if
        if (abs(model_args%b1 - config%cached%args%b1) > tol) then
            config%needtrans = .true.
        end if
        if (abs(model_args%b2 - config%cached%args%b2) > tol) then
            config%needtrans = .true.
        end if
        if (abs(model_args%Anorm - config%cached%args%Anorm) > tol) then
            config%needtrans = .true.
        end if

        ! Now check if frequency range and frequency grid have changed
        if( config%nf .ne. config%cached%nf ) then
            config%needtrans = .true.
        else if( abs(1.- (config%fhi / config%cached%fhi) ) > dtol) then
            config%needtrans = .true.
        else if( abs(1. - (config%flo - config%cached%flo) ) > dtol) then
            config%needtrans = .true.
        end if

        ! Now for needconv
        if( config%needtrans ) config%needconv = .true.
        if( model_args%cp .ne. config%cached%args%cp ) config%needconv = .true.
        if (abs(model_args%afe - config%cached%args%afe) > tol) then
            config%needconv = .true.
        end if
        if (abs(model_args%cutoff_obs - config%cached%args%cutoff_obs) > tol) then
            config%needconv = .true.
        end if
        if (abs(model_args%cutoff_s - config%cached%args%cutoff_s) > tol) then
            config%needconv = .true.
        end if
    end subroutine need_check

    ! Read in environment variables that configure reltrans
    subroutine read_environment_variables(config)
        use env_variables, only: adensity, idum
        use xillver_tables, only: path_tables, xillver, xillverDCp,pathname_xillver, pathname_xillverDCp
        type(t_config), intent(inout) :: config
        integer :: get_env_int, temp_test
        character (len=200) :: get_env_char
        config%me = get_env_int("MU_ZONES", 1)
        config%xe = get_env_int("ION_ZONES", 20)
        ! verbose:
        ! 0: XSPEC output only
        ! 1: Print quantities to terminal + 0
        ! 2: Model components, radial scalings, impulse responses written to
        ! file + 1
        config%verbose = get_env_int("REV_VERB", 0)
        ! include pivoting reflection
        config%refvar = get_env_int("REF_VAR", 1)
        ! include ionisation changes
        config%ionvar = get_env_int("ION_VAR", 1)

        ! Whether to disable the parameter cache system.
        if (0 .ne. get_env_int("REV_NOSAV", 0)) then
            config%disable_cache = .true.
        end if

        ! these are set in `env_variables`
        ! seed for simulation
        idum = get_env_int("SEED_SIM", -2851043)
        ! decide between zone A density profile or constant density profile
        adensity = max(min(get_env_int("A_DENSITY", 0), 1), 0)

        temp_test = get_env_int("TEST_RUN", 0)
        if (temp_test .eq. 0) then
           config%test = .false.
        else
           config%test = .true.
        end if

        ! this is from xillver_tables, sets the paths where the tables are read
        ! from
        path_tables = get_env_char("RELTRANS_TABLES", './')
        write(pathname_xillver, '(A, A, A)') trim(path_tables), '/', trim(xillver)
        write(pathname_xillverDCp, '(A, A, A)') trim(path_tables), '/', trim(xillverDCp)

        write(*,*)"----------------------------------------------------"
        write(*,*)" *** ENVIRONMENT VARIABLES *** "
        write(*,*)
        write(*,*) 'RADIAL ZONES', config%xe
        write(*,*) 'ANGLE ZONES', config%me
        if (adensity .eq. 0.0) then
            write(*,*) 'A_DENSITY:', adensity, 'Density profile is constant'
        else
            write(*,*) 'A_DENSITY:', adensity, 'Density profile is zone A SS73'
        endif
        write(*,*) 'VERBOSE is ', config%verbose
        write(*,*) 'REFVAR is ', config%refvar
        write(*,*) 'IONVAR is ', config%ionvar
        write(*,*)"----------------------------------------------------"

      end subroutine read_environment_variables

    ! Initialise all of the configuration fields that can be derived after
    ! `read_environment_variables` has been called, and allocate the arrays
    ! in `arrays`
    subroutine setup_arrays(config, arrays, nlp)
        use conv_mod, only: nex
        type(t_config), intent(inout) :: config
        type(t_arrays), intent(inout) :: arrays
        integer, intent(in) :: nlp
        integer :: i

        if (allocated(arrays%earx  )) deallocate(arrays%earx  )
        if (allocated(arrays%ReGbar)) deallocate(arrays%ReGbar)
        if (allocated(arrays%ImGbar)) deallocate(arrays%ImGbar)
        allocate(arrays%earx(0:nex))
        allocate(arrays%ReGbar(nex))
        allocate(arrays%ImGbar(nex))

        if (allocated(arrays%contx    )) deallocate(arrays%contx    )
        if (allocated(arrays%contx_int)) deallocate(arrays%contx_int)
        allocate(arrays%contx(nex,nlp))
        allocate(arrays%contx_int(nlp))

        config%dloge = log10(config%Emax / config%Emin) / float(nex)

        ! populate the energy array
        do i = 0,nex
           arrays%earx(i) = config%Emin                                        &
               * (config%Emax/config%Emin)**(float(i)/float(nex))
        end do

    end subroutine setup_arrays

    ! Reallocate arrays depending on whether they need to be resized
    subroutine realloc_arrays(config, model_args, arrays, prev_nf, flosave, fhisave)
        use conv_mod, only: nex
        type(t_config), intent(in) :: config
        type(t_model_arguments), intent(in) :: model_args
        type(t_arrays), intent(inout) :: arrays
        integer, intent(in) :: prev_nf
        integer :: i
        logical :: needs_allocating
        double precision :: fhisave, flosave
        double precision :: fhicheck, flocheck
        double precision, parameter   :: dtol = 1e-5

        ! TODO: magic constants should be named constants
        fhicheck = fhisave /(4.92695275718945d-06 * model_args%Mass)
        flocheck = flosave /(4.92695275718945d-06 * model_args%Mass)

        if ( abs(1. - model_args%floHz/flocheck) .gt. dtol ) then
            needs_allocating = .true.
        else if ( abs(1. - model_args%fhiHz/fhicheck) .gt. dtol ) then
            needs_allocating = .true.
        else if (allocated(arrays%fix)) then
            needs_allocating = prev_nf .ne. config%nf
        else
            needs_allocating = .true.
        endif

        if (needs_allocating) then
            if (allocated(arrays%fix)) deallocate(arrays%fix)
            allocate(arrays%fix(0:config%nf))

            ! populate the frequency array
            do i = 0, config%nf
                arrays%fix(i) = model_args%floHz                               &
                    *(model_args%fhiHz                                         &
                    / model_args%floHz)**(real(i) / real(config%nf))
            end do

            ! reallocate the transfer function arrays
            if (allocated(arrays%ker_W0)) deallocate(arrays%ker_W0)
            allocate(arrays%ker_W0(model_args%nlp,nex,config%nf,config%me,config%xe))

            if (allocated(arrays%ker_W1)) deallocate(arrays%ker_W1)
            allocate(arrays%ker_W1(model_args%nlp,nex,config%nf,config%me,config%xe))

            if (allocated(arrays%ker_W2)) deallocate(arrays%ker_W2)
            allocate(arrays%ker_W2(model_args%nlp,nex,config%nf,config%me,config%xe))

            if (allocated(arrays%ker_W3)) deallocate(arrays%ker_W3)
            allocate(arrays%ker_W3(model_args%nlp,nex,config%nf,config%me,config%xe))

            if (allocated(arrays%ReW0)) deallocate(arrays%ReW0)
            allocate(arrays%ReW0(model_args%nlp,nex,config%nf))

            if (allocated(arrays%ImW0)) deallocate(arrays%ImW0)
            allocate(arrays%ImW0(model_args%nlp,nex,config%nf))

            if (allocated(arrays%ReW1)) deallocate(arrays%ReW1)
            allocate(arrays%ReW1(model_args%nlp,nex,config%nf))

            if (allocated(arrays%ImW1)) deallocate(arrays%ImW1)
            allocate(arrays%ImW1(model_args%nlp,nex,config%nf))

            if (allocated(arrays%ReW2)) deallocate(arrays%ReW2)
            allocate(arrays%ReW2(model_args%nlp,nex,config%nf))

            if (allocated(arrays%ImW2)) deallocate(arrays%ImW2)
            allocate(arrays%ImW2(model_args%nlp,nex,config%nf))

            if (allocated(arrays%ReW3)) deallocate(arrays%ReW3)
            allocate(arrays%ReW3(model_args%nlp,nex,config%nf))

            if (allocated(arrays%ImW3)) deallocate(arrays%ImW3)
            allocate(arrays%ImW3(model_args%nlp,nex,config%nf))

            if (allocated(arrays%ReSraw)) deallocate(arrays%ReSraw)
            allocate(arrays%ReSraw(nex,config%nf))

            if (allocated(arrays%ImSraw)) deallocate(arrays%ImSraw)
            allocate(arrays%ImSraw(nex,config%nf))

            if (allocated(arrays%ReSrawa)) deallocate(arrays%ReSrawa)
            allocate(arrays%ReSrawa(nex,config%nf))

            if (allocated(arrays%ImSrawa)) deallocate(arrays%ImSrawa)
            allocate(arrays%ImSrawa(nex,config%nf))

            if (allocated(arrays%ReGrawa)) deallocate(arrays%ReGrawa)
            allocate(arrays%ReGrawa(nex,config%nf))

            if (allocated(arrays%ImGrawa)) deallocate(arrays%ImGrawa)
            allocate(arrays%ImGrawa(nex,config%nf))

            if (allocated(arrays%ReG)) deallocate(arrays%ReG)
            allocate(arrays%ReG(nex,config%nf))

            if (allocated(arrays%ImG)) deallocate(arrays%ImG)
            allocate(arrays%ImG(nex,config%nf))
        end if

        ! (Re)allocate lensing/reflection fraction arrays if necessary:
        ! TODO: check the length of the arrays, as `nlp` is unlikely to have
        ! changed between calls.
        if (config%needtrans) then
           if (allocated(arrays%lens)) deallocate(arrays%lens)
           allocate (arrays%lens(model_args%nlp))
           if (allocated(arrays%frobs)) deallocate(arrays%frobs)
           allocate (arrays%frobs(model_args%nlp))
           if (allocated(arrays%frrel)) deallocate(arrays%frrel)
           allocate (arrays%frrel(model_args%nlp))
        endif
    end subroutine
end module common_types
