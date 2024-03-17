! Time-stamp: <sensemble.f90 18:43, 22 March 2018 by K Kumari>


!!! $Log: sensemble.f90,v $
!!! Revision 1.3  2004/03/17 02:05:24  pxs565
!!! Working version of a sample demonstration
!!!
!!! Revision 1.2  2004/03/16 21:48:12  pxs565
!!! minor variations, tried using an ntdone file,
!!!
!!! Revision 1.1  2004/03/16 07:01:14  pxs565
!!! Initial revision
!!!


Program chainsim_p
   Use Global_parameters_variables_and_types
   Use Flock_utils
   Use netcdf
   use Initial_Position_Utilities
   use random_numbers
   use Netcdf_Utilities
   use properties
   use Physics_subroutines
   use simulation_utilities
   use Spring_Force_calculatons
   Implicit None
   Include 'mpif.h'
   Integer mpierr, Nprocs, MyId, count

   !_____________________________________________________________
   !        Other Declarations
   !_____________________________________________________________


   Character (len=12) clk(3)
   Integer (k4b) nseed, stored_seed, stored_ix, stored_iy
   !!! netcdf variable
   Character (len = *), Parameter :: FILE_NAME = "trial.nc"
   Integer :: ncid, ncid_VR
   Integer, parameter :: Ndims = 2!, NX =6, NY=12
!   Integer, dimension(:,:), allocatable :: data_out(:,:)
!  Real(DBprec) time_save
   ! File i/o
   Character(10), parameter :: FormatVersion = "GAVG-1.0"
   !Character (10) :: fver
   Character (len=100) infile, outfile, gavgfile, param_name, param_value
   Character (len=200) posfile, ntrajfile, netfile, netfile_VR
   Integer, Parameter ::   inunit=10, outunit=11, gavunit=12, &
      resunit=15, cunit=16, inphi=20, ntrajunit=227, posunit =18, &
      reunit = 19, inangle = 483
   Integer :: read_status

   Logical :: Filexists, Flocked

   ! for gfortran
   Character (len=80) :: Format3, Format31, Format42, &
      Format43, Format47

   ! Trajectories related
   Integer, Parameter :: MaxNDT = 10
   Integer i, clok(8), k
   Integer nthis, nblock, iblock, ntrajdone, ntraj, Nsamples, &
      idelts, ndelts,  ntrajvals(MaxNDT), nthisvals(MaxNDT), ntrajout

   Real (DBprec) t1zimm, tol(MaxNDT), emax
   Real (DBprec) deltseq, deltsne, deltseqvals(10), tlongest, &
      deltsvals(MaxNDT), tmax, t1rouse

   Integer :: bending_potential_type
   Real (DBprec)  :: bending_stiffness
   Real (DBprec)  :: natural_angle_scalar_input
   Real (DBprec), Allocatable  :: natural_angles(:) ! (Nbeads - 2)

   Real (DBprec), Allocatable :: PosVecR(:,:), times(:), samples(:,:), &
      avgs(:,:), errs(:,:), &
      mpp_avgs(:,:), mpp_errs(:,:), &
      global_avgs(:,:), global_errs(:,:), &
      PosVecR_VR(:,:), samples_VR(:,:)
   !!!!! phi parameter declared here
   Real (DBprec), Allocatable :: phi(:,:)
   Real (DBprec) :: phi_data
   !!! Variable for netcdf
   Real (DBprec), Allocatable :: confi(:,:,:), grad(:,:,:), time_cdf(:), &
      confi_VR(:,:,:), grad_VR(:,:,:), time_cdf_VR(:)

   Integer SType,nsact, eqprops, neqprops, ntrajinput

   Real (DBprec) hstar, zstar, dstar, sqrtb, Q0s
   Character (len=10), parameter :: ErString = "Error"
   !Real (DBprec) rems, reerr, rgms, rgerr, xms, xerr, ntot, Conf_t
   Real (DBprec) :: conf_t

   Integer(kind=k4b) :: step, num_steps
   Real(kind=dbprec) :: time, tcur
   Integer :: isample, index_of_min(1)
   Integer, allocatable :: sample_indexes(:)
   Real(kind=dbprec), allocatable :: tcdf(:), sample_times(:)

   Integer :: NBeads, EV, phiFromFile, natural_angle_from_file, COM_update
   Integer FlowTypeProduction
   Integer :: bens, Netcd, Restart, variance_reduction, InitialConfiguration
   Integer :: LookupTableOpt, delSCalcMethod, EigenvalueCalcMethod, ChebUpdateMethod
   Real (DBprec) :: Teq, Tpr, no_chebyshev_terms_multiplier, fd_err_max
   real (DBprec) :: gdots, max_gamma, lookup_table_tol
   Real (DBprec) :: trapOneXPosition, trapTwoXPosition
   Real (DBprec) :: trapTwoXVelocity
   Real (Dbprec) :: trapOneStrength, trapTwoStrength

   ! Simulation structure inputs
   type(physical_parameters) :: phys_params
   type(simulation_parameters) :: sim_params
   type(calculated_variables) :: calc_vars, calc_vars_VR

   Real(Dbprec), dimension(Ndim) :: dummy_cofm


   !-------------------------------------------------------------c
   !        Initialization of MPI                                c
   !-------------------------------------------------------------c
   ! Initializes the parallelization
   Call MPI_INIT(mpierr)
!   write(*,*) "my mpierr", mpierr
   ! Gets the total number of processors
   Call MPI_COMM_SIZE (MPI_COMM_WORLD, Nprocs, mpierr)
!   write(*,*) "Nprocs", Nprocs
   ! Gets the rank or ID of the current processor
   Call MPI_COMM_RANK (MPI_COMM_WORLD, MyId, mpierr)
!  write(*,*) "MyId", MyId

   !_____________________________________________________________
   !        Get input data                                       c
   !_____________________________________________________________

   infile="inputc.dat"

   SType=1
   COM_update=1
   NBeads=5
   hstar=0
   EV=0
   zstar=0
   dstar=1
   bending_potential_type = 0
   bending_stiffness = 0
   natural_angle_from_file = 0
   natural_angle_scalar_input = 0
   sqrtb=7
   Q0s=1
   gdots=0
   emax=0
   bens=0
   Netcd=0
   Restart=0
   Nsamples=10
   phiFromFile=0
   phi_data=0
   Teq=3
   Tpr=5
   ntrajinput=0
   ndelts=1
   variance_reduction=0
   FlowTypeProduction=1
   InitialConfiguration=1
   LookupTableOpt = 0
   lookup_table_tol = 1d-5
   max_gamma = 1000
   ntrajdone = 0
   

   Open (unit=inunit,file=infile,status="old")
   read_status=0
   i=0
   do while (read_status .EQ. 0)
      read(inunit,*,IOSTAT=read_status) param_name, param_value
      IF (read_status > 0) Then
         WRITE(*,*) "Invalid input file. inputc.dat must contain parameter name/value pairs."
         STOP
      END IF
      IF ((read_status .EQ. 0) .AND. (param_name(1:1) .ne. '!')) THEN
         SELECT CASE (param_name)
          CASE ("SpType","SType")
            Read (param_value,*,IOSTAT=read_status) SType
          CASE ("LookupTableOpt")
            Read (param_value,*,IOSTAT=read_status) LookupTableOpt
          CASE ("LookupTableTolerance")
            Read (param_value,*,IOSTAT=read_status) lookup_table_tol
          CASE ("MaxGamma")
            Read (param_value,*,IOSTAT=read_status) max_gamma
          CASE ("NBeads")
            Read (param_value,*,IOSTAT=read_status) NBeads
          CASE ("hstar")
            Read (param_value,*,IOSTAT=read_status) hstar
          CASE ("EV")
            Read (param_value,*,IOSTAT=read_status) EV
          CASE ("zstarByEpsstar","zstar")
            Read (param_value,*,IOSTAT=read_status) zstar
          CASE ("dstar")
            Read (param_value,*,IOSTAT=read_status) dstar
          CASE ("contour_dist_for_EV")
            Read (param_value,*,IOSTAT=read_status) phys_params%EV_inputs%contour_dist_for_EV
          CASE ("min_EV_cutoff")
            Read (param_value,*,IOSTAT=read_status) phys_params%EV_inputs%minimum_EV_cutoff
          CASE ("max_EV_cutoff")
            Read (param_value,*,IOSTAT=read_status) phys_params%EV_inputs%maximum_EV_cutoff
          CASE ("BendingPotential", "bending_potential", "BendingPotentialType")
            Read (param_value,*,IOSTAT=read_status) bending_potential_type
          CASE ("BendingStiffness", "bending_stiffness")
            Read (param_value,*,IOSTAT=read_status) bending_stiffness
          CASE ("natural_angle_data", "natural_angle_scalar_input")
            Read (param_value,*,IOSTAT=read_status) natural_angle_scalar_input
          CASE ("natural_angle_from_file", "NaturalAngleFromFile")
            Read (param_value,*,IOSTAT=read_status) natural_angle_from_file
          CASE ("TrapOneInitialPosition")
            Read (param_value,*,IOSTAT=read_status) trapOneXPosition
          CASE ("TrapTwoInitialPosition")
            Read (param_value,*,IOSTAT=read_status) trapTwoXPosition
          CASE ("TrapOneStrength")
            Read (param_value,*,IOSTAT=read_status) trapOneStrength
          CASE ("TrapTwoStrength")
            Read (param_value,*,IOSTAT=read_status) trapTwoStrength
          CASE ("TrapTwoVelocity")
            Read (param_value,*,IOSTAT=read_status) trapTwoXVelocity
          CASE ("COM_update_on")
            Read (param_value,*,IOSTAT=read_status) COM_update
          CASE ("L0star", "dQ", "sqrtb")
            Read (param_value,*,IOSTAT=read_status) sqrtb
          CASE ("Q0star","Q0s","sigma")
            Read (param_value,*,IOSTAT=read_status) Q0s
          CASE ("FlowCessationTime","flow_cessation_time")
            Read (param_value,*,IOSTAT=read_status) phys_params%Flow_inputs%flow_cessation_time
          CASE ("FlowPeriod","flow_period")
            Read (param_value,*,IOSTAT=read_status) phys_params%Flow_inputs%flow_period
          CASE ("Gdot","gdots")
            Read (param_value,*,IOSTAT=read_status) gdots
          CASE ("Bens","bens")
            Read (param_value,*,IOSTAT=read_status) bens
          CASE ("NetCDF","Netcd")
            Read (param_value,*,IOSTAT=read_status) Netcd
          CASE ("Restart")
            Read (param_value,*,IOSTAT=read_status) Restart
          CASE ("emax")
            Read (param_value,*,IOSTAT=read_status) emax
          CASE ("Nsamp","Nsamples")
            Read (param_value,*,IOSTAT=read_status) Nsamples
          CASE ("phiFile","phiFromFile")
            Read (param_value,*,IOSTAT=read_status) phiFromFile
          CASE ("phiSDK","phi_data")
            Read (param_value,*,IOSTAT=read_status) phi_data
          CASE ("Teq")
            Read (param_value,*,IOSTAT=read_status) Teq
          CASE ("Tpr")
            Read (param_value,*,IOSTAT=read_status) Tpr
          CASE ("Ntrajdone","ntrajinput")
            Read (param_value,*,IOSTAT=read_status) ntrajinput
          CASE ("ndelts")
            Read (param_value,*,IOSTAT=read_status) ndelts
            If (ndelts .GT. MaxNDT) Then
               Write(*,*) 'Number of delt exceeded ', MaxNDT
               Stop
            End If
          CASE ("dtseq")
            i=i+1
            Read (param_value,*,IOSTAT=read_status) deltseqvals(i)
          CASE ("dtsne")
            Read (param_value,*,IOSTAT=read_status) deltsvals(i)
          CASE ("nblock")
            Read (param_value,*,IOSTAT=read_status) nthisvals(i)
          CASE ("ntot")
            Read (param_value,*,IOSTAT=read_status) ntrajvals(i)
          CASE ("tol")
            Read (param_value,*,IOSTAT=read_status) tol(i)
          CASE ("variance_reduction", "VarianceReduction")
            Read (param_value,*,IOSTAT=read_status) variance_reduction
          CASE ("delSCalcMethod", "delS_calc_method")
            Read (param_value,*,IOSTAT=read_status) delSCalcMethod
          CASE ("EigsCalcMethod", "eigs_calc_method")
            Read (param_value,*,IOSTAT=read_status) EigenvalueCalcMethod
          CASE ("nchebMultiplier", "ncheb_multiplier")
            Read (param_value,*,IOSTAT=read_status) no_chebyshev_terms_multiplier
          CASE ("fd_err_max", "FdErrMax")
            Read (param_value,*,IOSTAT=read_status) fd_err_max
          CASE ("cheb_update_method", "ChebUpdateMethod")
            Read (param_value,*,IOSTAT=read_status) ChebUpdateMethod
          CASE ("FlowType")
            Read (param_value,*,IOSTAT=read_status) FlowTypeProduction
          CASE ("InitialConfiguration", "Initial_config")
            Read (param_value,*,IOSTAT=read_status) InitialConfiguration
          CASE DEFAULT
            WRITE(*,*) "Unrecognized parameter: "//trim(param_name)
            STOP
         END SELECT
         IF (read_status .NE. 0) THEN
            WRITE(*,*) "Invalid value for parameter: "//trim(param_name)
            STOP
         END IF
      END IF

   end do

   !_____________________________________________________________
   !        Initialization                                       c
   !_____________________________________________________________
   If (Nsamples.Lt.2) Nsamples = 2

   Allocate(PosVecR(Ndim,NBeads))
   if (variance_reduction .eq. 1) then
      allocate(PosVecR_VR(Ndim,Nbeads))
   end if

   Allocate(phi(NBeads,NBeads))
   if (NBeads.eq.2)  then
      allocate(natural_angles(1))
   else
      allocate(natural_angles(Nbeads-2))
   end if
   !Allocate(data_out(NY,NX))
!  Allocate(confi(Nsamples,Ndim, NBeads), &
!           grad(Nsamples,Ndim, NBeads))

   If (phiFromFile == 1) Then!   phi values read from the file "phi.txt"
      Open (unit=inphi,file="phi.txt",status="old")
      Do k = 1, NBeads
         Read (inphi,*) phi(k,:)
      End Do
      close(unit=inphi)
   Else if (phiFromFile==0) Then!!  phi read from the file "inputc.dat"
      phi = phi_data
   End if
!  write(*,*) 'phi', phi

   If (natural_angle_from_file == 1) Then
      Open (unit=inangle,file="natural_angles.txt",status="old")
      Read (inangle,*) natural_angles(:)
      close(unit=inangle)
   Else if (natural_angle_from_file==0) Then
      natural_angles = natural_angle_scalar_input
   End if

   ! set some simulation parameters
   phys_params%bend_inputs%bending_potential_type = bending_potential_type
   phys_params%bend_inputs%bending_stiffness = bending_stiffness
   phys_params%bend_inputs%natural_angles = natural_angles

   phys_params%EV_inputs%excluded_volume_type = EV
   phys_params%EV_inputs%dimensionless_EV_energy = zstar
   phys_params%EV_inputs%dimensionless_EV_radius = dstar
   phys_params%EV_inputs%bead_attractive_interaction_strengths = phi

   phys_params%Flow_inputs%flow_strength = gdots
   !phys_params%Flow_inputs%flow_type = Flow
   phys_params%Flow_inputs%trapOneInitialPosition = (/trapOneXPosition, 0.d0, 0.d0/)
   phys_params%Flow_inputs%trapTwoInitialPosition = (/trapTwoXPosition, 0.d0, 0.d0/)
   phys_params%Flow_inputs%trapOneStrength = trapOneStrength
   phys_params%Flow_inputs%trapTwoStrength = trapTwoStrength
   phys_params%Flow_inputs%trapTwoVelocity = (/trapTwoXVelocity, 0.d0, 0.d0/)

   phys_params%HI_params%hstar = hstar
   phys_params%HI_params%EigenvalueCalcMethod = EigsFixman
   phys_params%HI_params%delSCalcMethod = Chebyshev
   phys_params%HI_params%nchebMultiplier = no_chebyshev_terms_multiplier
   phys_params%HI_params%ChebUpdateMethod = ChebUpdateMethod
   phys_params%HI_params%fd_err_max = fd_err_max

   phys_params%spring_inputs%finite_extensibility_parameter = sqrtb
   phys_params%spring_inputs%natural_length = Q0s
   phys_params%spring_inputs%spring_type = SType

   phys_params%number_of_beads = NBeads

   sim_params%update_center_of_mass = COM_update

   ! actual production run time

   tmax = Tpr
!  ! longest relaxation times Rouse and Zimm
   t1rouse = 0.5/Sin(PI/2/Nbeads)**2
   t1zimm = lam1_th(hstar,NBeads) ! use thurstons formula for Zimm

   tlongest = t1rouse

!  teqbm1 = Teq*tlongest ! without EV
!  teqbm2 = Teq*tlongest
!!  If (zstar .Ne. 0) Then
!!     teqbm = 3*tlongest ! with excluded volume, it takes roughly 3 times
!     ! to attain equilibrium even with init dist
!!  End If

!   write (*,*) "trouse time is", t1rouse, "tmax is", tmax
   !_____________________________________________________________
   !        Initialization variable format expressions
   !         <> language extension not available in gfortran
   !_____________________________________________________________

   Write(Format3,"(a,I3,a)") "('#',", 2*NProps+4, "(A11,1X))"
   ! short cut way to get (nearly) left justification
   Write(Format31,"(a,I3,a)") "('#',", 2*NProps+4,"(G2.0,10X))"


   If (phys_params%Flow_inputs%flow_strength == 0) Then
      eqprops = 4
   Else
      neqprops = 4
   End if

!!$  Write(Format4 ,"(a,I3,a)") "(", 2*NProps+4,"(G11.4,1X))"
   Write(Format42,"(a,I3,a)") "(", 2*eqprops+1, "(G24.17,1X))"
   Write(Format47,"(a,I3,a)") "(", eqprops+2, "(G24.17,1X))"
   Write(Format43,"(a,I3,a)") "(", 2*neqprops+2, "(G24.17,1X))"

!!$  Write(Format5, "(a,I3,a)") "('#',A11,1X,", nsact, "(G11.4,2X))"
!!$  Write(Format51, "(a,I3,a)") "('#',", nsact+1, "(I3,10X))"
!!$
!!$  Write(Format6, "(a,I3,a)") "(", nsact+1, "(G11.4,2X))"

   !_____________________________________________________________
   !        Initialization of variables
   !_____________________________________________________________

   ! the standard error of mean obtained from t-distribution
   ! for sample mean and sample standard deviation
   Conf_t = 2 ! 95% confidence for degrees of freedom > 20

   Prop_names(1) = "R^2"         ! square of end-to-end vector
   Prop_names(2) = "S11"         ! 1,1 component of shape tensor
   Prop_names(3) = "S12"         ! 1,2 component of shape tensor
   Prop_names(4) = "S13"         ! 1,3 component of shape tensor
   Prop_names(5) = "S22"         ! 2,2 component of shape tensor
   Prop_names(6) = "S23"         ! 2,3 component of shape tensor
   Prop_names(7) = "S33"         ! 3,3 component of shape tensor
   Prop_names(8) = "N1"          ! First normal stress difference
   Prop_names(9) = "N2"          ! Second normal stress difference
   Prop_names(10) = "T12"        ! 1,2 component of stress tensor
   Prop_names(11) = "X1"         ! Stretch in 1 direction
   Prop_names(12) = "X2"         ! stretch in 2 direction
   Prop_names(13) = "X3"         ! stretch in 3 direction
   Prop_names(14) = "Rg^2"       ! sq. Radius of gyration
   Prop_names(15) = "PC error"   ! Predictor-Corrector Error
   Prop_names(16) = "<SxySxy>"   ! Auto correlation fn for shear stress
   Prop_names(17) = "Diffusvty"  ! Center of mass diffusivity
   Prop_names(18) = "PC count"   ! Number of loops taken for predictor-corrector loop
   Prop_names(19) = "ncheb count"   ! Number of chebyshev terms used in approx
   Prop_names(20) = "fd err"   ! Fluctuation-dissipation theorem error for Chebyshev approximation
   Prop_names(21) = "poly err"   ! Error in Chebyshev approximation vs exact representation

!   If (phys_params%Flow_inputs%flow_strength == 0 ) Then
!      Open (unit = resunit, file = "result.dat", status="unknown")
!      Write (resunit,200) "NBeads      ", &
!         "sqrtb       ", &
!         "h*          ", &
!         "z*          ", &
!         "d*          ", &
!         "dtsne       ", &
!         "<R^2>       ", &
!         "<Rg^2>      ", &
!         "Diffsvty    ", &
!         "l_eta       ", &
!         "<x>         "
!      Write (resunit,202) (i, i=1,16)
!200   Format ('#',6(G12.0,2X),5(G12.0,2X, 'Error       ',2X))
!202   Format ('#',16(I2,12X)) !approx LJustification
!   end If

   !_____________________________________________________________
   !     Generate the seeds based on the current time
   !_____________________________________________________________

   Call Date_and_time(clk(1), clk(2), clk(3), clok)
   !       ms              sec         min        hr
   nseed = clok(8)*100000+clok(7)*1000+clok(6)*10+clok(5)

   !-- an unique seed incase clok is the same
   nseed = (nseed + 201271)

   ! Add MPI ID to seed in case two processors reach this point at the exact same time
   nseed = nseed + MyId

   ! Set RNG seed using set_seed
   call set_seed(nseed)

   !_____________________________________________________________
   !    Begin the loop for time step sizes                       c
   !_____________________________________________________________

   timesteps: Do  idelts = 1,ndelts

      deltseq = deltseqvals(idelts)
      deltsne = deltsvals(idelts)
      sim_params%implicit_loop_exit_tolerance = tol(idelts)

      nsact = Nsamples


      ! If lookup table is used, generate the lookup table
      if (LookupTableOpt .eq. 1) then
         if (deltseq .ne. deltsne) then
            print *, "code does not currently support lookup table with different dt for EQ and prod"
            stop
         end if
         call generate_lookup_table(phys_params%spring_inputs, deltseq, max_gamma, lookup_table_tol)
      end if

      ! For NetCDF
      Allocate(confi(nsact,Ndim,NBeads), &
         grad(nsact,Ndim,NBeads), &
         time_cdf(nsact), &
         confi_VR(nsact,Ndim,NBeads), &
         grad_VR(nsact,Ndim,NBeads), &
         time_cdf_VR(nsact))
      confi = 0
      grad = 0
      time_cdf = 0
      confi_VR = 0
      grad_VR = 0
      time_cdf_VR = 0

      allocate(calc_vars%chain_configuration_at_sample_points(nsact,Ndim,Nbeads), &
               calc_vars%total_force_at_sample_points(nsact,Ndim,Nbeads),&
               calc_vars%true_times_at_sample_points(nsact), &
               calc_vars%samples(Nprops,nsact), &
               calc_vars_VR%chain_configuration_at_sample_points(nsact,Ndim,Nbeads), &
               calc_vars_VR%total_force_at_sample_points(nsact,Ndim,Nbeads),&
               calc_vars_VR%true_times_at_sample_points(nsact), &
               calc_vars_VR%samples(Nprops,nsact))

      Allocate(&
         mpp_avgs(NProps,nsact), &
         mpp_errs(NProps,nsact), &
         global_avgs(NProps,nsact), &
         global_errs(NProps,nsact))

      global_avgs = 0
      global_errs = 0
      mpp_avgs = 0
      mpp_errs = 0

      Allocate ( &
         times(nsact),  &
         samples(NProps,nsact), &
         avgs(NProps,nsact), &
         errs(NProps,nsact))

      if (variance_reduction.eq.1) then
         allocate(samples_VR(NProps,nsact))
      end if

      ! nsact = Nsamples
!      Do i = 1,nsact
!         times(i) = (i-1)*mmult*deltsne
!      End Do
      ! fixes an overflow bug for certain values of timestep and samples
      !times(nsact+1) = tmax*10.d0

!      print *, "nsact is: ", nsact
!      print *, "times for samples are : ", times

      avgs = 0.d0
      errs = 0.d0

      ! Obtain the number of traj completed from the disk, if present

!      Write (gavgfile, '("gavgs.",I2.2)') idelts
!
!      Inquire (file=gavgfile, exist=Filexists)
!
!      If (Filexists) Then
!         open (unit=gavunit,file=gavgfile,status='old')
!
!         Read (gavunit,*) fver  ! Format version
!         if (fver /= FormatVersion) Then
!            Write (*,*) 'Incompatible Format version in ', gavgfile
!            close(gavunit)
!            Go to 99999
!         end if
!
!         Read (gavunit,*) ntrajdone
!         close(gavunit)
!
!      Else
!         ntrajdone = 0
!      End If

      !--when several procs start at the same time it
      !  can lead to the total trajs done being > ntrajvals
      ntraj = Max(0,ntrajvals(idelts)-ntrajdone)

      !--number of trajectories for this run
      nthis = Min(ntraj, nthisvals(idelts))

      nblock = nthis

      Write (ntrajfile, '("ntrajdone.",I2.2)') idelts

      ! set up netCDF files
      If (netcd .eq. 1) Then
         write(netfile, '("net_dt",I2.2,"_proc",I3.3,".nc")') idelts, MyId
         write(netfile_VR, '("net_VR_dt",I2.2,"_proc",I3.3,".nc")') idelts, MyId

         call create_netcdf_file(netfile, nsact, Ndim, NBeads, ncid)
         If (variance_reduction .eq. 1) then
            call create_netcdf_file(netfile_VR, nsact, Ndim, NBeads, ncid_VR)
         End If
      End If

      ! set simulation inputs
      !sim_params%implicit_loop_exit_tolerance = Imploop_tol
      sim_params%simulation_seed = nseed
      sim_params%time_at_simulation_start = 0.d0
      !sim_params%times_to_take_samples = times

      !_____________________________________________________________
      !    Begin the loop for the blocks                            c
      !_____________________________________________________________

      trajectories: Do iblock = 1, nblock
         samples = 0.d0
         if (variance_reduction .eq. 1) then
            samples_VR = 0.d0
         end if
         PosVecR = 0.d0
         ntrajout = ntrajinput + (iblock)
         !ntrajout = ntrajout + (MyId*nblock) + (iblock -1)
!       write(*,*) ntrajout, MyId, nblock, iblock
!!!!_______________________________________________________
!           Initial configuration
!!!!_______________________________________________________

         If (Restart .eq. 0) Then
            if (SType.eq.WLC_bounded) then
               Call Initial_position_FENE_Fraenkel(SType,Nbeads,(sqrtb-Q0s),Q0s,PosVecR)
            else if (SType.ne.FENEFraenkel) then
               Call Initial_position(SType,Nbeads,sqrtb,Q0s,PosVecR)
            else
               if (InitialConfiguration.eq.XAxisAligned) then
                  Call FENE_Fraenkel_Aligned_x_axis(SType,Nbeads,sqrtb,Q0s,PosVecR)
               else if (InitialConfiguration.eq.RandomSpherical) then
                  Call Initial_position_FENE_Fraenkel(SType,Nbeads,sqrtb,Q0s,PosVecR)
               end if
            end if
            if (variance_reduction .eq. 1) then
               posVecR_VR = PosVecR
            end if
            ! Restart may be broken - please test before using then delete this comment!
         Else
            call get_config_and_time_from_netcdf(netfile, PosVecR, sim_params%time_at_simulation_start, &
                                                 sim_params%initial_center_of_mass, ntrajout, Restart, Ndim, Nbeads)
            if (variance_reduction .eq. 1) then
               call get_config_and_time_from_netcdf(netfile_VR, PosVecR_VR, sim_params%time_at_simulation_start, &
                                                 dummy_cofm, ntrajout, Restart, Ndim, Nbeads)
            end if
         End If

         phys_params%Flow_inputs%flow_type = FlowTypeProduction
         sim_params%simulation_timestep = deltsne
         sim_params%time_at_simulation_end = tmax
         sim_params%number_of_samples_to_take = Nsamples

         ! calculate correct steps to take samples
         tcur = sim_params%time_at_simulation_start
         num_steps = ceiling((tmax-tcur+sim_params%simulation_timestep/2.d0)/sim_params%simulation_timestep) + 1
         allocate(tcdf(num_steps), sample_times(Nsamples), &
                  sample_indexes(Nsamples), sim_params%sample_indexes(Nsamples))
         time = tcur
         tcdf(1) = tcur
         step = 1
         do while (time .le. (tmax+sim_params%simulation_timestep/2.d0))
            time = time + sim_params%simulation_timestep
            step = step + 1
            tcdf(step) = time
         end do

         do i=0,Nsamples-1
            sample_times(i+1) = tcur + real(i, kind=dbprec)*(tmax-tcur)/(Nsamples-1)
         end do

         sample_indexes(1) = 1
         sample_indexes(Nsamples) = num_steps
         do i=2,Nsamples-1
            index_of_min = minloc(abs(tcdf-sample_times(i)))
            sample_indexes(i) = index_of_min(1)
         end do

         sim_params%sample_indexes = sample_indexes

!         phys_params%Flow_inputs%flow_type = EQ
!         sim_params%simulation_timestep = deltseq
!         sim_params%time_at_simulation_end = teq
!         sim_params%number_of_samples_to_take = nsact
!         if (COM_update .eq. 1) then
!            Call Time_Integrate_Chain(PosVecR, phys_params,sim_params,calc_vars)
!         else
!            Call Time_Integrate_Chain_No_COM_update(PosVecR, phys_params,sim_params,calc_vars)
!         end if

         If (variance_reduction .eq. 1) then
            call get_all_parameters(stored_seed, stored_ix, stored_iy)
            !PosVecR_VR = PosVecR
         End If


         Call Time_Integrate_Chain(PosVecR, phys_params,sim_params,calc_vars)
         confi = calc_vars%chain_configuration_at_sample_points
         samples = calc_vars%samples
         grad = calc_vars%total_force_at_sample_points
         time_cdf = calc_vars%true_times_at_sample_points

         ! This way of implementing VR is fairly dangerous - if another call to the RNG is added
         ! which happens differently depending on the presence of flow or some other effect, then
         ! the stream of random numbers won't be identical, and variance reduction won't work.
         ! Should probably change it to use an identical stream of numbers specifically for the
         ! random vector X_0
         If (variance_reduction .eq. 1) then
            phys_params%Flow_inputs%flow_type = PR
            Call reset_RNG_with_seed(stored_seed, stored_ix, stored_iy)
            Call Time_Integrate_Chain(PosVecR_VR, phys_params,sim_params,calc_vars_VR)
            confi_VR = calc_vars_VR%chain_configuration_at_sample_points
            samples_VR = calc_vars_VR%samples
            grad_VR = calc_vars_VR%total_force_at_sample_points
            time_cdf_VR = calc_vars_VR%true_times_at_sample_points
         End If

         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !!!!!!           Block-ensemble and NetCDF
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!         !_________________Block-ensemble_____________________________
!         If (bens .eq. 1) then
!            write(posfile, '("traj_dt",I2.2,"_proc",I3.3,"_traj",I6.6".nc")') idelts, MyId, ntrajout
!            !Write (posfile, '("traj_", I3.3".txt")') ntrajout + 101
!            Open(unit = posunit, file = posfile, status = 'unknown')
!
!            If (phys_params%Flow_inputs%flow_strength .Eq. 0) Then
!               eqprops = 4
!               Write (posunit,Format47)(times(i), &
!                  samples(1,i), &
!                  samples(14,i), &
!                  samples(17,i), &
!                  samples(11,i), i ,  &
!                  i = 1,nsact )
!
!            Else
!               neqprops = 4
!               Write (posunit,Format3) &
!                  "Time      ",  &
!                  "Strain    ",  &
!                  Prop_names(1), Erstring, &
!                  Prop_names(14), Erstring, &
!                  "psi1     ", Erstring, &
!                  "etap     ", Erstring
!               Write (posunit,Format31) (i, i=1,2+neqprops*2)
!               Write (posunit,Format43) (times(i), times(i)*gdots,  &
!                  global_avgs(1,i), global_errs(1,i), &
!                  global_avgs(14,i), global_errs(14,i), &
!                  -global_avgs(8,i)/gdots/gdots, global_errs(8,i)/gdots/gdots, &
!                  -global_avgs(10,i)/gdots, global_errs(10,i)/gdots, &
!                  i = 1,nsact )
!            End If
!            !_________________Block-ensemble____________________________
!
!            Close (posunit)
!         end if ! for block ensemble

         !___________________NetCDF________________________________
         If (netcd .eq. 1) Then
            call write_to_netcdf(calc_vars, ntrajout, Restart+1, ncid)
            If (variance_reduction .eq. 1) then
               call write_to_netcdf(calc_vars_VR, ntrajout, Restart+1, ncid_VR)
            End If
         End If ! netcdf
         !_________________NetCDF_____________________________________

         !!!!!!!______________________________________________________!!!!!!!!!!!
         !!!!!!!    End Block-ensemble and NetCDF
         !!!!!!!______________________________________________________!!!!!!!!!!


         avgs = avgs + samples

         !
         errs = errs + samples*samples

         deallocate(tcdf, sample_times, sample_indexes, sim_params%sample_indexes)

!        if (variance_reduction .eq. 1) then
!            avgs
!        end if
      End Do trajectories

      ! close netcdf files
      If (netcd .eq. 1) Then
         call close_netcdf_file(ncid)
         If (variance_reduction .eq. 1) then
            call close_netcdf_file(ncid_VR)
         End If
      End If

      !_____________________________________________________________
      !    Consolidate and save final results                       c
      !_____________________________________________________________

      !Create barrier to synchronize all parallel tasks
      Call MPI_BARRIER (MPI_COMM_WORLD, mpierr)

      !-- Wait for a specified time to synchronise, before aborting
      ! call MPE_Timed_Aborting_Barrier(Abortsecs,MPI_Comm_World,mpierr)

      count = NProps*nsact

      Call mpi_reduce (avgs, mpp_avgs, count,      &
         MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, mpierr)

      Call mpi_reduce (errs, mpp_errs, count, &
         MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, mpierr)

      if (MyID == 0) then
         global_avgs = 0
         global_errs = 0
         ntrajdone = 0

         global_avgs = global_avgs + mpp_avgs
         global_errs = global_errs + mpp_errs
         ntrajdone = ntrajdone + (nblock*Nprocs)

         global_avgs = global_avgs/ntrajdone
         If (ntrajdone.Gt.2) Then
            global_errs = Conf_t*(Abs(global_errs-ntrajdone*global_avgs* &
               global_avgs)/ntrajdone/(ntrajdone-1))**0.5
         End If

         ! Write to standard output the average number of semi-implicit loops
         print *, "Average number of loops is: ", sum(global_avgs(18,1:nsact))/nsact

         print *, "Nsamples actually is : ", nsact

         print *, "final time is: ", calc_vars%true_times_at_sample_points(Nsamples)

         print *, "Restart should be: ", Nsamples+Restart

         print *, "total trajectories completed per processor is: ", ntrajout
      end if
!
!
!      If (MyId == 0) Then
!         ! obtain global values from the disk
!         Inquire (file=gavgfile, exist=Filexists)
!         If (Filexists) Then
!            !-- check if it is locked
!            Flocktest: do
!               call islocked(gavgfile,Flocked)
!               if(.not. Flocked) Then
!                  call lockfile(gavunit,gavgfile)
!                  !-- to be unlocked after new data is written
!
!                  open (unit=gavunit,file=gavgfile,status='old')
!                  Read (gavunit,*) fver  ! Format version
!                  if (fver /= FormatVersion) Then
!                     Write (*,*) 'Incompatible Format version in ', gavgfile
!                     call unlockfile(gavunit,gavgfile)
!                     Go to 99999
!                  end if
!
!                  Read (gavunit,*) ntrajdone
!                  Read (gavunit,*) global_avgs(1,:) , global_errs(1,:)
!                  Read (gavunit,*) global_avgs(8,:) , global_errs(8,:)
!                  Read (gavunit,*) global_avgs(10,:), global_errs(10,:)
!                  Read (gavunit,*) global_avgs(11,:), global_errs(11,:)
!                  Read (gavunit,*) global_avgs(14,:), global_errs(14,:)
!
!                  Close (gavunit)
!                  Exit Flocktest
!               Else  ! when File is locked
!                  ! Sleep for a while and try again until it is unlocked,
!                  ! requires -Vaxlib in ifc
!                  call sleep(1)
!               end if
!            end do Flocktest
!         Else ! if Globalavgs does .not. Filexists
!
!            ! lock the file before opening a new one for writing
!            call lockfile(gavunit,gavgfile) ! to be unlocked after new data is written
!            global_avgs = 0
!            global_errs = 0
!            ntrajdone = 0
!
!         end If
!
!         !-- add the current runs to that in the disk
!         global_avgs = global_avgs + mpp_avgs
!         global_errs = global_errs + mpp_errs
!         ntrajdone = ntrajdone + (nblock*Nprocs)
!
!         ! take the sum-total of time-ensemble averages for
!         ! steady equilibrium measurments before averaging
!         ! over ensembles
!
!         ntot = nsact * ntrajdone
!         rems  = Sum(global_avgs(1,1:nsact))/ntot
!         reerr = Conf_t *  &
!            ((Sum(global_errs(1,1:nsact)) - &
!            ntot * rems*rems)/ntot/(ntot-1))**0.5
!
!         rgms  = Sum(global_avgs(14,1:nsact))/ntot
!         rgerr = Conf_t *  &
!            ((Sum(global_errs(14,1:nsact)) - &
!            ntot * rgms*rgms )/ntot/(ntot-1))**0.5
!
!         xms  = Sum(global_avgs(11,1:nsact))/ntot
!         xerr = Conf_t *  &
!            ((Sum(global_errs(11,1:nsact)) - &
!            ntot * xms*xms )/ntot/(ntot-1))**0.5
!
!
!         ntot = (nsact-1)*ntrajdone
!
!         dfvty  = Sum(global_avgs(17,2:nsact))/ntot
!         derr = Conf_t *  &
!            ((Sum(global_errs(17,2:nsact)) - &
!            ntot * dfvty*dfvty)/ntot/(ntot-1))**0.5
!
!         ! average and errors over the ensemble
!         global_avgs = global_avgs/ntrajdone
!         If (ntrajdone.Gt.2) Then
!            global_errs = Conf_t*(Abs(global_errs-ntrajdone*global_avgs* &
!               global_avgs)/ntrajdone/(ntrajdone-1))**0.5
!         End If
!
!         ! Write to standard output the average number of semi-implicit loops
!         print *, "Average number of loops is: ", sum(global_avgs(18,1:nsact))/nsact
!
!         Write (outfile, '("output.",I2.2)') idelts
!         Open (unit = outunit, file = outfile, status="unknown")
!         Write (outunit,1) "SpringTyp ", "NBeads    ", "sqrtb     ", "sigma     ", &
!            "h*        ", "z*        ", "d*        ", &
!            "gdot*     ", "dt*_eq     ", "dt*_ne     ", &
!            "imp-tol   ", "Actual Samples", &
!            "Tmax      ", "NTraj     " , "FlowType  "
!         Write (outunit,2) SType, NBeads, sqrtb, Q0s, hstar, zstar, dstar, &
!            gdots, deltseq, deltsne, sim_params%implicit_loop_exit_tolerance, nsact, &
!            tmax, ntrajdone, phys_params%Flow_inputs%flow_type
!
!         If (gdots .Eq. 0) Then
!            eqprops = 4
!            Write (outunit,Format3) "Time       ",  &
!               Prop_names(1), Erstring, &
!               Prop_names(14), Erstring, &
!               Prop_names(17), Erstring, &
!               Prop_names(11), Erstring
!            Write (outunit,Format31) (i, i=1,1+eqprops*2)
!            Write (outunit,Format42) (times(i), &
!               global_avgs(1,i), global_errs(1,i), &
!               global_avgs(14,i), global_errs(14,i), &
!               global_avgs(17,i), global_errs(17,i), &
!               global_avgs(11,i), global_errs(11,i), &
!               i = 1,nsact )
!
!         Else
!            neqprops = 4
!            Write (outunit,Format3) &
!               "Time      ",  &
!               "Strain    ",  &
!               Prop_names(1), Erstring, &
!               Prop_names(14), Erstring, &
!               "psi1     ", Erstring, &
!               "etap     ", Erstring
!            Write (outunit,Format31) (i, i=1,2+neqprops*2)
!            Write (outunit,Format43) (times(i), times(i)*gdots,  &
!               global_avgs(1,i), global_errs(1,i), &
!               global_avgs(14,i), global_errs(14,i), &
!               -global_avgs(8,i)/gdots/gdots, global_errs(8,i)/gdots/gdots, &
!               -global_avgs(10,i)/gdots, global_errs(10,i)/gdots, &
!               i = 1,nsact )
!         End If
!
!
!         If (gdots .Eq. 0) Then
!            ! The integral of the ensenble averaged Stress-stress correlation
!            ! function gives the intrinsic viscosity at zero shear rate
!            cidx = 16
!
!            Do ord=1,3
!               Call numint(global_avgs(cidx,:), deltsne*mmult, 1, nsact, &
!                  ord, intCss)
!               Call numint(global_errs(cidx,:), deltsne*mmult, 1, nsact, &
!                  ord, csserr)
!               Write (outunit,101)  ord , intCss, csserr + 1./nsact**ord
!101            Format ('# zero sh rate intrinsic visc, ord = ',I2, ':', &
!                  F10.6,' +/- ',F10.6)
!            End Do
!
!            ! lambda_eta from Rouse model
!            l1rouse = 0
!            Do i=1,Nbeads-1
!               l1rouse = l1rouse +  0.5/Sin(i*PI/2/Nbeads)**2
!            End Do
!            Write (outunit,*) '# from Rouse model = ', l1rouse
!            Write (outunit,219) t1rouse,t1zimm
!219         Format('#longest Rouse = ',G10.3, ', Zimm = ',G10.3)
!
!            Write (resunit,210) NBeads, sqrtb, hstar, zstar, dstar, deltsne, &
!               rems, reerr, rgms, rgerr, dfvty, derr, intCss, csserr, xms, xerr
!210         Format (I3,  11X, 16(G12.5,2X))
!
!         end If ! gdots = 0
!
!         Close (unit = outunit)
!
!
!         !_____________________________________________________________
!         ! Save the global averages data on to the disk
!         ! But before that revert to the unnormalised values
!         ! Note that the gavunit file is still locked
!         !_____________________________________________________________
!         Open(unit=gavunit,file=gavgfile,status='replace')
!
!         If (ntrajdone.Gt.2) Then
!            global_errs = (global_errs/Conf_t)**2 * ntrajdone * (ntrajdone-1)&
!               + ntrajdone * global_avgs**2
!         end If
!
!         global_avgs = global_avgs * ntrajdone
!
!         Write (gavunit,*) FormatVersion
!         Write (gavunit,*) ntrajdone
!         Write (gavunit,*) global_avgs(1,:) , global_errs(1,:)
!         Write (gavunit,*) global_avgs(8,:) , global_errs(8,:)
!         Write (gavunit,*) global_avgs(10,:), global_errs(10,:)
!         Write (gavunit,*) global_avgs(11,:), global_errs(11,:)
!         Write (gavunit,*) global_avgs(14,:), global_errs(14,:)
!         Close (gavunit)
!
!         !-- global avgs file is released for use with other processors
!         call unlockfile(gavunit,gavgfile)
!
!         !_____________________________________________________________
!         !           Write info to disk for continuation
!         !_____________________________________________________________
!         ntraj = ntrajvals(idelts)-ntrajdone
!
!         Write (contfile, '("continue.",I2.2)') idelts
!         open(unit=cunit,file=contfile,status='replace')
!
!         if (ntraj > 0) Then
!            Write (cunit,*) ntraj, '  more trajectories to be completed'
!            close(cunit)
!         Else
!            close(cunit,status='delete')
!         end if
!
!      End If

      deallocate(&
         mpp_avgs, &
         mpp_errs, &
         global_avgs, &
         global_errs)

      deallocate ( &
         PosVecR, &
         times,  &
         samples, &
         avgs, &
         errs)

      if (variance_reduction.eq.1) then
         deallocate(PosVecR_VR, samples_VR)
      end if

      deallocate(phi)

      deallocate(confi, grad, time_cdf)
      deallocate(confi_VR, grad_VR, time_cdf_VR)

      deallocate(phys_params%EV_inputs%bead_attractive_interaction_strengths)

      deallocate(calc_vars%chain_configuration_at_sample_points, &
               calc_vars%total_force_at_sample_points,&
               calc_vars%true_times_at_sample_points, &
               calc_vars%samples, &
               calc_vars_VR%chain_configuration_at_sample_points, &
               calc_vars_VR%total_force_at_sample_points,&
               calc_vars_VR%true_times_at_sample_points, &
               calc_vars_VR%samples)

   End Do timesteps

!   Close(resunit)

   deallocate(natural_angles)

   ! Common termination statements
!99999   Stop
   Call MPI_FINALIZE(mpierr)
   !_____________________________________________________________
   !                    Format statements                        c
   !_____________________________________________________________

!1  Format ('#',16(A10,2X))
!   !2  Format ('#',2(I10,2X), 12(G10.3,2X), (I10,2X)  )
!2  Format ('#',2(G8.0,2X), 12(G10.2,2X), (G10.2,2X), G8.0,2X)

!3 Format ('#', <2*NProps+4>(A11,1X))
!3  Format (Format3)
   ! short cut way to get (nearly) left justification
!31 Format ('#', <2*NProps+4>(G2.0,10X))
!31 Format (Format31)

!4 Format (<2*NProps+4>(G11.4,1X))
!42 Format (<2*eqprops+1>(G11.4,1X))
!43 Format (<2*neqprops+2>(G11.4,1X))
!4  Format (Format4)
!42 Format (Format42)
!43 Format (Format43)

!5 Format ('#',A11,1X,<nsact>(G11.4,2X))
!51 Format ('#', <nsact+1>(I3,10X))
!5  Format (Format5)
!51 Format (Format51)

!6 Format (<nsact+1>(G11.4,2X))
!6  Format (Format6)

End Program chainsim_p


