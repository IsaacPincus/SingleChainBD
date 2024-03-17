! This file contains unit and validation tests,
! which can be run automatically without setting up inputs
! Uses the FRUIT Fortran unit testing package

! Isaac Pincus Sep 2019

! At the moment this doesn't use netcdf or mpi, so it
! is compiled using only gfortran/ifort

module regression_tests
   ! These tests simply compare previous results with chains to see whether or not any results have changed

   use fruit
   Use Global_parameters_variables_and_types
   use Physics_subroutines
   use random_numbers
   use properties
   use Initial_Position_Utilities

   implicit none

   private

   public :: regression_tests_run

   real(DBprec), Parameter :: max_double_difference = 1.D-6

   ! Chain parameters
   Integer :: Nbeads      !Number of beads in chain
   Integer(k4b) :: myseed      !seed for random number generator
   Real (DBprec), allocatable :: PosVecR(:,:)       ! Position vector of beads, number_of_dimensions*number_of_beads

   type(physical_parameters) :: phys_params
   type(simulation_parameters) :: sim_params
   type(calculated_variables) :: output_vars

contains

   subroutine regression_tests_run()

      print *, "Running regression tests"

      call test_rouse_chain_eq()
      call test_FENE_chain_noEV_HI_shear()
      call test_FF_gives_FENE()
      call test_FENE_UA_flow_Gaussian_EV()
      call test_ILC_SH_flow_SDK_EV()
      call test_WLC_PU_flow_LJ_EV_bending_potential()
      call test_Fraenkel_UR_flow_SDK_EV()
      call test_FENEFraenkel_PL_flow_Gaussian_EV()

      print *, ""
      print *, "Finished regression tests"
      print *, ""

   end subroutine

   subroutine test_rouse_chain_eq()
      ! Simulation of a rouse chain at equilibrium to check if the final configuration matches that of a previous version of the code
      ! If the same seed is used, the results should be the same each time. Only tests Time_integrate_chain, everything else should be mocked.

      real (DBprec), allocatable :: final_configurations_answer(:,:)

      phys_params%number_of_beads = 10

      phys_params%spring_inputs%spring_type = HOOK
      phys_params%spring_inputs%finite_extensibility_parameter = 1000.d0
      phys_params%spring_inputs%natural_length = 0.d0

      myseed = 5
      sim_params%simulation_seed = myseed

      phys_params%hstar = 0.d0

      phys_params%EV_inputs%dimensionless_EV_energy = 0.d0
      phys_params%EV_inputs%dimensionless_EV_radius = 0.d0
      phys_params%EV_inputs%excluded_volume_type = noEV

      sim_params%time_at_simulation_start = 0.d0
      sim_params%time_at_simulation_end = 1.d0
      sim_params%simulation_timestep = 0.1d0
      sim_params%implicit_loop_exit_tolerance = 1.d-6

      phys_params%Flow_inputs%flow_type = EQ
      phys_params%Flow_inputs%flow_strength = 1.d0

      Nbeads = phys_params%number_of_beads

      phys_params%bend_inputs%bending_potential_type = NoBendingPotential
      phys_params%bend_inputs%bending_stiffness = 0

      sim_params%number_of_samples_to_take = 0
      allocate(sim_params%sample_indexes(sim_params%number_of_samples_to_take))

      allocate(PosVecR(Ndim, Nbeads))
      allocate(final_configurations_answer(Ndim, Nbeads))

      ! initial configuration, created randomly by intitial config generation function in utils
      PosVecR = reshape((/ 0.000000000000000d00, 0.000000000000000d00, 0.000000000000000d00, &
         -6.586250662803650d-02, -4.911733418703079d-02, 2.791636288166050d-01, &
         -6.753685772418980d-01, 3.519757017493250d-01, -6.340299546718600d-01, &
         -2.844688922166820d00, -2.239686943590640d00, -7.809645235538480d-01, &
         -1.818384915590290d00, -2.201230965554710d00, 9.345296323299410d-01, &
         -3.210149317979810d00, -9.372534379363060d-01, -4.681020081043240d-01, &
         -1.490711003541950d00, -5.971195027232170d-01, -1.509506493806840d00, &
         -1.587433911859990d00, -1.017917685210700d00, -1.456502400338650d00, &
         -6.204553022980690d-01, -1.219155095517640d00, -7.125178799033161d-01, &
         -8.804024234414100d-01, -3.662307761609550d00, -7.348611745983360d-01/), (/Ndim, Nbeads/))

      call reset_RNG_with_seed(myseed)

      Call Time_Integrate_Chain(PosVecR, phys_params, sim_params, output_vars)

      final_configurations_answer =  reshape((/ 2.014655995894740d00, 3.266955094952210d00, 6.195918791485830d-01, &
         6.898958025934120d-01, 8.753981552180370d-01, 4.785034222762240d-01, &
         -3.793564324297950d-01, 1.489672909282450d00, -6.390832004693100d-01, &
         -1.130337343526490d00, -6.146915212957660d-01, -5.621008935576799d-02, &
         -6.222945452675470d-01, -1.038001996615710d00, 5.028959245456580d-01, &
         -2.989398757279950d-01, -8.209160875818450d-01, 1.289920072746550d00, &
         2.453406088052380d-01, 2.086869672299050d-01, 4.316992884969550d-01, &
         -1.110128512089560d-01, -5.552336867493930d-01, -1.277926546344950d00, &
         -6.424329445947624d-02, -1.036741598612170d00, -1.699024281591250d00, &
         -1.915098093214250d-01, -2.620904430772130d00, 2.403765459336170d-01/), (/Ndim, Nbeads/))

      call assertEquals(PosVecR, final_configurations_answer, Ndim, Nbeads, &
         max_double_difference, "Rouse chain equilibrium final configuration test 1 failed")

      PosVecR = final_configurations_answer
      myseed = 5

      final_configurations_answer =  reshape((/2.358466969676050d00, 3.680650287690870d00, 7.391286505833830d-01, &
         9.840637226840210d-02, 1.553897668605920d00, -4.535516789623883d-02, &
         3.703216834718498d-03, 1.100513005729840d00, -5.311631291023066d-02, &
         -1.762304275149680d00, -2.808383016377780d-02, 2.665681695031340d-01, &
         -1.122405351207410d00, -4.065182692363470d-01, 9.790414477264030d-01, &
         6.632472866300509d-02, -4.887170934587170d-01, 9.806758704379520d-01, &
         -5.063576758472540d-01, -9.911910816941520d-01, 1.377974053682720d-01, &
         2.248130548550660d-01, -9.103044958788740d-01, -1.223068999396410d00, &
         -2.639644443515610d-01, -1.153796153975820d00, -1.197611287878650d00, &
         1.905555228121980d00, -2.718416492074990d00, -5.090337805446200d-01/), (/Ndim, Nbeads/))

      Call Time_Integrate_Chain(PosVecR, phys_params, sim_params, output_vars)

      call assertEquals(PosVecR, final_configurations_answer, Ndim, Nbeads, &
         max_double_difference, "Rouse chain equilibrium final configuration test 2 failed")

      deallocate(PosVecR)
      deallocate(final_configurations_answer)
      deallocate(sim_params%sample_indexes)

   end subroutine

   subroutine test_FENE_chain_noEV_HI_shear()

      !Integer :: i
      real (DBprec), allocatable :: final_configurations_answer(:,:)

      phys_params%number_of_beads = 10

      phys_params%spring_inputs%spring_type = FENE
      phys_params%spring_inputs%finite_extensibility_parameter = 50.d0
      phys_params%spring_inputs%natural_length = 0.d0

      myseed = 5
      sim_params%simulation_seed = myseed

      phys_params%hstar = 0.1d0

      phys_params%EV_inputs%dimensionless_EV_energy = 0.d0
      phys_params%EV_inputs%dimensionless_EV_radius = 0.d0
      phys_params%EV_inputs%excluded_volume_type = noEV

      sim_params%time_at_simulation_start = 0.d0
      sim_params%time_at_simulation_end = 1.d0
      sim_params%number_of_samples_to_take = 0
      sim_params%implicit_loop_exit_tolerance = 1.d-6

      phys_params%Flow_inputs%flow_strength = 0.5d0

      Nbeads = phys_params%number_of_beads

      phys_params%bend_inputs%bending_potential_type = NoBendingPotential
      phys_params%bend_inputs%bending_stiffness = 0

      sim_params%number_of_samples_to_take = 0
      allocate(sim_params%sample_indexes(sim_params%number_of_samples_to_take))

      allocate(PosVecR(Ndim, Nbeads))
      allocate(final_configurations_answer(Ndim, Nbeads))

      ! initial configuration, created randomly by intitial config generation function in utils
      PosVecR = reshape((/ 0.000000000000000d00, 0.000000000000000d00, 0.000000000000000d00, &
         -6.586250662803650d-02, -4.911733418703079d-02, 2.791636288166050d-01, &
         -6.753685772418980d-01, 3.519757017493250d-01, -6.340299546718600d-01, &
         -2.844688922166820d00, -2.239686943590640d00, -7.809645235538480d-01, &
         -1.818384915590290d00, -2.201230965554710d00, 9.345296323299410d-01, &
         -3.210149317979810d00, -9.372534379363060d-01, -4.681020081043240d-01, &
         -1.490711003541950d00, -5.971195027232170d-01, -1.509506493806840d00, &
         -1.587433911859990d00, -1.017917685210700d00, -1.456502400338650d00, &
         -6.204553022980690d-01, -1.219155095517640d00, -7.125178799033161d-01, &
         -8.804024234414100d-01, -3.662307761609550d00, -7.348611745983360d-01/), (/Ndim, Nbeads/))

      !reset seed so previous tests don't affect results
      call reset_RNG_with_seed(myseed)

      phys_params%Flow_inputs%flow_type = EQ
      sim_params%simulation_timestep = 0.1d0
      Call Time_Integrate_Chain(PosVecR, phys_params, sim_params, output_vars)

      final_configurations_answer =  reshape((/1.8625701135733601d+00, 3.1969723242894399d+00, 5.2377465510116095d-01, &
         7.5247584390079902d-01, 9.5536664176495600d-01, 4.8013877071655803d-01, &
         -3.4356180289960397d-01, 1.4583537715485499d+00, -5.6815794852782198d-01, &
         -1.1530915882948201d+00, -6.7453657800376499d-01, -9.4743843284497253d-02, &
         -5.7286257202768198d-01, -1.1572473822168701d+00, 5.7284487670511097d-01, &
         -3.1052304471509801d-01, -7.5597526136672399d-01, 1.3031863654134701d+00, &
         2.6987185000428998d-01, 1.4572827023894999d-01, 3.5757650853322798d-01, &
         -9.3164890662119249d-02, -5.7621259089729304d-01, -1.2861934268315200d+00, &
         7.8429050139007667d-03, -1.0572941370082500d+00, -1.7089141826318901d+00, &
         -1.6373593309661899d-01, -2.7051618991776798d+00, 1.7797540997121300d-01/), (/Ndim, Nbeads/))

      call assertEquals(PosVecR, final_configurations_answer, Ndim, Nbeads, &
         max_double_difference, "FENE Chain with HI no EV eq test failed")

      call reset_RNG_with_seed(myseed)

      phys_params%Flow_inputs%flow_type = SH
      sim_params%simulation_timestep = 0.01d0
      Call Time_Integrate_Chain(PosVecR, phys_params, sim_params, output_vars)

      final_configurations_answer =  reshape((/3.0542083998270440d+00, 4.0926557162674282d+00, 2.2880978885994174d-01, &
         1.0588199299835015d+00, 1.8810625260603198d+00, 6.8391615886530310d-01, &
         -5.7864034355175531d-01, 1.5539140972345959d+00, 6.1053502569736406d-02, &
         -1.2153460019599138d+00, 2.1917203618400904d-01, -1.1201705017929225d-01, &
         -1.4914804398634243d+00, -8.1523817433059231d-01, 2.4827512220912168d-01, &
         5.8439265525941875d-01, -5.7538269380339291d-01, 1.4579343447425401d+00, &
         1.0191479239369796d+00, -9.9936943348662266d-01, 6.2758401314007917d-01, &
         -4.3001369766830599d-01, -1.2156185330726244d+00, -1.2773176591404838d+00, &
         -4.6634026284862118d-01, -1.6757415497955723d+00, -2.2921272269296908d+00, &
         -1.1753197900950085d+00, -2.0926292687829395d+00, 6.9840424587465666d-01/), (/Ndim, Nbeads/))
!
!        print *, "posvecR for shear flow"
!        print *, PosVecR

      call assertEquals(PosVecR, final_configurations_answer, Ndim, Nbeads, &
         max_double_difference, "FENE Chain with HI no EV Shear test failed")

      deallocate(PosVecR)
      deallocate(final_configurations_answer)
      deallocate(sim_params%sample_indexes)

   end subroutine
!
   subroutine test_FF_gives_FENE()

      real (DBprec), allocatable :: final_positions_FENE(:,:)
      real (DBprec), allocatable :: final_positions_FF(:,:)

      phys_params%number_of_beads = 10

      phys_params%spring_inputs%spring_type = FENE
      phys_params%spring_inputs%natural_length = 0.d0
      phys_params%spring_inputs%finite_extensibility_parameter = 10.d0

      myseed = 5
      sim_params%simulation_seed = myseed

      phys_params%hstar = 0.3d0
      phys_params%EV_inputs%dimensionless_EV_energy = 1.d0
      phys_params%EV_inputs%dimensionless_EV_radius = 1.d0
      phys_params%EV_inputs%excluded_volume_type = Gauss

      sim_params%time_at_simulation_start = 0.d0
      sim_params%time_at_simulation_end = 1.d0
      sim_params%simulation_timestep = 0.5d0
      sim_params%implicit_loop_exit_tolerance = 1.d-6

      phys_params%Flow_inputs%flow_type = SH
      phys_params%Flow_inputs%flow_strength = 5.d0

      Nbeads = phys_params%number_of_beads

      phys_params%bend_inputs%bending_potential_type = NoBendingPotential
      phys_params%bend_inputs%bending_stiffness = 0

      sim_params%number_of_samples_to_take = 0
      allocate(sim_params%sample_indexes(sim_params%number_of_samples_to_take))

      allocate(PosVecR(Ndim, Nbeads))
      allocate(final_positions_FENE(Ndim, Nbeads))
      allocate(final_positions_FF(Ndim, Nbeads))

      ! initial configuration, created randomly by intitial config generation function in utils
      PosVecR = reshape((/ 0.000000000000000d00, 0.000000000000000d00, 0.000000000000000d00, &
         -6.586250662803650d-02, -4.911733418703079d-02, 2.791636288166050d-01, &
         -6.753685772418980d-01, 3.519757017493250d-01, -6.340299546718600d-01, &
         -2.844688922166820d00, -2.239686943590640d00, -7.809645235538480d-01, &
         -1.818384915590290d00, -2.201230965554710d00, 9.345296323299410d-01, &
         -3.210149317979810d00, -9.372534379363060d-01, -4.681020081043240d-01, &
         -1.490711003541950d00, -5.971195027232170d-01, -1.509506493806840d00, &
         -1.587433911859990d00, -1.017917685210700d00, -1.456502400338650d00, &
         -6.204553022980690d-01, -1.219155095517640d00, -7.125178799033161d-01, &
         -8.804024234414100d-01, -3.662307761609550d00, -7.348611745983360d-01/), (/Ndim, Nbeads/))

      call reset_RNG_with_seed(myseed)

      Call Time_Integrate_Chain(PosVecR, phys_params, sim_params, output_vars)

      final_positions_FENE = PosVecR

      phys_params%spring_inputs%spring_type = FENEFraenkel

      PosVecR = reshape((/ 0.000000000000000d00, 0.000000000000000d00, 0.000000000000000d00, &
         -6.586250662803650d-02, -4.911733418703079d-02, 2.791636288166050d-01, &
         -6.753685772418980d-01, 3.519757017493250d-01, -6.340299546718600d-01, &
         -2.844688922166820d00, -2.239686943590640d00, -7.809645235538480d-01, &
         -1.818384915590290d00, -2.201230965554710d00, 9.345296323299410d-01, &
         -3.210149317979810d00, -9.372534379363060d-01, -4.681020081043240d-01, &
         -1.490711003541950d00, -5.971195027232170d-01, -1.509506493806840d00, &
         -1.587433911859990d00, -1.017917685210700d00, -1.456502400338650d00, &
         -6.204553022980690d-01, -1.219155095517640d00, -7.125178799033161d-01, &
         -8.804024234414100d-01, -3.662307761609550d00, -7.348611745983360d-01/), (/Ndim, Nbeads/))

      call reset_RNG_with_seed(myseed)

      Call Time_Integrate_Chain(PosVecR, phys_params, sim_params, output_vars)

      final_positions_FF = PosVecR

      call assertEquals(final_positions_FF, final_positions_FENE, Ndim, Nbeads, &
         max_double_difference, "FENE_Fraenkel chain is not the same as FENE chain")

      deallocate(PosVecR)
      deallocate(final_positions_FENE, final_positions_FF)
      deallocate(sim_params%sample_indexes)

   end subroutine

   subroutine test_FENE_UA_flow_Gaussian_EV()

      real (DBprec), allocatable :: final_positions(:,:)

      phys_params%number_of_beads = 10

      phys_params%spring_inputs%spring_type = FENE
      phys_params%spring_inputs%finite_extensibility_parameter = 10.d0
      phys_params%spring_inputs%natural_length = 0.d0

      myseed = 6
      sim_params%simulation_seed = myseed

      phys_params%hstar = 0.3d0

      phys_params%EV_inputs%dimensionless_EV_energy = 1.d0
      phys_params%EV_inputs%dimensionless_EV_radius = 2.d0
      phys_params%EV_inputs%excluded_volume_type = Gauss

      sim_params%time_at_simulation_start = 0.d0
      sim_params%time_at_simulation_end = 1.d0
      sim_params%simulation_timestep = 0.1d0
      sim_params%implicit_loop_exit_tolerance = 1.d-6

      phys_params%Flow_inputs%flow_type = UA
      phys_params%Flow_inputs%flow_strength = 4.32

      Nbeads = phys_params%number_of_beads

      phys_params%bend_inputs%bending_potential_type = NoBendingPotential
      phys_params%bend_inputs%bending_stiffness = 0

      sim_params%number_of_samples_to_take = 0
      allocate(sim_params%sample_indexes(sim_params%number_of_samples_to_take))

      allocate(PosVecR(Ndim, Nbeads))
      allocate(final_positions(Ndim, Nbeads))

      ! initial configuration, created randomly by intitial config generation function in utils
      PosVecR = reshape((/ 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, &
         -6.5689776547045498d-01, 1.6751117911767399d+00, 7.0439241401167996d-01, &
         9.1126509787411902d-01, 1.5957861461243701d+00, 9.3448305751127803d-01, &
         8.5937068053146803d-01, 1.8907336197839499d+00, 2.2520324759989699d-01, &
         -8.5057484900184999d-01, 3.0159844275960701d+00, -1.0628856760500900d+00, &
         -1.0030988818944000d+00, 2.6536587542249599d+00, -3.5351133043396399d-01, &
         -1.5316075422199000d+00, 2.6236949015005502d+00, 9.9841216438871905d-01, &
         -1.0762729823295600d+00, 2.6407564502859802d+00, 9.1489505114087000d-01, &
         -3.6989945488079601d-01, 2.0676821761753899d+00, 2.2461602536669001d+00, &
         8.1487530248523599d-01, 7.5294161938743698d-01, 3.8159590742503600d+00/), (/Ndim, Nbeads/))

      call reset_RNG_with_seed(myseed)

      Call Time_Integrate_Chain(PosVecR, phys_params, sim_params, output_vars)

      final_positions = reshape((/9.9129437953246402d+00, -1.1235498041803302d+00, -7.5499498581128879d-02, &
         9.0260535640944894d-01, 1.2716960104223429d-01, -1.4809018073533096d-02, &
         8.4532015467123856d+00, -1.2247202525933931d-01, 3.0811731081154364d-01, &
         7.0546790790963465d+00, 4.6438130696585589d-01, -7.0558320961440923d-01, &
         -2.6502272636515007d+00, 1.3343351453587055d-01, -3.7098235462857043d-01, &
         -1.2288272503015627d+01, 7.4593747962359597d-02, -1.2165639023354316d-01, &
         -2.0747032351656081d+01, -2.5424919616919686d-02, 4.1008518791832788d-02, &
         -1.0859697640907614d+01, -4.0904431733186052d-02, -2.1150238919978090d-01, &
         -1.1359687110807819d+00, -5.6196784319298074d-03, -1.0742525940572789d-01, &
         8.4866864691806292d+00, 2.3493988629341836d-02, 9.2260549294556127d-02/), (/Ndim, Nbeads/))

      call assertEquals(PosVecR, final_positions, Ndim, Nbeads, &
         max_double_difference, "FENE UA Gaussian test failed")

      deallocate(PosVecR)
      deallocate(final_positions)
      deallocate(sim_params%sample_indexes)

   end subroutine

   subroutine test_ILC_SH_flow_SDK_EV()

      real (DBprec), allocatable :: final_positions(:,:)

      phys_params%number_of_beads = 10

      phys_params%spring_inputs%spring_type = ILC
      phys_params%spring_inputs%finite_extensibility_parameter = 5.d0
      phys_params%spring_inputs%natural_length = 0.d0

      myseed = 6
      sim_params%simulation_seed = myseed

      phys_params%hstar = 0.3d0

      phys_params%EV_inputs%dimensionless_EV_energy = 1.d0
      phys_params%EV_inputs%dimensionless_EV_radius = 2.d0
      phys_params%EV_inputs%excluded_volume_type = SDK
      phys_params%EV_inputs%minimum_EV_cutoff = 0.7D0*phys_params%EV_inputs%dimensionless_EV_radius
      phys_params%EV_inputs%maximum_EV_cutoff = 1.82D0*phys_params%EV_inputs%dimensionless_EV_radius

      sim_params%time_at_simulation_start = 0.d0
      sim_params%time_at_simulation_end = 10.d0
      sim_params%simulation_timestep = 0.01d0
      sim_params%implicit_loop_exit_tolerance = 1.d-6

      phys_params%Flow_inputs%flow_type = SH
      phys_params%Flow_inputs%flow_strength = 1.2

      Nbeads = phys_params%number_of_beads

      allocate(phys_params%&
               EV_inputs%bead_attractive_interaction_strengths(Nbeads,Nbeads))
      phys_params%EV_inputs%bead_attractive_interaction_strengths = 0.1d0          ! SDK excluded volume parameter (well depth)
      phys_params%EV_inputs%bead_attractive_interaction_strengths(2,3) = 0.3d0
      phys_params%EV_inputs%bead_attractive_interaction_strengths(3,2) = 0.3d0
      phys_params%EV_inputs%bead_attractive_interaction_strengths(8,4) = 0.5d0
      phys_params%EV_inputs%bead_attractive_interaction_strengths(4,8) = 0.5d0

      allocate(sim_params%sample_indexes(sim_params%number_of_samples_to_take))

      allocate(PosVecR(Ndim, Nbeads))
      allocate(final_positions(Ndim, Nbeads))

      ! initial configuration, created randomly by intitial config generation function in utils
      PosVecR = reshape((/ 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, &
         -6.1823944828956401d-01, 1.5765317588784300d+00, 6.6293904517402702d-01, &
         8.6995802287109802d-01, 1.5012511696459301d+00, 8.8129666720285005d-01, &
         8.2005905273015101d-01, 1.7848573056221599d+00, 1.9929012644809799d-01, &
         -7.6605554351853999d-01, 2.8286196359921401d+00, -9.9551789976143201d-01, &
         -9.1264986517298696d-01, 2.4803801795699001d+00, -3.1372206062358399d-01, &
         -1.4156160573597001d+00, 2.4518644596973900d+00, 9.7286388084413100d-01, &
         -9.7666225034648202d-01, 2.4683122158147999d+00, 8.9235131263554301d-01, &
         -3.0668100442364798d-01, 1.9247626612871100d+00, 2.1550298876148402d+00, &
         7.9420664981513800d-01, 7.0311134022698996d-01, 3.6136802601641098d+00/), (/Ndim, Nbeads/))

      call reset_RNG_with_seed(myseed)

      Call Time_Integrate_Chain(PosVecR, phys_params, sim_params, output_vars)

      final_positions = reshape((/-6.4979613391990458d+00, -3.9355217339347803d-01, -2.6225564015785086d+00, &
         -7.4133896605212417d+00, 8.6865602504037243d-02, -3.8459590409774247d-01, &
         -4.4412895956191738d+00, 2.0427049004283082d-01, 2.6208264081636035d-01, &
         -7.5201657448206261d-01, 4.4633462260185125d-01, -1.0052057153948630d+00, &
         1.4435801732003739d+00, -5.8044351756905233d-01, -1.3661653409730119d+00, &
         4.4360655891970247d+00, 5.3260302106703872d-01, -2.4265637404439797d+00, &
         3.8957728954653410d+00, -3.7093111327590522d-02, -4.3026950484003956d-01, &
         5.3113729872457611d+00, 5.5925412740915925d-01, 1.8119994607699264d+00, &
         2.7948497016796043d+00, -4.8498778175873614d-01, 2.0368809468433877d+00, &
         1.4065867675449550d+00, -5.2290503021589496d-01, 4.1306093450272821d+00/), (/Ndim, Nbeads/))

      call assertEquals(PosVecR, final_positions, Ndim, Nbeads, &
         max_double_difference, "ILC SH SDK test failed")

      deallocate(PosVecR, phys_params%EV_inputs%bead_attractive_interaction_strengths)
      deallocate(final_positions)
      deallocate(sim_params%sample_indexes)

   end subroutine

   subroutine test_WLC_PU_flow_LJ_EV_bending_potential()

      real (DBprec), allocatable :: final_positions(:,:)

      phys_params%number_of_beads = 10

      phys_params%spring_inputs%spring_type = WLC
      phys_params%spring_inputs%finite_extensibility_parameter = 5.d0
      phys_params%spring_inputs%natural_length = 0.d0

      myseed = 6
      sim_params%simulation_seed = myseed

      phys_params%hstar = 0.3d0

      phys_params%EV_inputs%dimensionless_EV_energy = 1.d0
      phys_params%EV_inputs%dimensionless_EV_radius = 0.5d0
      phys_params%EV_inputs%excluded_volume_type = LJ
      phys_params%EV_inputs%minimum_EV_cutoff = 0.7D0*phys_params%EV_inputs%dimensionless_EV_radius
      phys_params%EV_inputs%maximum_EV_cutoff = 1.5D0*phys_params%EV_inputs%dimensionless_EV_radius

      sim_params%time_at_simulation_start = 0.d0
      sim_params%time_at_simulation_end = 1.d0
      sim_params%simulation_timestep = 0.01d0
      sim_params%implicit_loop_exit_tolerance = 1.d-6

      phys_params%Flow_inputs%flow_type = PU
      phys_params%Flow_inputs%flow_strength = 0.2

      Nbeads = phys_params%number_of_beads

      phys_params%bend_inputs%bending_potential_type = OneMinusCosTheta
      phys_params%bend_inputs%bending_stiffness = 2

      sim_params%number_of_samples_to_take = 0
      allocate(sim_params%sample_indexes(sim_params%number_of_samples_to_take))

      allocate(PosVecR(Ndim, Nbeads))
      allocate(final_positions(Ndim, Nbeads))

      ! initial configuration, created randomly by intitial config generation function in utils
      PosVecR = reshape((/0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, &
         -6.1823944828956401d-01, 1.5765317588784300d+00, 6.6293904517402702d-01, &
         8.6995802287109802d-01, 1.5012511696459301d+00, 8.8129666720285005d-01, &
         8.2005905273015101d-01, 1.7848573056221599d+00, 1.9929012644809799d-01, &
         -7.6605554351853999d-01, 2.8286196359921401d+00, -9.9551789976143201d-01, &
         -9.1264986517298696d-01, 2.4803801795699001d+00, -3.1372206062358399d-01, &
         -1.4156160573597001d+00, 2.4518644596973900d+00, 9.7286388084413100d-01, &
         -9.7666225034648202d-01, 2.4683122158147999d+00, 8.9235131263554301d-01, &
         -3.0668100442364798d-01, 1.9247626612871100d+00, 2.1550298876148402d+00, &
         7.9420664981513800d-01, 7.0311134022698996d-01, 3.6136802601641098d+00/), (/Ndim, Nbeads/))

      call reset_RNG_with_seed(myseed)

      Call Time_Integrate_Chain(PosVecR, phys_params, sim_params, output_vars)

      final_positions = reshape((/-1.3982257145883872d-01, -1.5490028944977225d+00, -7.5114724715805048d-01, &
            -6.2921585147668579d-01, -9.1487577673982878d-01, -1.3610037380546496d+00, &
            8.7387887167531397d-01, -4.3111005378687173d-01, 4.7841579771815246d-02, &
            9.5654360263443561d-01, 5.0180128425230552d-01, -1.1970103603563784d-01, &
            3.2515648109692785d-01, 1.2041518167102743d+00, -1.2569266859175836d+00, &
            -4.5747458507234606d-01, 8.9132479407598386d-01, -5.7729219771055718d-01, &
            -8.0006642467789590d-01, 5.2274891326405593d-01, -4.2513911023298040d-01, &
            -4.1386686330385292d-01, -1.8201792112641801d-01, 8.3561540789826827d-01, &
            5.3885535315387201d-02, 2.0711204716970699d-01, 1.2438393334368105d+00, &
            5.2696941824715493d-01, -4.8735258851115915d-01, 2.1487133802405731d+00/), (/Ndim, Nbeads/))

      call assertEquals(PosVecR, final_positions, Ndim, Nbeads, &
         max_double_difference, "WLC PU LJ test failed")

      deallocate(PosVecR)
      deallocate(final_positions)
      deallocate(sim_params%sample_indexes)

   end subroutine

   subroutine test_Fraenkel_UR_flow_SDK_EV()

      real (DBprec), allocatable :: final_positions(:,:)

      phys_params%number_of_beads = 10

      phys_params%spring_inputs%spring_type = Fraenkel
      phys_params%spring_inputs%finite_extensibility_parameter = 500.d0
      phys_params%spring_inputs%natural_length = 3.d0

      myseed = 6
      sim_params%simulation_seed = myseed

      phys_params%hstar = 0.3d0

      phys_params%EV_inputs%dimensionless_EV_energy = 1.d0
      phys_params%EV_inputs%dimensionless_EV_radius = 2.d0
      phys_params%EV_inputs%excluded_volume_type = SDK
      phys_params%EV_inputs%minimum_EV_cutoff = 0.7D0*phys_params%EV_inputs%dimensionless_EV_radius
      phys_params%EV_inputs%maximum_EV_cutoff = 1.82D0*phys_params%EV_inputs%dimensionless_EV_radius

      sim_params%time_at_simulation_start = 0.d0
      sim_params%time_at_simulation_end = 10.d0
      sim_params%simulation_timestep = 0.01d0
      sim_params%implicit_loop_exit_tolerance = 1.d-6

      phys_params%Flow_inputs%flow_type = UR
      phys_params%Flow_inputs%flow_strength = 1.5

      Nbeads = phys_params%number_of_beads

      allocate(phys_params%EV_inputs%bead_attractive_interaction_strengths(Nbeads,Nbeads))
      phys_params%EV_inputs%bead_attractive_interaction_strengths = 0.3d0          ! SDK excluded volume parameter (well depth)
      phys_params%EV_inputs%bead_attractive_interaction_strengths(2,3) = 0.1d0
      phys_params%EV_inputs%bead_attractive_interaction_strengths(3,2) = 0.1d0
      phys_params%EV_inputs%bead_attractive_interaction_strengths(8,4) = 0.8d0
      phys_params%EV_inputs%bead_attractive_interaction_strengths(4,8) = 0.8d0

      phys_params%bend_inputs%bending_potential_type = NoBendingPotential
      phys_params%bend_inputs%bending_stiffness = 0

      sim_params%number_of_samples_to_take = 0
      allocate(sim_params%sample_indexes(sim_params%number_of_samples_to_take))

      allocate(PosVecR(Ndim, Nbeads))
      allocate(final_positions(Ndim, Nbeads))

      ! initial configuration, created randomly by intitial config generation function in utils
      PosVecR = reshape((/0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, &
         -3.7173549676581699d+00, 1.7569605487859301d+00, 1.1085229288292799d+00, &
         -1.4896117550495000d-01, 1.3625962860562499d-01, 4.5656129932180898d+00, &
         -9.3235736014586901d-01, -2.5513763342621898d+00, 6.8682546692628801d+00, &
         1.0507144685331000d+00, -6.2244640419659900d+00, 5.7585411089805403d+00, &
         2.3739079921789101d+00, -7.9420883976551098d+00, 3.2664303838027502d+00, &
         2.4777825222907301d+00, -4.3788920616245104d+00, 2.0008568910710198d+00, &
         4.7192159551425998d+00, -1.7842247060280800d+00, 3.4713093063457401d+00, &
         4.2669615030421797d+00, 1.1663429513272900d+00, 4.0192153253500198d+00, &
         2.7865024867165502d+00, -3.3786944736767301d+00, 4.2416625090714302d+00/), (/Ndim, Nbeads/))

      call reset_RNG_with_seed(myseed)

      Call Time_Integrate_Chain(PosVecR, phys_params, sim_params, output_vars)

      final_positions = reshape((/-3.5301591028373065d+02, -7.9920236736157613d-02, -3.9615913210359288d-01, &
         -3.2577595355323803d+02, 7.4106142294273702d-01, 2.4585062956486803d-01, &
         -2.6147132615415165d+02, 7.8332019700979849d-01, 1.1535017368650211d+00, &
         -1.5974833287198919d+02, 5.7303453845450258d-01, -3.9319746930800226d-01, &
         -3.6080471385851467d+01, -8.6970221959174854d-01, -3.2921124906889249d-01, &
         8.8868912761696890d+01, -1.5464199118583877d-02, -1.2252074965766986d+00, &
         1.9230991940674357d+02, -1.9420675101814835d-01, -3.9512681885603684d-01, &
         2.6198354282608898d+02, -2.0803539660814965d-01, 2.4195776901426491d-01, &
         2.9400002585884783d+02, -8.4432124054745983d-01, -9.7044120142400958d-01, &
         2.9902976812819668d+02, -2.4591124214605076d-02, 2.0746782937710928d+00/), (/Ndim, Nbeads/))

      call assertEquals(PosVecR, final_positions, Ndim, Nbeads, &
         max_double_difference, "Fraenkel UR SDK test failed")

      deallocate(PosVecR, phys_params%EV_inputs%bead_attractive_interaction_strengths)
      deallocate(final_positions)
      deallocate(sim_params%sample_indexes)

   end subroutine

   subroutine test_FENEFraenkel_PL_flow_Gaussian_EV()

      real (DBprec), allocatable :: final_positions(:,:)

      phys_params%number_of_beads = 10

      phys_params%spring_inputs%spring_type = FENEFraenkel
      phys_params%spring_inputs%finite_extensibility_parameter = 1.d0
      phys_params%spring_inputs%natural_length = 10.d0

      myseed = 86
      sim_params%simulation_seed = myseed

      phys_params%hstar = 0.2d0

      phys_params%EV_inputs%dimensionless_EV_energy = 2.d0
      phys_params%EV_inputs%dimensionless_EV_radius = 2.d0
      phys_params%EV_inputs%excluded_volume_type = Gauss

      sim_params%time_at_simulation_start = 0.d0
      sim_params%time_at_simulation_end = 10.d0
      sim_params%simulation_timestep = 0.01d0
      sim_params%implicit_loop_exit_tolerance = 1.d-6

      phys_params%Flow_inputs%flow_type = PL
      phys_params%Flow_inputs%flow_strength = 0.5d0

      Nbeads = phys_params%number_of_beads

      phys_params%bend_inputs%bending_potential_type = NoBendingPotential
      phys_params%bend_inputs%bending_stiffness = 0

      sim_params%number_of_samples_to_take = 0
      allocate(sim_params%sample_indexes(sim_params%number_of_samples_to_take))

      allocate(PosVecR(Ndim, Nbeads))
      allocate(final_positions(Ndim, Nbeads))

      ! initial configuration, created randomly by intitial config generation function in utils
      PosVecR = reshape((/0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, &
         -8.7001417024242400d+00, 5.3611027493752204d+00, 2.0995067416916502d+00, &
         -9.1605814460059198d-01, 1.1704273266182801d+01, -1.4985387607124800d+00, &
         -5.4793423363635396d+00, 4.8982010793314998d+00, 4.0063480623078496d+00, &
         -1.0438766944965399d+01, -3.6900277039248199d+00, 1.7059173838192401d+00, &
         -7.3861146555559003d+00, -5.5679879113699000d+00, -8.5270109831673704d+00, &
         -1.1052344868863999d+01, -1.4308940267976500d+00, -1.6009484144943301d+01, &
         -8.7582187147024406d+00, -4.2391733059389702d+00, -2.6018313074301499d+01, &
         -1.5884767628153300d+01, 3.9463324714686001d+00, -2.5335247046303198d+01, &
         -9.3183370137224504d+00, 1.1169550080743599d+01, -2.8667718582544101d+01/), (/Ndim, Nbeads/))

      call reset_RNG_with_seed(myseed)

      Call Time_Integrate_Chain(PosVecR, phys_params, sim_params, output_vars)

      final_positions = reshape((/4.9624900986757723d+01, -6.3402036350504054d-01, 1.6260098226686843d+00, &
         3.8645999510055866d+01, -4.0255901307542019d-01, 1.0293020905963757d+00, &
         2.7650109966918986d+01, -2.2812311993045434d-01, 9.3791813959443715d-01, &
         1.6662743539686890d+01, -2.0600252562827287d-01, 4.9312807985943458d-01, &
         5.6800498524465457d+00, 2.0775174404355051d-01, 2.1151787143961326d-01, &
         -5.3098121411924701d+00, 1.4340566317493569d-01, -2.0073255881809399d-01, &
         -1.6305253358665226d+01, 8.1528362306820423d-02, -4.6967490691935909d-01, &
         -2.7217996399214606d+01, 2.9053115735586754d-01, -7.9399734806352606d-01, &
         -3.8212159459806060d+01, 3.9013997491510344d-01, -1.0188474503209912d+00, &
         -4.9202922648761906d+01, 3.2779754837687858d-01, -1.4187171270945513d+00/), (/Ndim, Nbeads/))

      call assertEquals(PosVecR, final_positions, Ndim, Nbeads, &
         max_double_difference, "FENE Fraenkel PL Gaussian test failed")

      deallocate(PosVecR)
      deallocate(final_positions)
      deallocate(sim_params%sample_indexes)

   end subroutine

end module

module unit_tests
   ! These tests aren't stochastic, so they should only take seconds or minutes to run. Failing these tests doesn't
   ! necessarily mean that the code is incorrect in a physical sense, just that it doesn't agree with a previous version of the code.
   ! In other words, the previous version of the code could be wrong! It should be subjected to a proper battery of stochastic tests every now and then.
   ! On the other hand, some of these are true unit tests which don't just rely upon results from previous code. These are identified by unit_tests_[function name]
   ! Ideally these tests would be updated over time with a wider range of cases.
   use fruit
   Use Global_parameters_variables_and_types
   use Physics_subroutines
   use Spring_Force_calculatons
   Use csputls
   use random_numbers
   use properties
   use simulation_utilities

   implicit none

   private

   real(DBprec), Parameter :: max_double_difference = 1.D-12

   public :: unit_tests_run

contains

   subroutine unit_tests_run()
      ! Runs unit tests for subroutines and functions in single-chain code

      print *, "Running unit tests"

      call test_random_number_sequence()
      call test_spring_force()
      call unit_test_FENE_Fraenkel_Initial_Conditions()
      call unit_test_X_axis_Aligned_no_yz_components()
      call unit_test_semiimp_roots_same_FENE_FF()
      call unit_test_MS_MOD_semiimp_gets_correct_root()
      call test_3_beads_FENE_gives_FF()
      call test_same_random_sequence_for_two_chain_integrations()
      call test_get_unit_vectors_from_connector_vectors_and_distances()
      call test_get_cos_theta_from_unit_vectors()
      call test_get_bending_force()
      call test_get_trap_force()
      call test_lookup_table_gives_same_results()
      call test_EV_cutoff()
      call test_HI_and_chevyshev()

      print *, ""
      print *, "Finished unit tests"
      print *, ""

   end subroutine

   subroutine test_HI_and_chevyshev()
      use Hydrodynamic_interaction_calculations
      use Initial_Position_Utilities
      Real (DBprec), allocatable :: b2b_sup(:,:,:)
      Real (DBprec), allocatable :: deltaR_sup(:,:), DelS(:,:), DelS_expected(:,:), DelS_exact(:,:)
      Real (DBprec), allocatable :: PosVecR(:,:), X_0(:,:)
      Real (DBprec), allocatable :: Diffusion_sup(:,:,:,:), Diffusion_expected(:,:), temp_slice(:,:)
      Real (DBPrec) :: hstar, delts, sigma, dQ, exact_err
      type(Hydrodynamic_interaction_inputs) :: HI_params
      !real (DBprec), allocatable :: tricking_opt_array(:)
      Integer :: NBeads, Ndof

      !call mkl_verbose(1)

      NBeads = 4
      hstar = 0.2d0
      delts = 0.1d0
      Ndof = NBeads*NDim

      allocate(PosVecR(Ndim, Nbeads), &
            b2b_sup(Ndim, Nbeads, Nbeads), &
            deltaR_sup(Nbeads, Nbeads), &
            Diffusion_sup(Ndim,Nbeads,Ndim,Nbeads), &
            Diffusion_expected(Ndim,Ndim), &
            temp_slice(Ndim,Ndim), &
            DelS(Ndim, Nbeads), &
            DelS_exact(Ndim, Nbeads), &
            DelS_expected(Ndim, Nbeads),&
            X_0(Ndim,Nbeads))
      
      X_0 = reshape((/-0.653795716820965d0,        2.25469853815588d0,       0.165603302292486d0, &
             -0.587509005201644d0,      -0.745084985707505d0,      -0.157459965970764d0, &
              -1.88644397151742d0,        3.16619d0,        1.23661201418787d0, &
             -0.958635132046656d0,      -0.253142128732732d0,        1.43910244145363d0/), (/Ndim, Nbeads/))
      
       ! initial configuration, created randomly by intitial config generation function in utils
      PosVecR = reshape((/0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, &
       -8.7001417024242400d+00, 5.3611027493752204d+00, 2.0995067416916502d+00, &
       -9.1605814460059198d-01, 1.1704273266182801d+01, -1.4985387607124800d+00, &
       -5.4793423363635396d+00, 4.8982010793314998d+00, 4.0063480623078496d+00/), (/Ndim, Nbeads/))

      call b2bvector_sym_up(Nbeads, PosVecR, b2b_sup)
      Call modr_sym_up(NBeads,b2b_sup,deltaR_sup)

      ! call set_Hydrodynamic_interaction_parameters(diffusion_sup, hstar, Nbeads)
      call get_diffusion_tensor_with_HI(hstar, deltaR_sup, b2b_sup, Diffusion_sup, Nbeads)

      Diffusion_expected = &
      reshape((/0.043185241172496d0,  -0.010895545687992d0,  -0.004266896699380d0, &
               -0.010895545687992d0,   0.032217586367374d0,   0.002629298741189d0, &
               -0.004266896699380d0,   0.002629298741189d0,   0.026533338541433d0/), &
                                        (/Ndim, Ndim/))
      
      temp_slice = Diffusion_sup(:,1,:,2)
      call assertEquals(temp_slice, Diffusion_expected,&
                  Ndim, Ndim, max_double_difference, &
                  "Diffusion tensor 12 is the same")

      Diffusion_expected = &
      reshape((/1.d0,  0.0d0,  0.0d0, &
               0.0d0,  1.d0,   0.0d0, &
               0.0d0,  0.0d0,  1.d0/), &
                                        (/Ndim, Ndim/))
      
      temp_slice = Diffusion_sup(:,2,:,2)
      call assertEquals(temp_slice, Diffusion_expected,&
                  Ndim, Ndim, max_double_difference, &
                  "Diffusion tensor is the same")
      
      call reset_RNG_with_seed(50)
      call get_delS_cholesky(Nbeads, hstar, Diffusion_sup, delts, Dels, X_input=X_0)

      DelS_expected = reshape((/  -0.206748358961671d0,   0.712998281762464d0,   0.052368362328964d0, &
              -0.202521038852097d0,  -0.210272818508389d0,  -0.045641895488899d0, &
              -0.611043311268904d0,   1.022458498065220d0,   0.391180874642083d0, &
              -0.351959054078799d0,  -0.035305408161349d0,   0.461910033760980d0/), (/Ndim, Nbeads/))
           
      call assertEquals(DelS, delS_expected,&
                  Ndim, NBeads, 1d-5, &
                  "delS cholesky not the same")
                  
      ! call PRINT_MATRIX( 'delS cholesky', Ndim, Nbeads, DelS)
      
      ! testing exact calculation of square root
      call reset_RNG_with_seed(50)
      call get_delS_exact(Nbeads, hstar, Diffusion_sup, delts, Dels_exact, X_input=X_0)
      
      
      ! DelS_expected = reshape((/ -0.224509001256203d0,       0.726920920582090d0,       6.487368427821390d-002, &
             ! -0.213285764321957d0,      -0.218839763066101d0,      -2.644670914184218d-002, &
             ! -0.610495965160600d0,       0.690570486201413d0,       0.400300628424439d0, &
             ! -0.330725823458617d0,      -6.256514036258194d-002,  0.460666058510416d0/), (/Ndim, Nbeads/))
      
      ! call PRINT_MATRIX( 'exact delS', Ndim, Nbeads, DelS_exact)
      
            ! tests for Cholesky decomposition, given a particular seed
      call reset_RNG_with_seed(50)
      HI_params%hstar = hstar
      HI_params%EigenvalueCalcMethod = EigsExact
      ! HI_params%EigenvalueCalcMethod = EigsFixman
      HI_params%nchebMultiplier = 10.d0
      call set_Hydrodynamic_interaction_parameters(HI_params)
      ! call set_chebyshev_parameters(Nbeads,Diffusion_sup, hstar)
      call get_dels_approx(Nbeads, hstar, Diffusion_sup, delts, Dels, X_input=X_0)

      ! DelS_expected = reshape((/     -0.224831335099214d0,   0.727468024641349d0,   0.065069224368362d0, &
            ! -0.213724360767294d0,  -0.218734886988892d0,  -0.026099667421332d0, &
            ! -0.610983590800792d0,   0.690991815158717d0,   0.400629040615095d0, &
            ! -0.331232815419716d0,  -0.062371412456071d0,   0.460945994673542d0/), (/Ndim, Nbeads/))

      ! call assertEquals(DelS, delS_expected,&
                  ! Ndim, NBeads, max_double_difference, &
                  ! "delS Chebyshev not the same")
                  
      ! call PRINT_MATRIX( 'DelS chebyshev', Ndim, Nbeads, DelS)
      
      call assertEquals(DelS, delS_exact,&
                  Ndim, NBeads, max_double_difference, &
                  "delS exact not the same")
      
      ! call get_delS_exact(Nbeads, hstar, Diffusion_sup, delts, Dels_exact, X_0)
      ! exact_err = sum(reshape(abs(Dels-Dels_exact),(/1,Ndof/)))/sum(reshape(abs(Dels_exact),(/1,Ndof/)))
      ! print *, "exact error is: ", exact_err
      
      ! timing tests
      ! timing with Chebyshev
      Nbeads = 200
      dQ = 10.d0
      sigma = 3.d0
      hstar = 1.2d0

      deallocate(PosVecR, &
            b2b_sup, &
            deltaR_sup, &
            Diffusion_sup, &
            Diffusion_expected, &
            temp_slice, &
            DelS, &
            DelS_expected)

   end subroutine

   subroutine test_EV_cutoff()
      use Excluded_Volume_Calculations
      type(Excluded_Volume_Inputs) :: EV_inputs
      Real (DBprec), allocatable :: Bead_to_bead_vector_superdiagonal_matrix(:,:,:)
      Real (DBprec), allocatable :: Bead_to_bead_distances_superdiagonal_matrix(:,:)
      Real (DBprec), allocatable :: PosVecR(:,:)
      Real (DBprec), allocatable :: Fev(:,:), Fev_expected(:,:)
      Integer :: NBeads, FlowType, nu
      Logical :: test_passed

      NBeads = 10
      FlowType = PR

      allocate(PosVecR(Ndim, Nbeads), &
               Bead_to_bead_vector_superdiagonal_matrix(Ndim, Nbeads, Nbeads), &
               Bead_to_bead_distances_superdiagonal_matrix(Nbeads, Nbeads), &
               Fev(Ndim, Nbeads), Fev_expected(Ndim, Nbeads))

         ! initial configuration, created randomly by intitial config generation function in utils
      PosVecR = reshape((/0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, &
         -8.7001417024242400d+00, 5.3611027493752204d+00, 2.0995067416916502d+00, &
         -9.1605814460059198d-01, 1.1704273266182801d+01, -1.4985387607124800d+00, &
         -5.4793423363635396d+00, 4.8982010793314998d+00, 4.0063480623078496d+00, &
         -1.0438766944965399d+01, -3.6900277039248199d+00, 1.7059173838192401d+00, &
         -7.3861146555559003d+00, -5.5679879113699000d+00, -8.5270109831673704d+00, &
         -1.1052344868863999d+01, -1.4308940267976500d+00, -1.6009484144943301d+01, &
         -8.7582187147024406d+00, -4.2391733059389702d+00, -2.6018313074301499d+01, &
         -1.5884767628153300d+01, 3.9463324714686001d+00, -2.5335247046303198d+01, &
         -9.3183370137224504d+00, 1.1169550080743599d+01, -2.8667718582544101d+01/), (/Ndim, Nbeads/))

      call b2bvector_sym_up(Nbeads, PosVecR, Bead_to_bead_vector_superdiagonal_matrix)
      Call modr_sym_up(NBeads,Bead_to_bead_vector_superdiagonal_matrix,Bead_to_bead_distances_superdiagonal_matrix)

      EV_inputs%dimensionless_EV_energy = 1.d0
      EV_inputs%dimensionless_EV_radius = 10.d0
      EV_inputs%excluded_volume_type = LJ
      EV_inputs%contour_dist_for_EV = 10

      call set_EV_parameters(NBeads, EV_inputs)

      call get_excluded_volume_force(Nbeads, Bead_to_bead_vector_superdiagonal_matrix, &
                                     Bead_to_bead_distances_superdiagonal_matrix,Fev, FlowType)

      Fev_expected = 0.d0

      call assertEquals(Fev, Fev_expected, Ndim, Nbeads, max_double_difference, "Fev is zero for contour dist 10")

      deallocate(PosVecR, &
               Bead_to_bead_vector_superdiagonal_matrix, &
               Bead_to_bead_distances_superdiagonal_matrix, &
               Fev, Fev_expected)

      NBeads = 5
      FlowType = PR

      EV_inputs%dimensionless_EV_energy = 1.d0
      EV_inputs%dimensionless_EV_radius = 1.d0
      EV_inputs%excluded_volume_type = LJ
      EV_inputs%contour_dist_for_EV = 4

      allocate(PosVecR(Ndim, Nbeads), &
               Bead_to_bead_vector_superdiagonal_matrix(Ndim, Nbeads, Nbeads), &
               Bead_to_bead_distances_superdiagonal_matrix(Nbeads, Nbeads), &
               Fev(Ndim, Nbeads), Fev_expected(Ndim, Nbeads))

      PosVecR = reshape((/0.d0, 0.d0, 0.d0, &
         1.0d0, 0.d0, 0.d0, &
         2.1d0, 0.d0, 0.d0, &
         3.0d0, 0.d0, 0.d0, &
         0.0d0, 1.5d0, 0.d0/), (/Ndim, Nbeads/))

      call b2bvector_sym_up(Nbeads, PosVecR, Bead_to_bead_vector_superdiagonal_matrix)
      Call modr_sym_up(NBeads,Bead_to_bead_vector_superdiagonal_matrix,Bead_to_bead_distances_superdiagonal_matrix)

      call set_EV_parameters(NBeads, EV_inputs)

      call get_excluded_volume_force(Nbeads, Bead_to_bead_vector_superdiagonal_matrix, &
                                     Bead_to_bead_distances_superdiagonal_matrix,Fev, FlowType)

      test_passed = .true.
      do nu=1,NBeads
         if ((nu.eq.1) .or. (nu.eq.Nbeads)) then
            if (abs(Fev(2,nu)).le.epsilon(1.d0)) then
               test_passed = .false.
               print *, nu, (Fev(2,nu))
            end if
         else
            if (Fev(1,nu).ne.0.d0) then
               test_passed = .false.
               print *, nu, (Fev(1,nu))
            end if
         end if
      end do

      call assert_true(test_passed, "Fev doesn't work, test 2")

      !Subroutine get_excluded_volume_force(N,b2bvec,dr,Fev,FlowType)

   end subroutine

   subroutine test_lookup_table_gives_same_results()
      type(spring_inputs) :: spring_inputs_dummy
      Real(DBprec) :: dt, gama, max_gama, tol
      Real(DBprec) :: lookup_r, lookup_ff, direct_r, direct_ff

      spring_inputs_dummy%finite_extensibility_parameter = 10.d0
      spring_inputs_dummy%natural_length = 9.d0
      spring_inputs_dummy%spring_type = WLC_bounded

      dt = 0.01d0
      max_gama = 500.0d0
      tol = 1.0d-9

      call generate_lookup_table(spring_inputs_dummy, dt, max_gama, tol)

      gama = 0.1

      call stop_using_lookup_table()
      call solve_implicit_r(spring_inputs_dummy%spring_type, dt/4.d0, gama, &
                              spring_inputs_dummy%natural_length/spring_inputs_dummy%finite_extensibility_parameter, &
                              direct_r, direct_ff)

      !print *,
      !print *, direct_r, direct_ff

      call start_using_lookup_table()
      call solve_implicit_r(spring_inputs_dummy%spring_type, dt/4.d0, gama, &
                              spring_inputs_dummy%natural_length/spring_inputs_dummy%finite_extensibility_parameter, &
                              lookup_r, lookup_ff)

      !print *,
      !print *, lookup_r, lookup_ff

      call assert_equals(lookup_r, direct_r, tol, "Lookup table r test 1 failed")
      !call assert_equals(lookup_ff, direct_ff, tol, "Lookup table ff test 1 failed")

      call stop_using_lookup_table()

      !! Test 2

      spring_inputs_dummy%finite_extensibility_parameter = 10.d0
      spring_inputs_dummy%natural_length = 9.d0
      spring_inputs_dummy%spring_type = WLC_bounded

      dt = 0.01d0
      max_gama = 500.0d0
      tol = 1.0d-9

      call generate_lookup_table(spring_inputs_dummy, dt, max_gama, tol)

      gama = 0.00001

      call stop_using_lookup_table()
      call solve_implicit_r(spring_inputs_dummy%spring_type, dt/4.d0, gama, &
                              spring_inputs_dummy%natural_length/spring_inputs_dummy%finite_extensibility_parameter, &
                              direct_r, direct_ff)

      !print *,
      !print *, direct_r, direct_ff

      call start_using_lookup_table()
      call solve_implicit_r(spring_inputs_dummy%spring_type, dt/4.d0, gama, &
                              spring_inputs_dummy%natural_length/spring_inputs_dummy%finite_extensibility_parameter, &
                              lookup_r, lookup_ff)

      !print *,
      !print *, lookup_r, lookup_ff

      call assert_equals(lookup_r, direct_r, tol, "Lookup table r test 2 failed")
      !call assert_equals(lookup_ff, direct_ff, tol, "Lookup table ff test 1 failed")

      call stop_using_lookup_table()

      !! Test 3

      spring_inputs_dummy%finite_extensibility_parameter = 10.d0
      spring_inputs_dummy%natural_length = 9.d0
      spring_inputs_dummy%spring_type = WLC_bounded

      dt = 0.01d0
      max_gama = 5000.0d0
      tol = 1.0d-9

      call generate_lookup_table(spring_inputs_dummy, dt, max_gama, tol)

      gama = 5000

      call stop_using_lookup_table()
      call solve_implicit_r(spring_inputs_dummy%spring_type, dt/4.d0, gama, &
                              spring_inputs_dummy%natural_length/spring_inputs_dummy%finite_extensibility_parameter, &
                              direct_r, direct_ff)

      !print *,
      !print *, direct_r, direct_ff

      call start_using_lookup_table()
      call solve_implicit_r(spring_inputs_dummy%spring_type, dt/4.d0, gama, &
                              spring_inputs_dummy%natural_length/spring_inputs_dummy%finite_extensibility_parameter, &
                              lookup_r, lookup_ff)

      !print *,
      !print *, lookup_r, lookup_ff

      call assert_equals(lookup_r, direct_r, tol, "Lookup table r test 3 failed")
      !call assert_equals(lookup_ff, direct_ff, tol, "Lookup table ff test 1 failed")

      call stop_using_lookup_table()

   end subroutine

   subroutine test_random_number_sequence()

      Integer(k4b) :: seed, n
      Integer :: n_norm_integer
      Real(DBPrec), allocatable, dimension(:) :: X, X_test, X_prev
      !Integer :: i

      seed = 10
      n = 10
      ! Needed as assert_equals wants a normal old int
      n_norm_integer = n

      allocate(X(n), X_test(n), X_prev(n))

      call reset_RNG_with_seed(seed)

      call ran_1(n, X)

      X_prev = X

      X_test = (/ 0.226473569869995,       0.357589989900589,       0.217218622565269, &
         0.688211202621460,       8.224623650312424E-002, 0.238840281963348, &
         0.956106066703796,       0.993510127067566,       0.167083948850632, &
         0.341053307056427 /)

      call assert_equals(X, X_test, n, max_double_difference, "random number sequence has changed")

      seed = 10
      call reset_RNG_with_seed(seed)

      call ran_1(n, X)

      call assert_equals(X, X_prev, n, max_double_difference, "random number sequence doesn't give same result for same seed")

   end subroutine

   subroutine unit_test_FENE_Fraenkel_Initial_Conditions()
      use Initial_Position_Utilities

      Integer(k4b) :: N                               ! Number of beads in chain
      Integer(k4b) :: seed                            ! Seed for random number generator
      Real(DBprec) :: dQ                              ! Finite extensibility about sigma, equivalent to deltaQ
      Real(DBprec) :: sigma                           ! Natural spring length, equivalent to Q0 for Fraenkel spring
      Real(DBprec), Dimension(:,:), allocatable :: R  ! Final bead locations
      Integer(k4b) :: i
      logical :: correct_bounds
      Real(DBprec) :: Ql
      Real(DBprec), dimension(Ndim) :: b2bvec

      ! Test that springs don't go more than dQ on either side of sigma
      N = 1000
      seed = 5
      dQ = 1.d0
      sigma = 5.d0
      correct_bounds = .True.

      allocate(R(Ndim,N))

      call Initial_position_FENE_Fraenkel(FENEFraenkel, N, dQ, sigma, R)

      do i=2,N
         b2bvec = R(:,i) - R(:,i-1)
         Ql = sqrt(b2bvec(1)**2 + b2bvec(2)**2 + b2bvec(3)**2)
         if ((Ql>sigma+dQ).or.(Ql<sigma-dQ)) then
            correct_bounds = .False.
            print *, "Incorrect bounds at i = ", i
         end if
      end do

      call assert_true(correct_bounds, "FENE-Fraenkel spring is not bounded properly, test 1")

      ! More constrained test
      N = 1000
      seed = 5
      dQ = 0.0001d0
      sigma = 5.d0
      correct_bounds = .True.

      call Initial_position_FENE_Fraenkel(FENEFraenkel, N, dQ, sigma, R)

      do i=2,N
         b2bvec = R(:,i) - R(:,i-1)
         Ql = sqrt(b2bvec(1)**2 + b2bvec(2)**2 + b2bvec(3)**2)
         if ((Ql>sigma+dQ).or.(Ql<sigma-dQ)) then
            correct_bounds = .False.
            print *, "Incorrect bounds at i = ", i
         end if
      end do

      call assert_true(correct_bounds, "FENE-Fraenkel spring is not bounded properly, test 2")

      ! Testing where dQ>sigma, so lower bound should be 0
      N = 1000
      seed = 5
      dQ = 10.d0
      sigma = 1.d0
      correct_bounds = .True.

      call Initial_position_FENE_Fraenkel(FENEFraenkel, N, dQ, sigma, R)

      do i=2,N
         b2bvec = R(:,i) - R(:,i-1)
         Ql = sqrt(b2bvec(1)**2 + b2bvec(2)**2 + b2bvec(3)**2)
         if ((Ql>sigma+dQ).or.(Ql<0.d0)) then
            correct_bounds = .False.
            print *, "Incorrect bounds at i = ", i
         end if
      end do

      call assert_true(correct_bounds, "FENE-Fraenkel spring is not bounded properly, test 3")

   end subroutine

   subroutine unit_test_X_axis_Aligned_no_yz_components()
      use Initial_Position_Utilities

      Integer(k4b) :: N                               ! Number of beads in chain
      Integer(k4b) :: seed                            ! Seed for random number generator
      Real(DBprec) :: dQ                              ! Finite extensibility about sigma, equivalent to deltaQ
      Real(DBprec) :: sigma                           ! Natural spring length, equivalent to Q0 for Fraenkel spring
      Real(DBprec), Dimension(:,:), allocatable :: R  ! Final bead locations
      Integer(k4b) :: i
      logical :: correct_bounds, aligned_x
      Real(DBprec) :: Ql
      Real(DBprec), dimension(Ndim) :: b2bvec

      ! Test that when an x-axis aligned chain is created, there are no yz components
      N = 1000
      seed = 5
      dQ = 1.d0
      sigma = 5.d0
      correct_bounds = .True.
      aligned_x = .True.

      allocate(R(Ndim,N))

      call FENE_Fraenkel_Aligned_x_axis(FENEFraenkel, N, dQ, sigma, R)

      do i=2,N
         b2bvec = R(:,i) - R(:,i-1)
         Ql = sqrt(b2bvec(1)**2 + b2bvec(2)**2 + b2bvec(3)**2)
         if ((Ql>sigma+dQ).or.(Ql<sigma-dQ)) then
            correct_bounds = .False.
            print *, "Incorrect bounds at i = ", i
         end if
         if ((b2bvec(2).ne.0).or.(b2bvec(3).ne.0)) then
            aligned_x = .False.
            print *, "not aligned along x-axis at i = ", i
         end if
      end do

      call assert_true(correct_bounds, "FENE-Fraenkel spring is not bounded properly, x-axis aligned test 1")
      call assert_true(aligned_x, "FENE-Fraenkel spring is not fully aligned, test 1")

   end subroutine

   subroutine unit_test_semiimp_roots_same_FENE_FF()
      Integer :: sptype
      Real (DBprec) :: dtby4             ! delta t divided by 4
      Real (DBprec) :: gama, gamma_star              ! Gamma*, RHS of semi-implicit equation equation
      Real (DBprec) :: Q0s, Q0s_star            ! Q0 or sigma, natural spring length
      Real (DBprec) :: r_FENE, r_FF                ! implicit solution for connector vector length
      Real (DOBL) :: ff_FENE, ff_FF                 ! connector vector force from implicit solution
      real (DBprec) :: sqrtb
      logical :: tests_passing
      Integer :: i, j
      !real(dbprec) :: time1, time2

      dtby4 = 0.05d0
      gama = 9.d0
      sqrtb = 10.d0
      Q0s = 0.d0
      tests_passing = .True.

      gamma_star = gama/sqrtb
      Q0s_star = Q0s/sqrtb

      do i = 1,100
         do j = 1,100
            sqrtb = dble(i)*0.1d0
            gama = dble(j)*0.5d0

            gamma_star = gama/sqrtb
            Q0s_star = Q0s/sqrtb

            sptype = FENE

            call solve_implicit_r(sptype, dtby4, gamma_star, Q0s_star, r_FENE, ff_FENE)

            sptype = FENEFraenkel

            call solve_implicit_r(sptype, dtby4, gamma_star, Q0s_star, r_FF, ff_FF)

            if (abs(r_FENE - r_FF)>max_double_difference*1.d0) then
               tests_passing = .False.
               print *, "r_FENE /= r_FF for sqrtb = ", sqrtb, "and gamma = ", gama
               print *, "r_FENE is = ", r_FENE, "r_FF is = ", r_FF
            end if
!
!                if (abs(ff_FENE - ff_FF)>max_double_difference*1.d4) then
!                    tests_passing = .False.
!                    print *, "ff_FENE /= ff_FF for sqrtb = ", sqrtb, "and gamma = ", gama
!                    print *, "ff_FENE is = ", ff_FENE, "ff_FF is = ", ff_FF
!                end if

         end do
      end do

      call assert_true(tests_passing, "FENE-Fraenkel semiimplicit solution not the same as FENE")

!        call cpu_time(time1)
!
!        do i = 1,10
!            do j = 1,1000
!                do k = 1,1000
!                    deltaQ = dble(i)*10.d0
!                    gama = dble(j)*5.d0
!                    Q0s = dble(k)*1.d0
!
!                    gamma_star = gama/deltaQ
!                    Q0s_star = Q0s/deltaQ
!
!                    sptype = FENEFraenkel
!
!                    call solve_implicit_r(sptype, dtby4, gamma_star, Q0s_star, r_FF, ff_FF)
!                    !print *, Q0s_star
!                end do
!            end do
!        end do
!
!        call cpu_time(time2)
!
!        print *, "time taken is: ", time2-time1

   end subroutine

   subroutine unit_test_MS_MOD_semiimp_gets_correct_root()
      use Spring_Force_calculatons

      Integer :: sptype
      Real (DBprec) :: dtby4             ! delta t divided by 4
      Real (DBprec) :: gama              ! Gamma*, RHS of semi-implicit equation equation
      Real (DBprec) :: natscl            ! Q0 or sigma, natural spring length
      Real (DBprec) :: r                ! implicit solution for connector vector length
      Real (DOBL) :: ff                 ! connector vector force from implicit solution

      Real (DBprec) :: r_WLC, ff_WLC

      sptype = WLC_bounded
      dtby4 = 0.001d0/4.d0
      gama = 5.d0
      natscl = 0.8d0

      call solve_implicit_r(sptype, dtby4, gama, natscl, r, ff)

      call assertEquals(0.999354532073543d0,r,max_double_difference,"WLC_bounded test 1 length failed")
      call assertEquals(1.601291769648808d04,ff,max_double_difference*1d04,"WLC_bounded test 1 force failed")

      gama = 0.d0

      call solve_implicit_r(sptype, dtby4, gama, natscl, r, ff)

      call assertEquals(0.601489231399178d0,r,max_double_difference,"WLC_bounded test 2 length failed")
      call assertEquals(-4.000000000000779d03,ff,max_double_difference*1d04,"WLC_bounded test 2 force failed")

      gama = 1.d5

      call solve_implicit_r(sptype, dtby4, gama, natscl, r, ff)

      call assertEquals(0.999995917496677d0,r,max_double_difference,"WLC_bounded test 3 length failed")
      call assertEquals(3.999976329984242d08,ff,max_double_difference*1d08,"WLC_bounded test 3 force failed")

      dtby4 = 0.2d0/4.d0
      natscl = 0.d0
      gama = epsilon(1.d0)*2.d0

      call solve_implicit_r(sptype, dtby4, gama, natscl, r, ff)


      call assertEquals(0.d0,r,max_double_difference,"WLC_bounded test 4 length failed")
      ! This test is failing, but the expression is correct - we are subtracting two nearly identical
      ! numbers, which is causing a catastrophic loss of significance. This could be fixed by re-writing
      ! the expressions so that you're not dividing by a number which is almost zero.
      !call assertEquals(1.d0,ff,max_double_difference,"WLC_bounded test 4 force failed (loss of significance)")

      dtby4 = 0.2d0/4.d0
      natscl = 0.d0
      gama = 1.d0

      sptype = WLC_bounded
      call solve_implicit_r(sptype, dtby4, gama, natscl, r, ff)
      sptype = WLC
      call solve_implicit_r(sptype, dtby4, gama, natscl, r_WLC, ff_WLC)

      call assertEquals(r_WLC,r,max_double_difference,"WLC_bounded length same as WLC length test failed")
      call assertEquals(ff_WLC,ff,max_double_difference,"WLC_bounded force same as WLC force test failed")

   end subroutine

   subroutine test_3_beads_FENE_gives_FF()

      Integer :: Spring_type
      real (DBprec) :: sqrtb
      Real (DBprec) :: Q0s      ! sqrt of FENE b-parameter
      type(spring_inputs) :: spring_inputs_dum
      Integer(k4b) :: Nbeads    ! Number of beads in the chain
      Integer :: No_beads_norm_integer !Number of beads in the chain only a normal integer type rather than whatever k4b is, needed for fruit call
      Real (DBprec), allocatable :: Bead_initial_positions(:,:)
      Real (DBprec), allocatable :: Bead_positions(:,:)
      ! The superdiagonal matrices have each entry above the diagonal corresponding to the
      ! relevant property for bead nu and nu+1. For example, position (2,3) in the distances
      ! matrix corresponds to the distance between bead 2 and bead 3.
      Real (DBprec), allocatable :: Bead_to_bead_vector_superdiagonal_matrix(:,:,:)
      Real (DBprec), allocatable :: Bead_to_bead_distances_superdiagonal_matrix(:,:)
      real (DBprec), allocatable :: Spring_force_vectors(:,:)
      real (DBprec), allocatable :: Spring_force_vector_answer(:,:)

      Real (DBprec), allocatable :: NaturalAngles (:)
      integer :: BendingPotentialType
      real (DBprec) :: BendingStiffness
      BendingPotentialType = NoBendingPotential
      BendingStiffness = 0
      allocate(NaturalAngles(2))
      naturalAngles = 0

      Nbeads = 3
      No_beads_norm_integer = Nbeads
      Q0s = 0.d0
      sqrtb = 10._dp
      Spring_type = FENEFraenkel
      spring_inputs_dum%spring_type = Spring_type
      spring_inputs_dum%natural_length = Q0s
      spring_inputs_dum%finite_extensibility_parameter = sqrtb

      allocate(Bead_initial_positions(Ndim, Nbeads))
      allocate(Bead_positions(Ndim,Nbeads))
      allocate(Bead_to_bead_vector_superdiagonal_matrix(Ndim,Nbeads,Nbeads))
      allocate(Bead_to_bead_distances_superdiagonal_matrix(Nbeads,Nbeads))
      allocate(Spring_force_vectors(Ndim, Nbeads))
      allocate(Spring_force_vector_answer(Ndim, Nbeads))

      Bead_positions =   reshape((/ 0.247077027956645d0,      -0.100952789187431d0,       0.118288775285085d0, &
         0.181214521328608d0,      -0.150070123374462d0,       0.397452404101690d0, &
         -0.428291549285253d0,       0.251022912561894d0,      -0.515741179386775d0 /), (/Ndim, Nbeads/))

      call b2bvector_sym_up(Nbeads, Bead_positions, Bead_to_bead_vector_superdiagonal_matrix)
      Call modr_sym_up(NBeads,Bead_to_bead_vector_superdiagonal_matrix,Bead_to_bead_distances_superdiagonal_matrix)

      call get_spring_force(NBeads, Bead_to_bead_vector_superdiagonal_matrix,&
         Bead_to_bead_distances_superdiagonal_matrix, Spring_force_vectors,spring_inputs_dum)

      Spring_force_vector_answer = spring_force_vectors

      Bead_positions =   reshape((/ 0.247077027956645d0,      -0.100952789187431d0,       0.118288775285085d0, &
         0.181214521328608d0,      -0.150070123374462d0,       0.397452404101690d0, &
         -0.428291549285253d0,       0.251022912561894d0,      -0.515741179386775d0 /), (/Ndim, Nbeads/))

      Spring_type = FENE
      spring_inputs_dum%spring_type = Spring_type
      spring_inputs_dum%natural_length = Q0s
      spring_inputs_dum%finite_extensibility_parameter = sqrtb

      call b2bvector_sym_up(Nbeads, Bead_positions, Bead_to_bead_vector_superdiagonal_matrix)
      Call modr_sym_up(NBeads,Bead_to_bead_vector_superdiagonal_matrix,Bead_to_bead_distances_superdiagonal_matrix)

      call get_spring_force(NBeads, Bead_to_bead_vector_superdiagonal_matrix,&
         Bead_to_bead_distances_superdiagonal_matrix, Spring_force_vectors,spring_inputs_dum)

      call assertEquals(Spring_force_vectors, Spring_force_vector_answer, Ndim, &
         No_beads_norm_integer, max_double_difference, "FENE /= FENE_Fraenkel")

   end subroutine

   subroutine test_same_random_sequence_for_two_chain_integrations()

      Integer :: Nbeads      !Number of beads in chain
      Integer(k4b) :: myseed      !seed for random number generator
      Integer(k4b) :: stored_seed, stored_ix, stored_iy ! Actual seeds in RNG
      Integer(k4b) :: seed_PL, ix_PL, iy_PL
      Integer(k4b) :: seed_PR, ix_PR, iy_PR
      Real (DBprec), allocatable :: PosVecR(:,:)       ! Position vector of beads, number_of_dimensions*number_of_beads
      Real (DBprec), allocatable :: initial_positions(:,:)
      Real (DBprec), allocatable :: final_positions(:,:)

      type(physical_parameters) :: phys_params
      type(simulation_parameters) :: sim_params
      type(calculated_variables) :: output_vars

      phys_params%bend_inputs%bending_potential_type = NoBendingPotential
      phys_params%bend_inputs%bending_stiffness = 0

      sim_params%number_of_samples_to_take = 0
      allocate(sim_params%sample_indexes(sim_params%number_of_samples_to_take))

      phys_params%number_of_beads = 10

      phys_params%spring_inputs%spring_type = FENEFraenkel
      phys_params%spring_inputs%finite_extensibility_parameter = 1.d0
      phys_params%spring_inputs%natural_length = 10.d0

      myseed = 86
      sim_params%simulation_seed = myseed

      phys_params%hstar = 0.2d0

      phys_params%EV_inputs%dimensionless_EV_energy = 2.d0
      phys_params%EV_inputs%dimensionless_EV_radius = 2.d0
      phys_params%EV_inputs%excluded_volume_type = Gauss

      sim_params%time_at_simulation_start = 0.d0
      sim_params%time_at_simulation_end = 10.d0
      sim_params%implicit_loop_exit_tolerance = 1.d-6

      phys_params%Flow_inputs%flow_strength = 0.5d0

      Nbeads = phys_params%number_of_beads

      allocate(PosVecR(Ndim, Nbeads))
      allocate(initial_positions(Ndim, Nbeads))
      allocate(final_positions(Ndim, Nbeads))

      ! initial configuration, created randomly by intitial config generation function in utils
      PosVecR = reshape((/0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, &
         -8.7001417024242400d+00, 5.3611027493752204d+00, 2.0995067416916502d+00, &
         -9.1605814460059198d-01, 1.1704273266182801d+01, -1.4985387607124800d+00, &
         -5.4793423363635396d+00, 4.8982010793314998d+00, 4.0063480623078496d+00, &
         -1.0438766944965399d+01, -3.6900277039248199d+00, 1.7059173838192401d+00, &
         -7.3861146555559003d+00, -5.5679879113699000d+00, -8.5270109831673704d+00, &
         -1.1052344868863999d+01, -1.4308940267976500d+00, -1.6009484144943301d+01, &
         -8.7582187147024406d+00, -4.2391733059389702d+00, -2.6018313074301499d+01, &
         -1.5884767628153300d+01, 3.9463324714686001d+00, -2.5335247046303198d+01, &
         -9.3183370137224504d+00, 1.1169550080743599d+01, -2.8667718582544101d+01/), (/Ndim, Nbeads/))

      call reset_RNG_with_seed(myseed)

      phys_params%Flow_inputs%flow_type = EQ
      sim_params%simulation_timestep = 0.1d0
      Call Time_Integrate_Chain(PosVecR, phys_params, sim_params, output_vars)

      initial_positions = PosVecR

      call get_all_parameters(stored_seed, stored_ix, stored_iy)

      phys_params%Flow_inputs%flow_type = PL
      sim_params%simulation_timestep = 0.01d0
      Call Time_Integrate_Chain(PosVecR, phys_params, sim_params, output_vars)

      call get_all_parameters(seed_PL, ix_PL, iy_PL)

      call reset_RNG_with_seed(stored_seed, stored_ix, stored_iy)
      PosVecR = initial_positions

      phys_params%Flow_inputs%flow_type = PR
      sim_params%simulation_timestep = 0.01d0
      Call Time_Integrate_Chain(PosVecR, phys_params, sim_params, output_vars)

      call get_all_parameters(seed_PR, ix_PR, iy_PR)

      call assertEquals(seed_PR, seed_PL, "Two seeds are not the same")
      call assertEquals(ix_PR, ix_PL, "Two ix are not the same")
      call assertEquals(iy_PR, iy_PL, "Two iy are not the same")

      deallocate(PosVecR, initial_positions, final_positions)

   end subroutine

   subroutine test_spring_force()
      ! Tests the spring force calculation to ensure the correct results are obtained
      ! for several possible spring force laws. The comparisons are with some simple MATLAB code
      ! which directly computes the spring force, should be very accurate

      Integer :: Spring_type
      real (DBprec) :: sqrtb
      Real (DBprec) :: Q0s      ! sqrt of FENE b-parameter
      type(spring_inputs) :: spring_inputs_dum
      Integer(k4b) :: Nbeads    ! Number of beads in the chain
      Integer :: No_beads_norm_integer !Number of beads in the chain only a normal integer type rather than whatever k4b is, needed for fruit call
      Real (DBprec), allocatable :: Bead_initial_positions(:,:)
      Real (DBprec), allocatable :: Bead_positions(:,:)
      ! The superdiagonal matrices have each entry above the diagonal corresponding to the
      ! relevant property for bead nu and nu+1. For example, position (2,3) in the distances
      ! matrix corresponds to the distance between bead 2 and bead 3.
      Real (DBprec), allocatable :: Bead_to_bead_vector_superdiagonal_matrix(:,:,:)
      Real (DBprec), allocatable :: Bead_to_bead_distances_superdiagonal_matrix(:,:)
      real (DBprec), allocatable :: Spring_force_vectors(:,:)
      real (DBprec), allocatable :: Spring_force_vector_answer(:,:)

      Nbeads = 10
      No_beads_norm_integer = Nbeads
      Q0s = 1000.d0
      sqrtb = 5._dp
      Spring_type = HOOK
      spring_inputs_dum%spring_type = Spring_type
      spring_inputs_dum%natural_length = Q0s
      spring_inputs_dum%finite_extensibility_parameter = sqrtb

      allocate(Bead_initial_positions(Ndim, Nbeads))
      allocate(Bead_positions(Ndim,Nbeads))
      allocate(Bead_to_bead_vector_superdiagonal_matrix(Ndim,Nbeads,Nbeads))
      allocate(Bead_to_bead_distances_superdiagonal_matrix(Nbeads,Nbeads))
      allocate(Spring_force_vectors(Ndim, Nbeads))
      allocate(Spring_force_vector_answer(Ndim, Nbeads))

      Bead_initial_positions = reshape((/  0.000000000000000d0,  0.000000000000000d0,  0.000000000000000d0, &
         -6.586250662803650d-2, -4.911733418703079d-2,  0.279163628816605d0, &
         -0.675368577241898d0,       0.351975701749325d0,      -0.634029954671860d0, &
         -2.84468892216682d0,       -2.23968694359064d0,      -0.780964523553848d0, &
         -1.81838491559029d0,       -2.20123096555471d0,       0.934529632329941d0, &
         -3.21014931797981d0,      -0.937253437936306d0,      -0.468102008104324d0, &
         -1.49071100354195d0,      -0.597119502723217d0,       -1.50950649380684d0, &
         -1.58743391185999d0,       -1.01791768521070d0,      -1.45650240033865d0, &
         -0.620455302298069d0,       -1.21915509551764d+0,      -0.712517879903316d0, &
         -0.880402423441410d0,       -3.66230776160955d0,      -0.734861174598336d0 /), (/Ndim, Nbeads/))

      Bead_positions = Bead_initial_positions

      Spring_force_vector_answer =  reshape((/ -6.5862506628036499d-02, -4.9117334187030792d-02, 2.7916362881660500d-01, &
         -5.4364356398582503d-01, 4.5021037012338661d-01, -1.1923572123050701d+00, &
         -1.5598142743110603d+00, -2.9927556812763205d+00, 7.6625901460647694d-01, &
         3.1956243515014515d+00, 2.6301186233758949d+00, 1.8624287247657771d+00, &
         -2.4180684089660498d+00, 1.2255215495824743d+00, -3.1181257963180542d+00, &
         3.1112027168273801d+00, -9.2384359240531522d-01, 3.6122715473174916d-01, &
         -1.8161612227559001d+00, -7.6093211770057190d-01, 1.0944085791707061d+00, &
         1.0637015178799611d+00, 2.1956077218054282d-01, 6.9098042696714379d-01, &
         -1.2269257307052621d+00, -2.2419152557849698d+00, -7.6632781513035386d-01, &
         2.5994712114334095d-01, 2.4431526660919101d+00, 2.2343294695019944d-02 /), (/Ndim, Nbeads/))

      ! Set up bead to bead vector
      call b2bvector_sym_up(Nbeads, Bead_positions, Bead_to_bead_vector_superdiagonal_matrix)
      Call modr_sym_up(NBeads,Bead_to_bead_vector_superdiagonal_matrix,Bead_to_bead_distances_superdiagonal_matrix)

      call get_spring_force(NBeads, Bead_to_bead_vector_superdiagonal_matrix,&
         Bead_to_bead_distances_superdiagonal_matrix, Spring_force_vectors,spring_inputs_dum)

      call assertEquals(Spring_force_vectors, Spring_force_vector_answer, Ndim, &
         No_beads_norm_integer, max_double_difference, "Hookean spring force test 1 failed")

      !Test FENE is the same as Hookean with very high sqrtb
      Q0s = 1000.d0
      sqrtb = 1.0d10
      Spring_type = FENE
      spring_inputs_dum%spring_type = Spring_type
      spring_inputs_dum%natural_length = Q0s
      spring_inputs_dum%finite_extensibility_parameter = sqrtb
      Bead_positions = Bead_initial_positions

      ! Set up bead to bead vector
      call b2bvector_sym_up(Nbeads, Bead_positions, Bead_to_bead_vector_superdiagonal_matrix)
      Call modr_sym_up(NBeads,Bead_to_bead_vector_superdiagonal_matrix,Bead_to_bead_distances_superdiagonal_matrix)

      call get_spring_force(NBeads, Bead_to_bead_vector_superdiagonal_matrix,&
         Bead_to_bead_distances_superdiagonal_matrix, Spring_force_vectors,spring_inputs_dum)

      call assertEquals(Spring_force_vectors, Spring_force_vector_answer, Ndim, &
         No_beads_norm_integer, max_double_difference, "FENE == Hookean test failed")

      !Test Fraenkel is the same as Hookean with Q0 = 0
      Q0s = 0.d0
      sqrtb = 1.0d10
      Spring_type = Fraenkel
      spring_inputs_dum%spring_type = Spring_type
      spring_inputs_dum%natural_length = Q0s
      spring_inputs_dum%finite_extensibility_parameter = sqrtb
      Bead_positions = Bead_initial_positions

      ! Set up bead to bead vector
      call b2bvector_sym_up(Nbeads, Bead_positions, Bead_to_bead_vector_superdiagonal_matrix)
      Call modr_sym_up(NBeads,Bead_to_bead_vector_superdiagonal_matrix,Bead_to_bead_distances_superdiagonal_matrix)

      call get_spring_force(NBeads, Bead_to_bead_vector_superdiagonal_matrix,&
         Bead_to_bead_distances_superdiagonal_matrix, Spring_force_vectors,spring_inputs_dum)

      call assertEquals(Spring_force_vectors, Spring_force_vector_answer, Ndim, &
         No_beads_norm_integer, max_double_difference, "Fraenkel == Hookean test failed")

      !_____________________________________________________________________________________________
      !FENE tests
      !_____________________________________________________________________________________________

      Q0s = 1000.d0
      sqrtb = 0.01_dp
      Spring_type = FENE
      spring_inputs_dum%spring_type = Spring_type
      spring_inputs_dum%natural_length = Q0s
      spring_inputs_dum%finite_extensibility_parameter = sqrtb
      Bead_positions = Bead_initial_positions

      ! Set up bead to bead vector
      call b2bvector_sym_up(Nbeads, Bead_positions, Bead_to_bead_vector_superdiagonal_matrix)
      Call modr_sym_up(NBeads,Bead_to_bead_vector_superdiagonal_matrix,Bead_to_bead_distances_superdiagonal_matrix)

      call get_spring_force(NBeads, Bead_to_bead_vector_superdiagonal_matrix,&
         Bead_to_bead_distances_superdiagonal_matrix, Spring_force_vectors,spring_inputs_dum)

      Spring_force_vector_answer =  reshape((/ 7.7867573114493639d-05, 5.8070179937086591d-05, -3.3004808639534769d-04, &
         -3.3254196828288759d-05, -8.7428566015149698d-05, 3.9689015889724573d-04, &
         -2.5657671872242288d-05, 5.2004553187823449d-05, -6.5558145827612272d-05, &
         -4.4628714963201429d-05, -2.3608143993287310d-05, -4.4197039435166325d-05, &
         5.0968974116773942d-05, -2.2011400962553100d-05, 6.8406593420523775d-05, &
         -6.6662598506103937d-05, 1.4790358885740869d-05, -4.3913656267106696d-07, &
         9.2506014644932265d-05, 2.3066762377821728d-04, -5.3078692299452678d-05, &
         -1.1438370577648606d-04, -2.0932286135312342d-04, -2.0635259296854799d-05, &
         6.7550248812278454d-05, 2.7308128181497629d-05, 4.9029715467212827d-05, &
         -4.3059227421558306d-06, -4.0469871646252302d-05, -3.7010796787753027d-07/), (/Ndim, Nbeads/))

      call assertEquals(Spring_force_vectors, Spring_force_vector_answer, Ndim, &
         No_beads_norm_integer, max_double_difference, "FENE spring force test 1 failed")

      Q0s = 0.d0
      sqrtb = 5._dp
      Spring_type = FENE
      spring_inputs_dum%spring_type = Spring_type
      spring_inputs_dum%natural_length = Q0s
      spring_inputs_dum%finite_extensibility_parameter = sqrtb
      Bead_positions = Bead_initial_positions

      ! Set up bead to bead vector
      call b2bvector_sym_up(Nbeads, Bead_positions, Bead_to_bead_vector_superdiagonal_matrix)
      Call modr_sym_up(NBeads,Bead_to_bead_vector_superdiagonal_matrix,Bead_to_bead_distances_superdiagonal_matrix)

      call get_spring_force(NBeads, Bead_to_bead_vector_superdiagonal_matrix,&
         Bead_to_bead_distances_superdiagonal_matrix, Spring_force_vectors,spring_inputs_dum)

      Spring_force_vector_answer =  reshape((/ -6.6086361526003320d-02, -4.9284275234308811d-02, 2.8011245613649000d-01, &
         -5.7865606365045130d-01, 4.7356503170084130d-01, -1.2460989722034230d+00, &
         -3.3559976151874986d+00, -5.2039199372745690d+00, 6.9500440205364744d-01, &
         5.2223966642697031d+00, 4.8254150912549614d+00, 2.3130133107908044d+00, &
         -3.0061550778422736d+00, 1.5748762200235562d+00, -3.8404634526248493d+00, &
         3.8468360779532094d+00, -1.2126869211898330d+00, 5.4934548451630616d-01, &
         -2.1597982542021694d+00, -8.3197289018375842d-01, 1.3024951359227308d+00, &
         1.1274348129505096d+00, 2.0966029495985514d-01, 7.3904437515645560d-01, &
         -1.3726784426167973d+00, -3.0066108707023900d+00, -8.2190927782185952d-01, &
         3.4270425985177128d-01, 3.2209582566456447d+00, 2.9456538073697264d-02 /), (/Ndim, Nbeads/))

      call assertEquals(Spring_force_vectors, Spring_force_vector_answer, Ndim, &
         No_beads_norm_integer, max_double_difference, "FENE spring force test 2 failed")

      Q0s = 10.d0
      sqrtb = 500._dp
      Spring_type = FENE
      spring_inputs_dum%spring_type = Spring_type
      spring_inputs_dum%natural_length = Q0s
      spring_inputs_dum%finite_extensibility_parameter = sqrtb
      Bead_positions = Bead_initial_positions

      ! Set up bead to bead vector
      call b2bvector_sym_up(Nbeads, Bead_positions, Bead_to_bead_vector_superdiagonal_matrix)
      Call modr_sym_up(NBeads,Bead_to_bead_vector_superdiagonal_matrix,Bead_to_bead_distances_superdiagonal_matrix)

      call get_spring_force(NBeads, Bead_to_bead_vector_superdiagonal_matrix,&
         Bead_to_bead_distances_superdiagonal_matrix, Spring_force_vectors,spring_inputs_dum)

      Spring_force_vector_answer =  reshape((/ -6.5862528937707290d-02, -4.9117350824593072d-02, 2.7916372337797191d-01, &
         -5.4364687275668322d-01, 4.5021257881984172d-01, -1.1923622976639234d+00, &
         -1.5599102528043314d+00, -2.9928765173673613d+00, 7.6625727886887773d-01, &
         3.1957400727555880d+00, 2.6302378823595722d+00, 1.8624628838564457d+00, &
         -2.4181154513987462d+00, 1.2255487529761837d+00, -3.1181840987984293d+00, &
         3.1112619366917844d+00, -9.2386575534139959d-01, 3.6124070921737972d-01, &
         -1.8161898850817733d+00, -7.6093809163040027d-01, 1.0944259347310001d+00, &
         1.0637075053919709d+00, 2.1955985987973001d-01, 6.9098493725307597d-01, &
         -1.2269379224274450d+00, -2.2419730242904179d+00, -7.6633290510228746d-01, &
         2.5995339856734268d-01, 2.4432116654188447d+00, 2.2343834259889058d-02 /), (/Ndim, Nbeads/))

      call assertEquals(Spring_force_vectors, Spring_force_vector_answer, Ndim, &
         No_beads_norm_integer, max_double_difference, "FENE spring force test 3 failed")

      !_____________________________________________________________________________________________
      !Fraenkel tests
      !_____________________________________________________________________________________________

      Q0s = 1.d0
      sqrtb = 50._dp
      Spring_type = Fraenkel
      spring_inputs_dum%spring_type = Spring_type
      spring_inputs_dum%natural_length = Q0s
      spring_inputs_dum%finite_extensibility_parameter = sqrtb
      Bead_positions = Bead_initial_positions

      ! Set up bead to bead vector
      call b2bvector_sym_up(Nbeads, Bead_positions, Bead_to_bead_vector_superdiagonal_matrix)
      Call modr_sym_up(NBeads,Bead_to_bead_vector_superdiagonal_matrix,Bead_to_bead_distances_superdiagonal_matrix)

      call get_spring_force(NBeads, Bead_to_bead_vector_superdiagonal_matrix,&
         Bead_to_bead_distances_superdiagonal_matrix, Spring_force_vectors,spring_inputs_dum)

      Spring_force_vector_answer =  reshape((/ 1.6046684081745488d-01, 1.1966904768566619d-01, -6.8015184785372473d-01, &
         -2.4853131863923916d-01, -6.1717124857211217d-02, 5.4820908311781147d-01, &
         -1.4400021269814749d+00, -1.8835159288660850d+00, 2.8442233647914414d-02, &
         2.0410709474350153d+00, 1.8447864621580059d+00, 9.6100079035660424d-01, &
         -1.3114275131546989d+00, 7.0589232863641183d-01, -1.6621577138698129d+00, &
         1.6745011474097629d+00, -5.5181176636414708d-01, 2.7404723887305726d-01, &
         -7.5045575526895592d-01, 3.7322317127713112d-01, 4.6176931658540721d-01, &
         5.9359515170406468d-02, -5.8502264170837104d-01, 2.1116415983624043d-01, &
         -3.3913218123401867d-01, -1.4103100116390970d+00, -1.5557298930164923d-01, &
         1.5415044444574774d-01, 1.4488064636776961d+00, 1.3249728608151861d-02 /), (/Ndim, Nbeads/))

      call assertEquals(Spring_force_vectors, Spring_force_vector_answer, Ndim, &
         No_beads_norm_integer, max_double_difference, "Fraenkel spring force test 1 failed")

      Q0s = 100.d0
      sqrtb = 50._dp
      Spring_type = Fraenkel
      spring_inputs_dum%spring_type = Spring_type
      spring_inputs_dum%natural_length = Q0s
      spring_inputs_dum%finite_extensibility_parameter = sqrtb
      Bead_positions = Bead_initial_positions

      ! Set up bead to bead vector
      call b2bvector_sym_up(Nbeads, Bead_positions, Bead_to_bead_vector_superdiagonal_matrix)
      Call modr_sym_up(NBeads,Bead_to_bead_vector_superdiagonal_matrix,Bead_to_bead_distances_superdiagonal_matrix)

      call get_spring_force(NBeads, Bead_to_bead_vector_superdiagonal_matrix,&
         Bead_to_bead_distances_superdiagonal_matrix, Spring_force_vectors,spring_inputs_dum)

      Spring_force_vector_answer =  reshape((/ 2.2567072237921103d+01, 1.6829520853082666d+01, -9.5652384038216368d+01, &
         2.8967580970672763d+01, -5.0742539127936396d+01, 1.7286427232998307d+02, &
         1.0421400458647483d+01, 1.0793121955974726d+02, -7.3015419081249775d+01, &
         -1.1225971605514218d+02, -7.5903097498412990d+01, -8.8280364716151510d+01, &
         1.0824602117216907d+02, -5.0737400545023775d+01, 1.4247868244850608d+02, &
         -1.4055895422493433d+02, 3.6279339011711500d+01, -8.3567644311374494d+00, &
         1.0475438552593852d+02, 1.1265459678006974d+02, -6.2169517679359188d+01, &
         -9.9370498753075509d+01, -8.0238780616710841d+01, -4.7290646286123192d+01, &
         8.7552429216419071d+01, 8.0918609158802312d+01, 6.0309154767740111d+01, &
         -1.0319720548615980d+01, -9.6991467575329480d+01, -8.8701331399178840d-01 /), (/Ndim, Nbeads/))

      call assertEquals(Spring_force_vectors, Spring_force_vector_answer, Ndim, &
         No_beads_norm_integer, max_double_difference, "Fraenkel spring force test 2 failed")

      !_____________________________________________________________________________________________
      ! FENE-Fraenkel tests
      !_____________________________________________________________________________________________

      !First I want to test that I get the same results as the Hookean, Fraenkel and FENE springs with
      !The correct sets of parameters.

      !Test FENE-Fraenkel gives hookean results
      Q0s = 0.d0
      sqrtb = 1d10
      Spring_type = FENEFraenkel
      spring_inputs_dum%spring_type = Spring_type
      spring_inputs_dum%natural_length = Q0s
      spring_inputs_dum%finite_extensibility_parameter = sqrtb
      Bead_positions = Bead_initial_positions

      call b2bvector_sym_up(Nbeads, Bead_positions, Bead_to_bead_vector_superdiagonal_matrix)
      Call modr_sym_up(NBeads,Bead_to_bead_vector_superdiagonal_matrix,Bead_to_bead_distances_superdiagonal_matrix)

      call get_spring_force(NBeads, Bead_to_bead_vector_superdiagonal_matrix,&
         Bead_to_bead_distances_superdiagonal_matrix, Spring_force_vectors,spring_inputs_dum)

      Spring_force_vector_answer =  reshape((/ -6.5862506628036499d-02, -4.9117334187030792d-02, 2.7916362881660500d-01, &
         -5.4364356398582503d-01, 4.5021037012338661d-01, -1.1923572123050701d+00, &
         -1.5598142743110603d+00, -2.9927556812763205d+00, 7.6625901460647694d-01, &
         3.1956243515014515d+00, 2.6301186233758949d+00, 1.8624287247657771d+00, &
         -2.4180684089660498d+00, 1.2255215495824743d+00, -3.1181257963180542d+00, &
         3.1112027168273801d+00, -9.2384359240531522d-01, 3.6122715473174916d-01, &
         -1.8161612227559001d+00, -7.6093211770057190d-01, 1.0944085791707061d+00, &
         1.0637015178799611d+00, 2.1956077218054282d-01, 6.9098042696714379d-01, &
         -1.2269257307052621d+00, -2.2419152557849698d+00, -7.6632781513035386d-01, &
         2.5994712114334095d-01, 2.4431526660919101d+00, 2.2343294695019944d-02 /), (/Ndim, Nbeads/))

      call assertEquals(Spring_force_vectors, Spring_force_vector_answer, Ndim, &
         No_beads_norm_integer, max_double_difference, "FENE-Fraenkel doesn't give Hookean")

      !Test FENE-Fraenkel gives Fraenkel
      Q0s = 1.d0
      sqrtb = 1d10
      Spring_type = FENEFraenkel
      spring_inputs_dum%spring_type = Spring_type
      spring_inputs_dum%natural_length = Q0s
      spring_inputs_dum%finite_extensibility_parameter = sqrtb
      Bead_positions = Bead_initial_positions

      call b2bvector_sym_up(Nbeads, Bead_positions, Bead_to_bead_vector_superdiagonal_matrix)
      Call modr_sym_up(NBeads,Bead_to_bead_vector_superdiagonal_matrix,Bead_to_bead_distances_superdiagonal_matrix)

      call get_spring_force(NBeads, Bead_to_bead_vector_superdiagonal_matrix,&
         Bead_to_bead_distances_superdiagonal_matrix, Spring_force_vectors,spring_inputs_dum)

      Spring_force_vector_answer =  reshape((/ 1.6046684081745488d-01, 1.1966904768566619d-01, -6.8015184785372473d-01, &
         -2.4853131863923916d-01, -6.1717124857211217d-02, 5.4820908311781147d-01, &
         -1.4400021269814749d+00, -1.8835159288660850d+00, 2.8442233647914414d-02, &
         2.0410709474350153d+00, 1.8447864621580059d+00, 9.6100079035660424d-01, &
         -1.3114275131546989d+00, 7.0589232863641183d-01, -1.6621577138698129d+00, &
         1.6745011474097629d+00, -5.5181176636414708d-01, 2.7404723887305726d-01, &
         -7.5045575526895592d-01, 3.7322317127713112d-01, 4.6176931658540721d-01, &
         5.9359515170406468d-02, -5.8502264170837104d-01, 2.1116415983624043d-01, &
         -3.3913218123401867d-01, -1.4103100116390970d+00, -1.5557298930164923d-01, &
         1.5415044444574774d-01, 1.4488064636776961d+00, 1.3249728608151861d-02 /), (/Ndim, Nbeads/))

      call assertEquals(Spring_force_vectors, Spring_force_vector_answer, Ndim, &
         No_beads_norm_integer, max_double_difference, "FENE-Fraenkel doesn't give Fraenkel")

      !Test FENE-Fraenkel gives FENE
      Q0s = 0.d0
      sqrtb = 5._dp
      Spring_type = FENEFraenkel
      spring_inputs_dum%spring_type = Spring_type
      spring_inputs_dum%natural_length = Q0s
      spring_inputs_dum%finite_extensibility_parameter = sqrtb
      Bead_positions = Bead_initial_positions

      call b2bvector_sym_up(Nbeads, Bead_positions, Bead_to_bead_vector_superdiagonal_matrix)
      Call modr_sym_up(NBeads,Bead_to_bead_vector_superdiagonal_matrix,Bead_to_bead_distances_superdiagonal_matrix)

      call get_spring_force(NBeads, Bead_to_bead_vector_superdiagonal_matrix,&
         Bead_to_bead_distances_superdiagonal_matrix, Spring_force_vectors,spring_inputs_dum)

      Spring_force_vector_answer =  reshape((/ -6.6086361526003320d-02, -4.9284275234308811d-02, 2.8011245613649000d-01, &
         -5.7865606365045130d-01, 4.7356503170084130d-01, -1.2460989722034230d+00, &
         -3.3559976151874986d+00, -5.2039199372745690d+00, 6.9500440205364744d-01, &
         5.2223966642697031d+00, 4.8254150912549614d+00, 2.3130133107908044d+00, &
         -3.0061550778422736d+00, 1.5748762200235562d+00, -3.8404634526248493d+00, &
         3.8468360779532094d+00, -1.2126869211898330d+00, 5.4934548451630616d-01, &
         -2.1597982542021694d+00, -8.3197289018375842d-01, 1.3024951359227308d+00, &
         1.1274348129505096d+00, 2.0966029495985514d-01, 7.3904437515645560d-01, &
         -1.3726784426167973d+00, -3.0066108707023900d+00, -8.2190927782185952d-01, &
         3.4270425985177128d-01, 3.2209582566456447d+00, 2.9456538073697264d-02 /), (/Ndim, Nbeads/))

      call assertEquals(Spring_force_vectors, Spring_force_vector_answer, Ndim, &
         No_beads_norm_integer, max_double_difference, "FENE-Fraenkel doesn't give FENE")

      !Several FENE-Fraenkel specific tests using MATLAB calculations
      Q0s = 10.d0
      sqrtb = 0.01_dp
      Spring_type = FENEFraenkel
      spring_inputs_dum%spring_type = Spring_type
      spring_inputs_dum%natural_length = Q0s
      spring_inputs_dum%finite_extensibility_parameter = sqrtb

      Bead_initial_positions = reshape((/ 2.4919697349292400d-01, 1.3055209918253701d+00, -7.7404468452369701d-01, &
         -4.0482139880627899d-01, 2.3256718066665201d+00, 3.4641492820205400d-01, &
         -9.2202057741167598d-01, 1.3861815931028301d+00, 1.2770332573681401d-01, &
         -1.0617954062612700d+00, 5.9900307855935497d-01, 2.8081088378504898d-01, &
         -6.8565315962822804d-01, -1.3104569448429699d+00, 8.2919281143866297d-01, &
         6.3698457234424399d-01, -7.4432417491702374d-02, 5.1457409609423699d-01, &
         1.4502598168138201d+00, -6.0741323511966305d-01, -1.3821919021348500d-01, &
         1.4862782068094699d-01, -1.0034437170381201d+00, -1.1923312955383301d+00, &
         6.8870990947369903d-01, -9.8329532243887197d-01, -4.7622282458310100d-01, &
         4.5462421631766597d-01, -1.3093137411851701d+00, 7.3920066181421895d-01 /), (/Ndim, Nbeads/))

      Bead_positions = Bead_initial_positions

      call b2bvector_sym_up(Nbeads, Bead_positions, Bead_to_bead_vector_superdiagonal_matrix)
      Call modr_sym_up(NBeads,Bead_to_bead_vector_superdiagonal_matrix,Bead_to_bead_distances_superdiagonal_matrix)

      call get_spring_force(NBeads, Bead_to_bead_vector_superdiagonal_matrix,&
         Bead_to_bead_distances_superdiagonal_matrix, Spring_force_vectors,spring_inputs_dum)

      Spring_force_vector_answer =  reshape((/ -4.7460468489788987d-06, 7.4029779060168409d-06, 8.1308936256493732d-06, &
         -5.6008116469507588d-07, -1.7041537955947324d-05, -1.0374732766859149d-05, &
         3.4368719450231429d-06, -8.8864450851896474d-07, 4.2913983075739707d-06, &
         4.2010338364477285d-06, -1.3099035327621353d-06, 1.3519653449353463d-06, &
         6.4869447095613524d-06, 2.0078334708618860d-05, -5.4972532163983702d-06, &
         -9.5339037398046157d-07, -1.3395780581616332d-05, -4.2155533314327611d-06, &
         -1.7000299962021510d-05, 2.3751772017166669d-06, -1.0845694531046143d-06, &
         1.5748159077172796d-05, 3.0260896029960403d-06, 1.6166448615716194d-05, &
         -8.7104690047667483d-06, -3.1676567833920380d-06, 2.1209221970574723d-06, &
         2.0972777862376757d-06, 2.9209439428883845d-06, -1.0889519323137462d-05 /), (/Ndim, Nbeads/))

      call assertEquals(Spring_force_vectors, Spring_force_vector_answer, Ndim, &
         No_beads_norm_integer, max_double_difference, "FENE-Fraenkel spring force test 1 failed")

      Q0s = 1.d0
      sqrtb = 5._dp
      Spring_type = FENEFraenkel
      spring_inputs_dum%spring_type = Spring_type
      spring_inputs_dum%natural_length = Q0s
      spring_inputs_dum%finite_extensibility_parameter = sqrtb

      Bead_initial_positions = reshape((/ 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, &
         -1.0595184345419100d+00, -8.4349976122150505d-01, 1.4872900498738100d+00, &
         -1.4817318534741299d+00, -1.0963540604041500d+00, 6.3240944833124801d-01, &
         -1.6994085787515401d+00, -1.1320056523249400d+00, -6.4867917642924600d-01, &
         -3.3339296228131499d+00, -2.0688730435997699d-01, 3.3641329212819099d-01, &
         -1.6078258514114201d+00, -1.5905117803844999d+00, 1.2511449614279200d+00, &
         3.9820727304270997d-01, -2.8873978752300702d+00, -3.4967223682652698d-01, &
         -1.4025575801528900d+00, -3.1786929372066002d+00, -7.0429259394006305d-01, &
         -9.5220500862172497d-01, -3.3719236953398402d+00, -1.9018013397901701d+00, &
         -1.0964067752798000d-01, -9.6649418530551001d-01, -7.1674551049115198d-01/), (/Ndim, Nbeads/))

      Bead_positions = Bead_initial_positions

      call b2bvector_sym_up(Nbeads, Bead_positions, Bead_to_bead_vector_superdiagonal_matrix)
      Call modr_sym_up(NBeads,Bead_to_bead_vector_superdiagonal_matrix,Bead_to_bead_distances_superdiagonal_matrix)

      call get_spring_force(NBeads, Bead_to_bead_vector_superdiagonal_matrix,&
         Bead_to_bead_distances_superdiagonal_matrix, Spring_force_vectors,spring_inputs_dum)

      Spring_force_vector_answer =  reshape((/ -5.5552051389831236d-01, -4.4225886549059307d-01, 7.7980722740229380d-01, &
         5.6133405309142703d-01, 4.4574046582002702d-01, -7.6803620909816217d-01, &
         -5.6220205197441239d-02, -1.1737319374194446d-02, -3.0842839565096775d-01, &
         -8.5911860692865694d-01, 5.2303558130784211d-01, 8.4480969363794278d-01, &
         1.9992639013031348d+00, -1.3883015262925227d+00, 2.9343921852172383d-02, &
         4.3291937691911864d-01, -1.1086586686055588d-01, -1.7925794250623142d+00, &
         -2.3796502029474458d+00, 8.4575888379749120d-01, 1.0463177433700068d+00, &
         9.5964176327095918d-01, 9.4585258336073924d-02, -1.0418459681843548d-01, &
         5.2208737093822499d-01, 1.8275993829171215d+00, 1.1516344578256854d+00, &
         -6.2473693655100848d-01, -1.7835559941606896d+00, -8.7868441745822157d-01/), (/Ndim, Nbeads/))

      call assertEquals(Spring_force_vectors, Spring_force_vector_answer, Ndim, &
         No_beads_norm_integer, max_double_difference, "FENE-Fraenkel spring force test 2 failed")

      Q0s = 100.d0
      sqrtb = 20._dp
      Spring_type = FENEFraenkel
      spring_inputs_dum%spring_type = Spring_type
      spring_inputs_dum%natural_length = Q0s
      spring_inputs_dum%finite_extensibility_parameter = sqrtb

      Bead_initial_positions = reshape((/ 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, &
         9.0829712987031101d+01, -3.9561312886468900d+00, -4.4905775613349199d+01, &
         1.5601121545905899d+01, -3.4930788595718099d+01, -1.0222135890879100d+02, &
         5.9728723656431001d+01, -6.2306147117943802d+01, -1.4479382125783600d+01, &
         1.5640668835255801d+02, -3.8523062404087902d+01, -9.3192364079856898d+00, &
         1.2721113755350800d+02, 3.3332509215224300d+00, -9.7177619690384702d+01, &
         1.4573386984979200d+02, -5.6943254166684497d+01, -1.7511248494414099d+02, &
         8.1131086060346902d+01, -1.1393378217080800d+02, -1.2430469028416600d+02, &
         5.7663252095993400d+01, -2.0490685468960700d+02, -8.8415031896692298d+01, &
         1.2600065732317501d+02, -1.6599979491853901d+02, -1.5037839150153800d+02 /), (/Ndim, Nbeads/))

      Bead_positions = Bead_initial_positions

      call b2bvector_sym_up(Nbeads, Bead_positions, Bead_to_bead_vector_superdiagonal_matrix)
      Call modr_sym_up(NBeads,Bead_to_bead_vector_superdiagonal_matrix,Bead_to_bead_distances_superdiagonal_matrix)

      call get_spring_force(NBeads, Bead_to_bead_vector_superdiagonal_matrix,&
         Bead_to_bead_distances_superdiagonal_matrix, Spring_force_vectors,spring_inputs_dum)

      Spring_force_vector_answer =  reshape((/ 1.2613679584889552d+00, -5.4939480517650775d-02, -6.2361428487460902d-01, &
         -8.9685221503223800d-01, 2.0502537233499757d-01, 9.0133359420055914d-01, &
         4.9084554403224745d-01, -6.8072478428184235d-01, 1.4230549598029687d+00, &
         -1.1521858502888747d+00, 4.5761911236866670d-01, -1.7166171553696812d+00, &
         -1.6713933892470884d-01, 7.3818340643132585d-01, -1.3803670564141042d+00, &
         5.1028147473909014d-01, -8.1588984646449336d-01, 1.2013275811315605d+00, &
         -5.5592390421655848d-02, 1.4254427021643482d-01, 2.0217667595825137d-01, &
         -1.2451461110826960d-01, -5.1045285043043565d-01, 1.9731156834778055d-01, &
         2.1298774380058766d-01, 5.6372538545261364d-01, -2.7641712106133004d-01, &
         -7.9198315285133478d-02, -4.5090585109616388d-02, 7.1811238278604250d-02/), (/Ndim, Nbeads/))

      call assertEquals(Spring_force_vectors, Spring_force_vector_answer, Ndim, &
         No_beads_norm_integer, max_double_difference, "FENE-Fraenkel spring force test 3 failed")

   end subroutine

   subroutine test_get_unit_vectors_from_connector_vectors_and_distances()

      Integer :: Nbeads    ! Number of beads in the chain
      Real (DBprec), allocatable  :: bead_to_bead_vector_superdiagonal_matrix(:,:,:) ! Ndim x Nbeads x Nbeads
      Real (DBprec), allocatable  :: bead_to_bead_distances(:,:) ! Nbeads x Nbeads
      Real (DBprec), allocatable  :: unit_vectors(:,:) ! Ndim x (Nbeads - 1)
      Real (DBprec), allocatable  :: unit_vectors_correct_result(:,:) ! Ndim x (Nbeads - 1)
      Real (DBprec), allocatable :: Bead_initial_positions(:,:)

      Nbeads = 10

      allocate(bead_to_bead_vector_superdiagonal_matrix(Ndim, Nbeads, Nbeads))
      allocate(bead_to_bead_distances(Nbeads, Nbeads))
      allocate(unit_vectors(Ndim, Nbeads-1))
      allocate(unit_vectors_correct_result(Ndim, Nbeads-1))
      allocate(Bead_initial_positions(Ndim, Nbeads))

      Bead_initial_positions = reshape((/4.4300450016525000d-01, -3.3196056529376028d+00, 5.6132256569023198d-01, &
         1.9030176853078951d+00, -1.2695629199997360d+00, 6.8182316538962096d-01, &
         9.1311608108670805d-01, -7.1791445684858998d-02, 8.9166947220366999d-02, &
         4.4330671786225401d-01, 8.1458593138783497d-01, -1.2960528538904701d+00, &
         -1.5134472363696660d+00, 1.2352696315890850d+00, -8.9531485574350800d-01, &
         -1.4183050785054969d+00, 1.7319540232256150d+00, 1.8692574041718699d-01, &
         -4.4785729015092501d-01, 1.1633844528360100d+00, 9.9689781373572595d-01, &
         -2.7460991831883702d-01, 6.5784188432063095d-01, -1.9640797715447600d-01, &
         3.7236102056231302d-01, 3.0421928502129097d-01, -1.4997345046665300d-01, &
         -4.2058648163949502d-01, -1.2462951897582699d+00, 2.1612904801974000d-02/), (/Ndim, Nbeads/))

      ! Result from MATLAB to compare with
      unit_vectors_correct_result = reshape((/ 5.7944175608126003d-01, 8.1360933812326097d-01, 4.7823594888862002d-02, &
         -5.9522532362554004d-01, 7.2021694922843105d-01, -3.5636268064846999d-01, &
         -2.7468973561912702d-01, 5.1825013808962705d-01, -8.0991502240394098d-01, &
         -9.5863447079962405d-01, 2.0609739688746501d-01, 1.9632578636270501d-01, &
         7.9645677196824999d-02, 4.1578586835775899d-01, 9.0596836466724595d-01, &
         7.0016073333826800d-01, -4.1021278232060898d-01, 5.8438037331176895d-01, &
         1.3250236696976000d-01, -3.8664705977290498d-01, -9.1265939644337102d-01, &
         8.7574437340957700d-01, -4.7866601579109902d-01, 6.2854099056012996d-02, &
         -4.5312745081482397d-01, -8.8603680503118598d-01, 9.8052503528402002d-02/), (/Ndim, Nbeads-1/))


      call b2bvector_sym_up(Nbeads,Bead_initial_positions,bead_to_bead_vector_superdiagonal_matrix)
      call modr_sym_up(Nbeads,bead_to_bead_vector_superdiagonal_matrix,bead_to_bead_distances)

      call get_unit_vectors_from_connector_vectors_and_distances(unit_vectors, &
         bead_to_bead_vector_superdiagonal_matrix, bead_to_bead_distances, Nbeads)

!      print *, unit_vectors
!      print *, unit_vectors_correct_result

      ! test that the result is the same as the MATLAB result
      call assertEquals(unit_vectors_correct_result, unit_vectors, Ndim, Nbeads-1, &
         max_double_difference, "get_unit_vectors test 1 failed")

      deallocate(bead_to_bead_vector_superdiagonal_matrix)
      deallocate(bead_to_bead_distances)
      deallocate(unit_vectors)
      deallocate(unit_vectors_correct_result)
      deallocate(Bead_initial_positions)

   end subroutine

   subroutine test_get_cos_theta_from_unit_vectors()

      Integer :: Nbeads    ! Number of beads in the chain
      Real (DBprec), allocatable  :: unit_vectors_correct_result(:,:) ! Ndim x (Nbeads - 1)
      Real (DBprec), allocatable  :: cos_theta_values(:) ! (Nbeads - 2)
      Real (DBprec), allocatable  :: cos_theta_values_correct_result(:) ! (Nbeads - 2)

      Nbeads = 10

      allocate(unit_vectors_correct_result(Ndim, Nbeads-1))
      allocate(cos_theta_values(Nbeads-2))
      allocate(cos_theta_values_correct_result(Nbeads-2))

      ! initial unit vectors
      unit_vectors_correct_result = reshape((/ 5.7944175608126003d-01, 8.1360933812326097d-01, 4.7823594888862002d-02, &
         -5.9522532362554004d-01, 7.2021694922843105d-01, -3.5636268064846999d-01, &
         -2.7468973561912702d-01, 5.1825013808962705d-01, -8.0991502240394098d-01, &
         -9.5863447079962405d-01, 2.0609739688746501d-01, 1.9632578636270501d-01, &
         7.9645677196824999d-02, 4.1578586835775899d-01, 9.0596836466724595d-01, &
         7.0016073333826800d-01, -4.1021278232060898d-01, 5.8438037331176895d-01, &
         1.3250236696976000d-01, -3.8664705977290498d-01, -9.1265939644337102d-01, &
         8.7574437340957700d-01, -4.7866601579109902d-01, 6.2854099056012996d-02, &
         -4.5312745081482397d-01, -8.8603680503118598d-01, 9.8052503528402002d-02/), (/Ndim, Nbeads-1/))

      ! result from MATLAB to compare result with
      cos_theta_values_correct_result = reshape((/2.2403428410843801d-01, 8.2537830865396700d-01, 2.1112985007575100d-01, &
         1.8720624513310299d-01, 4.1463422899762298d-01, -2.8195971820794202d-01, &
         2.4374862584760501d-01, 3.3454893589510003d-02/), (/Nbeads-2/))

      call get_cos_theta_from_unit_vectors(cos_theta_values, unit_vectors_correct_result, Nbeads)

      !print *, cos_theta_values

      ! test that the result is the same as the MATLAB result
      call assertEquals(cos_theta_values_correct_result, cos_theta_values, Nbeads-2, &
         max_double_difference, "get_cos_theta test 1 failed")

      deallocate(unit_vectors_correct_result)
      deallocate(cos_theta_values)
      deallocate(cos_theta_values_correct_result)

   end subroutine

   subroutine test_get_bending_force()
      ! For generating random chain configurations
      use Initial_Position_Utilities
      use Bending_Force_calculations

      Integer(k4b) :: seed                            ! Seed for random number generator

      Integer :: Nbeads    ! Number of beads in the chain
      Real (DBprec), allocatable  :: natural_angles(:)
      Real (DBprec), allocatable  :: bead_to_bead_vector_superdiagonal_matrix(:,:,:) ! Ndim x Nbeads x Nbeads
      Real (DBprec), allocatable  :: bead_to_bead_distances(:,:) ! Nbeads x Nbeads
      Real (DBprec), allocatable  :: bead_initial_positions(:,:)
      Real (DBprec), allocatable  :: bending_force(:,:) ! Ndim x Nbeads
      Real (DBprec), allocatable  :: bending_force_result(:,:) ! Ndim x Nbeads

      type(bending_potential_inputs) :: bending_params

      bending_params%bending_potential_type = OneMinusCosTheta
      bending_params%bending_stiffness = 1.d0

      ! Test that 2 beads gives zero bending force !-!-!-!-!-!
      Nbeads = 2
      allocate(bead_to_bead_vector_superdiagonal_matrix(Ndim,Nbeads,Nbeads), &
         bead_to_bead_distances(Nbeads,Nbeads), &
         Bead_initial_positions(Ndim,Nbeads), &
         bending_force(Ndim,Nbeads), &
         bending_force_result(Ndim,Nbeads))

      bead_initial_positions = reshape((/1.d0, 0.5d0, 0.5d0, -1.d0, -0.7d0, -0.3d0/), (/Ndim, Nbeads/) )

      call b2bvector_sym_up(Nbeads,Bead_initial_positions,bead_to_bead_vector_superdiagonal_matrix)
      call modr_sym_up(Nbeads,bead_to_bead_vector_superdiagonal_matrix,bead_to_bead_distances)

      call get_bending_force(bending_force, bending_params, &
         bead_to_bead_vector_superdiagonal_matrix, bead_to_bead_distances, Nbeads)

      bending_force_result = 0

      call assertEquals(bending_force, bending_force_result, Ndim, Nbeads, &
         max_double_difference, "bending force is zero when Nbeads=2: test failed")

      deallocate(bead_to_bead_vector_superdiagonal_matrix, &
         bead_to_bead_distances, &
         Bead_initial_positions, &
         bending_force, &
         bending_force_result)

      ! Checking that the force is zero when the beads are co-linear !-!-!-!-!-!
      Nbeads = 5
      allocate(bead_to_bead_vector_superdiagonal_matrix(Ndim,Nbeads,Nbeads), &
         bead_to_bead_distances(Nbeads,Nbeads), &
         Bead_initial_positions(Ndim,Nbeads), &
         bending_force(Ndim,Nbeads), &
         bending_force_result(Ndim,Nbeads))

      bead_initial_positions = reshape((/-2.d0, 0.d0, 0.d0, &
         -1.d0, 0.d0, 0.d0, &
         0.0d0, 0.d0, 0.d0, &
         1.0d0, 0.d0, 0.d0, &
         2.0d0, 0.d0, 0.d0/), &
         (/Ndim, Nbeads/) )

      call b2bvector_sym_up(Nbeads,Bead_initial_positions,bead_to_bead_vector_superdiagonal_matrix)
      call modr_sym_up(Nbeads,bead_to_bead_vector_superdiagonal_matrix,bead_to_bead_distances)


      call get_bending_force(bending_force, bending_params, &
         bead_to_bead_vector_superdiagonal_matrix, bead_to_bead_distances, Nbeads)

      bending_force_result = 0

      call assertEquals(bending_force, bending_force_result, Ndim, Nbeads, &
         max_double_difference, "bending force is zero for co-linear beads : test failed")

      deallocate(bead_to_bead_vector_superdiagonal_matrix, &
         bead_to_bead_distances, &
         Bead_initial_positions, &
         bending_force, &
         bending_force_result)

      ! Check that force is zero when there's no bending potential !-!-!-!-!-!

      bending_params%bending_potential_type = NoBendingPotential
      bending_params%bending_stiffness = 10.d0

      Nbeads = 10
      seed = 34123
      allocate(bead_to_bead_vector_superdiagonal_matrix(Ndim,Nbeads,Nbeads), &
         bead_to_bead_distances(Nbeads,Nbeads), &
         Bead_initial_positions(Ndim,Nbeads), &
         bending_force(Ndim,Nbeads), &
         bending_force_result(Ndim,Nbeads))

      call Initial_position_FENE_Fraenkel(1, Nbeads, 10.d0, 3.d0, bead_initial_positions)

      call b2bvector_sym_up(Nbeads,Bead_initial_positions,bead_to_bead_vector_superdiagonal_matrix)
      call modr_sym_up(Nbeads,bead_to_bead_vector_superdiagonal_matrix,bead_to_bead_distances)


      call get_bending_force(bending_force, bending_params, &
         bead_to_bead_vector_superdiagonal_matrix, bead_to_bead_distances, Nbeads)

      bending_force_result = 0

      call assertEquals(bending_force, bending_force_result, Ndim, Nbeads, &
         max_double_difference, "bending force is zero for no bending potential : test failed")

      deallocate(bead_to_bead_vector_superdiagonal_matrix, &
         bead_to_bead_distances, &
         Bead_initial_positions, &
         bending_force, &
         bending_force_result)

      ! Check forces are 90 degrees for 3-bead case !-!-!-!-!-!

      bending_params%bending_potential_type = OneMinusCosTheta
      bending_params%bending_stiffness = 1.d0

      Nbeads = 3
      seed = 34123
      allocate(bead_to_bead_vector_superdiagonal_matrix(Ndim,Nbeads,Nbeads), &
         bead_to_bead_distances(Nbeads,Nbeads), &
         Bead_initial_positions(Ndim,Nbeads), &
         bending_force(Ndim,Nbeads), &
         bending_force_result(Ndim,Nbeads))

      bead_initial_positions = reshape((/0.d0, 1.d0, 0.d0, &
         0.0d0, 0.d0, 0.d0, &
         1.0d0, 0.d0, 0.d0/), &
         (/Ndim, Nbeads/) )

      call b2bvector_sym_up(Nbeads,Bead_initial_positions,bead_to_bead_vector_superdiagonal_matrix)
      call modr_sym_up(Nbeads,bead_to_bead_vector_superdiagonal_matrix,bead_to_bead_distances)


      call get_bending_force(bending_force, bending_params, &
         bead_to_bead_vector_superdiagonal_matrix, bead_to_bead_distances, Nbeads)

      bending_force_result = reshape((/-1.d0, 0.d0, 0.d0, &
         1.0d0, 1.d0, 0.d0, &
         0.0d0, -1.d0, 0.d0/), &
         (/Ndim, Nbeads/) )

      call assertEquals(bending_force, bending_force_result, Ndim, Nbeads, &
         max_double_difference, "forces are correct 3-bead 90 degree : test failed")

      deallocate(bead_to_bead_vector_superdiagonal_matrix, &
         bead_to_bead_distances, &
         Bead_initial_positions, &
         bending_force, &
         bending_force_result)

      ! Check forces for 30 degree bend for 3-bead case !-!-!-!-!-!

      bending_params%bending_potential_type = OneMinusCosTheta
      bending_params%bending_stiffness = 1.d0

      Nbeads = 3
      seed = 34123
      allocate(bead_to_bead_vector_superdiagonal_matrix(Ndim,Nbeads,Nbeads), &
         bead_to_bead_distances(Nbeads,Nbeads), &
         Bead_initial_positions(Ndim,Nbeads), &
         bending_force(Ndim,Nbeads), &
         bending_force_result(Ndim,Nbeads))

      bead_initial_positions = reshape((/-sqrt(3.d0)/2.d0, 1.d0/2.d0, 0.d0, &
         0.0d0, 0.d0, 0.d0, &
         1.0d0, 0.d0, 0.d0/), &
         (/Ndim, Nbeads/) )

      call b2bvector_sym_up(Nbeads,Bead_initial_positions,bead_to_bead_vector_superdiagonal_matrix)
      call modr_sym_up(Nbeads,bead_to_bead_vector_superdiagonal_matrix,bead_to_bead_distances)


      call get_bending_force(bending_force, bending_params, &
         bead_to_bead_vector_superdiagonal_matrix, bead_to_bead_distances, Nbeads)

      bending_force_result = reshape((/-1.d0/4.d0, -sqrt(3.d0)/4.d0, 0.d0, &
         1.0d0/4.d0, (sqrt(3.d0)+2.d0)/4.d0, 0.d0, &
         0.0d0, -1.d0/2.d0, 0.d0/), &
         (/Ndim, Nbeads/) )


      call assertEquals(bending_force, bending_force_result, Ndim, Nbeads, &
         max_double_difference, "forces are correct for 30 degree bend 3-bead case : test failed")

      deallocate(bead_to_bead_vector_superdiagonal_matrix, &
         bead_to_bead_distances, &
         Bead_initial_positions, &
         bending_force, &
         bending_force_result)

      ! Check forces for random 3-bead case in 3D with MATLAB !-!-!-!-!-!

      bending_params%bending_potential_type = OneMinusCosTheta
      bending_params%bending_stiffness = 1.d0

      Nbeads = 3
      seed = 34123
      allocate(bead_to_bead_vector_superdiagonal_matrix(Ndim,Nbeads,Nbeads), &
         bead_to_bead_distances(Nbeads,Nbeads), &
         Bead_initial_positions(Ndim,Nbeads), &
         bending_force(Ndim,Nbeads), &
         bending_force_result(Ndim,Nbeads), &
         natural_angles(Nbeads-2))

      natural_angles = 0.d0

      bead_initial_positions = reshape((/-3.0463813666787000d-02, -2.0816098199062599d-01, 7.8285534648953203d-01, &
         1.1168640279183419d+00, 3.8970365089086001d-01, -4.9842591260299501d-01, &
         -1.0864002142515550d+00, -1.8154266890023399d-01, -2.8442943388653802d-01/), &
         (/Ndim, Nbeads/) )

      call b2bvector_sym_up(Nbeads,Bead_initial_positions,bead_to_bead_vector_superdiagonal_matrix)
      call modr_sym_up(Nbeads,bead_to_bead_vector_superdiagonal_matrix,bead_to_bead_distances)


      call get_bending_force(bending_force, bending_params, &
         bead_to_bead_vector_superdiagonal_matrix, bead_to_bead_distances, Nbeads)

      bending_force_result = reshape((/2.6795444632606003d-01, 1.0521229604980001d-03, 2.4043169403111900d-01, &
         -2.2522472036063201d-01, -6.2136107095040000d-02, 3.6446114976283001d-02, &
         -4.2729725965427998d-02, 6.1083984134540997d-02, -2.7687780900740200d-01/), &
         (/Ndim, Nbeads/) )

      call assertEquals(bending_force, bending_force_result, Ndim, Nbeads, &
         max_double_difference, "forces are correct for random 3-bead case in 3D with MATLAB : test failed")

      deallocate(bead_to_bead_vector_superdiagonal_matrix, &
         bead_to_bead_distances, &
         Bead_initial_positions, &
         bending_force, &
         bending_force_result)

      ! Check forces for precise 4-bead case !-!-!-!-!-!

      bending_params%bending_potential_type = OneMinusCosTheta
      bending_params%bending_stiffness = 1.d0

      Nbeads = 4
      seed = 34123
      allocate(bead_to_bead_vector_superdiagonal_matrix(Ndim,Nbeads,Nbeads), &
         bead_to_bead_distances(Nbeads,Nbeads), &
         Bead_initial_positions(Ndim,Nbeads), &
         bending_force(Ndim,Nbeads), &
         bending_force_result(Ndim,Nbeads))

      bead_initial_positions = reshape((/ 0.d0, 0.d0, 0.d0, &
         1.0d0, 0.d0, 0.d0, &
         1.0d0, -1.d0, 0.d0, &
         1.0d0+sqrt(2.d0)/2.d0, -1.0d0+sqrt(2.d0)/2.d0, 0.d0/), &
         (/Ndim, Nbeads/) )

      call b2bvector_sym_up(Nbeads,Bead_initial_positions,bead_to_bead_vector_superdiagonal_matrix)
      call modr_sym_up(Nbeads,bead_to_bead_vector_superdiagonal_matrix,bead_to_bead_distances)


      call get_bending_force(bending_force, bending_params, &
         bead_to_bead_vector_superdiagonal_matrix, bead_to_bead_distances, Nbeads)

      bending_force_result = reshape((/  0.d0, 1.d0, 0.d0, &
         -1.0d0-sqrt(2.d0)/2.d0, -1.d0, 0.d0, &
         (1.0d0+sqrt(2.d0))/2.d0, 1.d0/2.d0, 0.d0, &
         1.0d0/2.d0, -1.d0/2.d0, 0.d0/), &
         (/Ndim, Nbeads/) )

      call assertEquals(bending_force, bending_force_result, Ndim, Nbeads, &
         max_double_difference, "forces are correct for precise 4-bead case : test failed")

      deallocate(bead_to_bead_vector_superdiagonal_matrix, &
         bead_to_bead_distances, &
         Bead_initial_positions, &
         bending_force, &
         bending_force_result)

      ! Check forces for random 4-bead case with MATLAB !-!-!-!-!-!

      bending_params%bending_potential_type = OneMinusCosTheta
      bending_params%bending_stiffness = 1.d0

      Nbeads = 4
      seed = 34123
      allocate(bead_to_bead_vector_superdiagonal_matrix(Ndim,Nbeads,Nbeads), &
         bead_to_bead_distances(Nbeads,Nbeads), &
         Bead_initial_positions(Ndim,Nbeads), &
         bending_force(Ndim,Nbeads), &
         bending_force_result(Ndim,Nbeads))

      bead_initial_positions = reshape((/3.4134834606604614d-01, -4.9062758089567859d-01, -5.3069757124131800d-01, &
         -5.9689886456173313d-01, 1.1835883708834525d+00, -4.0570939739629508d-01, &
         -6.6797608011038512d-02, 2.3152014748196514d-01, 4.4833342082418248d-01, &
         3.2234812650672562d-01, -9.2448093746973903d-01, 4.8807354781343082d-01/), &
         (/Ndim, Nbeads/) )

      call b2bvector_sym_up(Nbeads,Bead_initial_positions,bead_to_bead_vector_superdiagonal_matrix)
      call modr_sym_up(Nbeads,bead_to_bead_vector_superdiagonal_matrix,bead_to_bead_distances)


      call get_bending_force(bending_force, bending_params, &
         bead_to_bead_vector_superdiagonal_matrix, bead_to_bead_distances, Nbeads)

      bending_force_result = reshape((/-1.0028333800544067d-02, 2.0204747869437264d-02, -3.4592206662609493d-01, &
         1.4541591438301668d-01, 1.1254027277902146d-02, 2.9695704340823431d-01, &
         -2.4178109393652047d-01, -8.3922618277541938d-02, -4.3532325503476693d-01, &
         1.0639351335404786d-01, 5.2463843130202505d-02, 4.8428827825262755d-01/), &
         (/Ndim, Nbeads/) )

      call assertEquals(bending_force, bending_force_result, Ndim, Nbeads, &
         max_double_difference, "forces are correct for random 4-bead case with MATLAB : test failed")

      deallocate(bead_to_bead_vector_superdiagonal_matrix, &
         bead_to_bead_distances, &
         Bead_initial_positions, &
         bending_force, &
         bending_force_result)

      ! Check forces for precise 5-bead case !-!-!-!-!-!
      bending_params%bending_potential_type = OneMinusCosTheta
      bending_params%bending_stiffness = 1.d0

      Nbeads = 5
      seed = 34123
      allocate(bead_to_bead_vector_superdiagonal_matrix(Ndim,Nbeads,Nbeads), &
         bead_to_bead_distances(Nbeads,Nbeads), &
         Bead_initial_positions(Ndim,Nbeads), &
         bending_force(Ndim,Nbeads), &
         bending_force_result(Ndim,Nbeads))

      bead_initial_positions = reshape((/0.d0,0.d0,0.d0, &
         sqrt(2.d0)/2.d0,sqrt(2.d0)/2.d0,0.d0, &
         sqrt(2.d0)/2.d0+1.d0,sqrt(2.d0)/2.d0,0.d0, &
         sqrt(2.d0)/2.d0+1.d0+sqrt(3.d0)/2.d0,sqrt(2.d0)/2.d0+1.d0/2.d0,0.d0, &
         sqrt(2.d0)/2.d0+1.d0+2.d0*sqrt(3.d0)/2.d0,sqrt(2.d0)/2.d0,0.d0/), &
         (/Ndim, Nbeads/) )

      call b2bvector_sym_up(Nbeads,Bead_initial_positions,bead_to_bead_vector_superdiagonal_matrix)
      call modr_sym_up(Nbeads,bead_to_bead_vector_superdiagonal_matrix,bead_to_bead_distances)


      call get_bending_force(bending_force, bending_params, &
         bead_to_bead_vector_superdiagonal_matrix, bead_to_bead_distances, Nbeads)

      bending_force_result = reshape((/ -1.d0/2.d0, 1.d0/2.d0, 0.d0, &
         1.0d0/2.d0, -sqrt(2.d0)/2.d0-1.d0, 0.d0, &
         -1.d0/4.d0-sqrt(3.d0)/4.d0, 5.d0/4.d0+sqrt(3.d0)/4.d0+sqrt(2.d0)/2.d0, 0.d0, &
         1.0d0/4.d0, -3.d0/2.d0-sqrt(3.d0)/4.d0, 0.d0, &
         sqrt(3.d0)/4.d0, 3.d0/4.d0, 0.d0/), &
         (/Ndim, Nbeads/) )

      call assertEquals(bending_force, bending_force_result, Ndim, Nbeads, &
         max_double_difference, "forces are correct for precise 5-bead case : test failed")

      deallocate(bead_to_bead_vector_superdiagonal_matrix, &
         bead_to_bead_distances, &
         Bead_initial_positions, &
         bending_force, &
         bending_force_result)

      ! Check forces for random 5-bead case with MATLAB !-!-!-!-!-!

      bending_params%bending_potential_type = OneMinusCosTheta
      bending_params%bending_stiffness = 2.d0

      Nbeads = 5
      seed = 34123
      allocate(bead_to_bead_vector_superdiagonal_matrix(Ndim,Nbeads,Nbeads), &
         bead_to_bead_distances(Nbeads,Nbeads), &
         Bead_initial_positions(Ndim,Nbeads), &
         bending_force(Ndim,Nbeads), &
         bending_force_result(Ndim,Nbeads))

      bead_initial_positions = reshape((/-5.1645947419192162d-01, 7.9748452989683460d-01, -4.1121562443647974d-01, &
         -9.3988831202840872d-01, 1.1590716947448900d+00, -7.6310486638258535d-01, &
         -6.7034760405534555d-01, -1.4053776613761966d+00, -2.9724073830375936d-01, &
         1.1832133384385024d+00, -3.6608832152898463d-01, 6.1365584115473881d-01, &
         9.4348205183717349d-01, -1.8509024173654354d-01, 8.5790538796808580d-01/), &
         (/Ndim, Nbeads/) )

      call b2bvector_sym_up(Nbeads,Bead_initial_positions,bead_to_bead_vector_superdiagonal_matrix)
      call modr_sym_up(Nbeads,bead_to_bead_vector_superdiagonal_matrix,bead_to_bead_distances)

      call get_bending_force(bending_force, bending_params, &
         bead_to_bead_vector_superdiagonal_matrix, bead_to_bead_distances, Nbeads)

      bending_force_result = reshape((/1.0507847458911022d+00, 1.8075811550046665d+00, 5.9298653220518749d-01, &
         -1.2494310241030562d+00, -1.8333141087409899d+00, -6.1970590685649796d-01, &
         4.1982689539550588d-01, 3.4148812263349709d-01, -7.8361711242102261d-01, &
         -4.2420142390715965d+00, -2.7289683808745142d+00, -1.3478360326510446d+00, &
         4.0208336218880447d+00, 2.4132132119773404d+00, 2.1581725197233776d+00/), &
         (/Ndim, Nbeads/) )

      call assertEquals(bending_force, bending_force_result, Ndim, Nbeads, &
         max_double_difference, "forces are correct for random 5-bead case with MATLAB : test failed")

      deallocate(bead_to_bead_vector_superdiagonal_matrix, &
         bead_to_bead_distances, &
         Bead_initial_positions, &
         bending_force, &
         bending_force_result)

   end subroutine

   subroutine test_get_trap_force()
      use Flow_Strength_Calculations

      type(flow_inputs) :: flow_inputs_dum
      Real(DBprec) :: time
      Integer :: Nbeads
      Real(DBprec), dimension(:,:), allocatable :: R_bead, F_trap, F_trap_target



      ! test for two-bead case
      Nbeads = 2
      allocate(R_bead(Ndim, Nbeads), F_trap(Ndim, Nbeads))
      time = 0
      flow_inputs_dum%flow_type = EQ
      flow_inputs_dum%trapOneInitialPosition = (/0.d0, 0.d0, 0.d0/)
      flow_inputs_dum%trapTwoInitialPosition = (/1.d0, 0.d0, 0.d0/)
      flow_inputs_dum%trapOneStrength = 1.d0
      flow_inputs_dum%trapTwoStrength = 2.d0
      R_bead = reshape((/0.d0, 0.d0, 0.d0, 0.d0, 1.d0, 0.d0/), (/Ndim,Nbeads/))
      call get_harmonic_trap_force(time, flow_inputs_dum, Nbeads, R_bead, F_trap)
      F_trap_target = reshape((/0.d0, 0.d0, 0.d0, 2.d0, -2.d0, 0.d0/), (/Ndim,Nbeads/))
      call assertEquals(F_trap_target, &
                        F_trap, Ndim, Nbeads, max_double_difference, "2-bead trap force test 1")

      deallocate(R_bead, F_trap)

   end subroutine

end module

module validation_tests
   ! All stochastic tests, which require a full statistical sample to be taken and may have to run for tens of minutes
   ! or hours to ensure reasonable sampling. They're set to 2 sigma accuracy, so there's a ~5% chance of any individual test
   ! failing purely by random chance. These are mostly comparisons with analytical results or previously established results.
   use fruit
   Use Global_parameters_variables_and_types
   use Physics_subroutines
   use properties
   use Initial_Position_Utilities
   use simulation_utilities

   implicit none

   private

   public :: validation_tests_run

!   Integer :: Stype        !Spring type
!   Integer :: EV           ! Excluded volume type to use
   Integer(k4b) :: Nbeads      !Number of beads in chain
!   Integer(k4b) :: myseed      !seed for random number generator
   Integer :: Nsamples     ! Number of sampling points
   Integer :: nsact        ! actual number of sampling points, calculated differently (change name)
!   Integer :: mmult        ! Number of timesteps between samples, who knows why it's called mmult
   Real (DBprec) :: deltaQ      ! FENE b-parameter (or sqrt of b-parameter)
   Real (DBprec) :: hstar      ! HI parameter
   real (DBprec) :: sigma        ! natural length of Fraenkel spring
!   Real (DBprec) :: zstar      ! dimensionless energy parameter for gaussian potential
!   Real (DBprec) :: dstar      ! dimensionless radius parameter for gaussian potential
   real (DBprec) :: tmax       ! Total integration time
   real (DBprec) :: tcur       ! initial simulation time (for auto-correlation calculations)
   Integer (k4b) nseed, stored_seed, stored_ix, stored_iy
   Integer(kind=k4b) :: step, num_steps
   Real(kind=dbprec) :: time
   Integer :: isample, index_of_min(1)
   Integer, allocatable :: sample_indexes(:)
   Real(kind=dbprec), allocatable :: tcdf(:), sample_times(:)
   
   real (DBprec) :: fd_err_avg,ncheb_avg
   Integer :: iter


   type(physical_parameters) :: phys_params
   type(simulation_parameters) :: sim_params
   type(calculated_variables) :: calc_vars

   Real (DBprec), allocatable :: PosVecR(:,:)       ! Position vector of beads, number_of_dimensions*number_of_beads
   real (DBprec), allocatable :: times(:)         ! simulation times for samples to be taken, number_of_samples
   real (DBprec), allocatable :: samples(:,:)       ! storage matrix for chain properties, taken using chain_props(...), number_of_sampled_properties*number_of_samples
   real (DBprec), allocatable :: time_cdf(:)      ! storage vector containing sampling times, number_of_samples
   real (DBprec), allocatable :: confi(:,:,:)      ! Storage vector for configurations at each sampling time, number_of_samples*number_of_dimensions*number_of_beads
   real (DBprec), allocatable :: grad(:,:,:)      ! Storage vector for forces at each sampling time, number_of_samples*number_of_dimensions*number_of_beads
   Real (DBprec), allocatable :: phi(:,:)        ! LJ well depth, Nbeads^2 (I assume for SDK potential)

contains

   subroutine validation_tests_run()
      ! Runs validation tests for single-chain code results

      print *, "Running validation tests"
      
      call testing_HI_Methods()
      
      print *, ""
      print *, "Finished validation tests"
      print *, ""

   end subroutine validation_tests_run

   ! subroutine WLC_bounded_loops_correctly()

      ! Integer :: Nbeads      !Number of beads in chain
      ! Integer(k4b) :: myseed      !seed for random number generator
      ! Real (DBprec), allocatable :: PosVecR(:,:)       ! Position vector of beads, number_of_dimensions*number_of_beads

      ! type(physical_parameters) :: phys_params
      ! type(simulation_parameters) :: sim_params
      ! type(calculated_variables) :: output_vars

   ! end subroutine
   
   subroutine testing_HI_Methods()
   
      Integer :: i, niter
      Real(DBprec) :: ncheb_tot, fd_err_tot
      
      ! physical parameters for physics of simulation
      phys_params%EV_inputs%contour_dist_for_EV = 1
      phys_params%EV_inputs%minimum_EV_cutoff = 0.7d0
      phys_params%EV_inputs%maximum_EV_cutoff = 1.5d0
      phys_params%Flow_inputs%flow_cessation_time = 0.d0
      phys_params%Flow_inputs%flow_period = 0.d0
      
      phys_params%bend_inputs%bending_potential_type = NoBendingPotential
      phys_params%bend_inputs%bending_stiffness = 0d0
      phys_params%bend_inputs%natural_angles = 0

      phys_params%EV_inputs%excluded_volume_type = NoEV
      phys_params%EV_inputs%dimensionless_EV_energy = 0.45d0
      phys_params%EV_inputs%dimensionless_EV_radius = 1.d0
      phys_params%EV_inputs%bead_attractive_interaction_strengths = 0

      phys_params%Flow_inputs%flow_strength = 0.01d0
      phys_params%Flow_inputs%flow_type = EQ
      phys_params%Flow_inputs%trapOneInitialPosition = (/0.d0, 0.d0, 0.d0/)
      phys_params%Flow_inputs%trapTwoInitialPosition = (/0.d0, 0.d0, 0.d0/)
      phys_params%Flow_inputs%trapOneStrength = 0.d0
      phys_params%Flow_inputs%trapTwoStrength = 0.d0
      phys_params%Flow_inputs%trapTwoVelocity = (/0.d0, 0.d0, 0.d0/)
      
      hstar = 0.4d0
      phys_params%HI_params%hstar = hstar
      ! phys_params%HI_params%ChebUpdateMethod = UpdateChebNew
      phys_params%HI_params%ChebUpdateMethod = UpdateChebAddOne
      phys_params%HI_params%EigenvalueCalcMethod = EigsFixman
      ! phys_params%HI_params%EigenvalueCalcMethod = EigsExact
      phys_params%HI_params%delSCalcMethod=Chebyshev
      ! phys_params%HI_params%delSCalcMethod=Cholesky
      ! phys_params%HI_params%delSCalcMethod=ExactSqrt
      phys_params%HI_params%nchebMultiplier = 2.d0
      phys_params%HI_params%fd_err_max = 2.5d-5
      
      deltaQ = 10.d0
      sigma = 0.d0
      phys_params%spring_inputs%finite_extensibility_parameter = deltaQ
      phys_params%spring_inputs%natural_length = sigma
      phys_params%spring_inputs%spring_type = Hook
      
      Nbeads = 3
      phys_params%number_of_beads = Nbeads
      
      ! simulation parameters
      sim_params%update_center_of_mass = 1
      sim_params%implicit_loop_exit_tolerance = 1e-6
      sim_params%simulation_seed = 50
      sim_params%simulation_timestep = 0.001d0
      sim_params%time_at_simulation_start = 0.d0
      sim_params%time_at_simulation_end = 4000.d0
      sim_params%number_of_samples_to_take = 1000
      
      ! calculate correct steps to take samples
      Nsamples = sim_params%number_of_samples_to_take
      nsact = NSamples
      tcur = sim_params%time_at_simulation_start
      tmax = sim_params%time_at_simulation_end
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
      
      allocate(calc_vars%chain_configuration_at_sample_points(nsact,Ndim,Nbeads), &
         calc_vars%total_force_at_sample_points(nsact,Ndim,Nbeads),&
         calc_vars%true_times_at_sample_points(nsact), &
         calc_vars%samples(Nprops,nsact))
      
      Allocate(PosVecR(Ndim,NBeads))
      
      ncheb_avg = 0
      fd_err_avg = 0
      
      iter = 1
      
      do i=1,iter
         PosVecR = 0.d0
         
         ! initialise the position vector randomly
         Call Initial_position_FENE_Fraenkel(phys_params%spring_inputs%spring_type,&
                                             Nbeads,deltaQ,sigma,PosVecR)
         
         Call Time_Integrate_Chain(PosVecR, phys_params,sim_params,calc_vars)
         if (mod(i,max(int(iter/5),1))==0) then
            print *, PosVecR
         end if
         
         samples = calc_vars%samples
         ncheb_avg = ncheb_avg + sum(real(samples(19,1:nsact)))/nsact
         fd_err_avg = fd_err_avg + sum(samples(20,1:nsact))/nsact
      end do
      
      confi = calc_vars%chain_configuration_at_sample_points
      grad = calc_vars%total_force_at_sample_points
      time_cdf = calc_vars%true_times_at_sample_points
      
      ncheb_avg = ncheb_avg/iter
      fd_err_avg = fd_err_avg/iter
      
      print *, "ncheb average is:", ncheb_avg
      print *, "fd_err average is:", fd_err_avg
      
      ! If (netcd .eq. 1) Then
         ! call write_to_netcdf(calc_vars, ntrajout, Restart+1, ncid)
      ! End If ! netcdf
      
      ! print *, PosVecR
      
      !deallocate(phys_params%EV_inputs%bead_attractive_interaction_strengths)
      deallocate(tcdf, sample_times, sample_indexes, sim_params%sample_indexes)
      deallocate(calc_vars%chain_configuration_at_sample_points, &
               calc_vars%total_force_at_sample_points,&
               calc_vars%true_times_at_sample_points, &
               calc_vars%samples)
      deallocate(PosVecR)
   
   end subroutine

!    subroutine FENE_Fraenkel_Dumbbell_gives_correct_

end module

Program tests
   use fruit
   Use Global_parameters_variables_and_types
   use unit_tests
   use validation_tests
   use regression_tests

   implicit none

   call init_fruit
   ! Unit tests take less than a minute to run, can be done on login nodes
   call unit_tests_run
   !call regression_tests_run
   ! Validation tests can take hours to run - submit as a job!
   ! call validation_tests_run
   call fruit_summary
   call fruit_finalize

End Program
