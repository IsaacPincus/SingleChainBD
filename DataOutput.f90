! Utlilities for Data ouptut

Module Netcdf_Utilities
   use netcdf
   use Global_parameters_variables_and_types
   use properties
   ! Creates and writes configs to netcdf files

   implicit none
   private

   public create_netcdf_file, write_to_netcdf, close_netcdf_file, get_config_and_time_from_netcdf

   Integer :: sample_dimid, dimensions_dimid, beads_dimid, dimids(4), traj_dimid, time_dimids(2), &
      cofm_dimids(3), cofm_count(3), cofm_start(3), time_start(2), time_count(2), &
      config_start(4), config_count(4), &
      varid_time, varid_confi, varid_grad, varid_cofm

   ! Using one set of varid for both VR and non-VR may lead to hard-to-diagnose bugs down the road...
   ! Keep a lookout for this issue. At the moment it seems they map to each other if the netCDF files
   ! are originally created in the exact same way.

contains

   subroutine create_netcdf_file(filePath, Nsamples, Ndimensions, Nbeads, netcdfID)
      Character(len=*), intent(in) :: filePath
      Integer, intent(in) :: Nsamples, Ndimensions, Nbeads
      Integer, intent(out) :: netcdfID
      logical :: fileExists
      integer :: cofm_status

      ! check whether the netCDF file already exists
      inquire(file=filePath, exist=fileExists)
      if (fileExists) then
         call check(nf90_open(filePath, NF90_WRITE, netcdfID))
         ! get the dimension IDs
         call check( nf90_inq_dimid(netcdfID, "No_of_samples", sample_dimid) )
         call check( nf90_inq_dimid(netcdfID, "Ndim", dimensions_dimid) )
         call check(nf90_inq_dimid(netcdfID, "NBeads", beads_dimid))
         call check(nf90_inq_dimid(netcdfID, "Ntraj", traj_dimid))
         call check(nf90_inq_varid(netcdfID, "Time", varid_time))
         call check( nf90_inq_varid(netcdfID, "configuration", varid_confi))
         call check(nf90_inq_varid(netcdfID, "Gradient", varid_grad))
         cofm_status = nf90_inq_varid(netcdfID, "cofm", varid_cofm)
         if (cofm_status == nf90_enotvar) then
            cofm_dimids = (/traj_dimid, sample_dimid, dimensions_dimid/)
            call check(nf90_def_var(netcdfID, "cofm", NF90_DOUBLE, cofm_dimids,varid_cofm))
         elseif (cofm_status /= nf90_noerr) then
            print *, trim(nf90_strerror(cofm_status))
            stop 2
         end if
      else
         Call check(nf90_create(filePath, NF90_NETCDF4, netcdfID))
         !   Define the dimensions. NetCDF will hand back an ID for each.
         call check( nf90_def_dim(netcdfID, "No_of_samples", NF90_UNLIMITED, sample_dimid) )
         call check( nf90_def_dim(netcdfID, "Ndim", Ndimensions, dimensions_dimid) )
         call check(nf90_def_dim(netcdfID, "NBeads", Nbeads, beads_dimid))
         call check(nf90_def_dim(netcdfID, "Ntraj", NF90_UNLIMITED, traj_dimid))
         dimids =  (/traj_dimid, sample_dimid, dimensions_dimid, beads_dimid /)
         time_dimids = (/traj_dimid, sample_dimid/)
         cofm_dimids = (/traj_dimid, sample_dimid, dimensions_dimid/)
         call check(nf90_def_var(netcdfID, "Time", NF90_DOUBLE,time_dimids, varid_time, chunksizes=(/32,32/)))
         call check( nf90_def_var(netcdfID, "configuration", NF90_DOUBLE, dimids,varid_confi, chunksizes=(/32,32,Ndim,Nbeads/)) )
         call check(nf90_def_var(netcdfID, "Gradient", NF90_DOUBLE, dimids,varid_grad, chunksizes=(/32,32,Ndim,Nbeads/)))
         call check(nf90_def_var(netcdfID, "cofm", NF90_DOUBLE, cofm_dimids,varid_cofm, chunksizes=(/32,32,Ndim/)))
         call check( nf90_enddef(netcdfID) )
      end if

      time_count = (/1, Nsamples/)
      config_count = (/1, Nsamples, Ndimensions, Nbeads/)
      cofm_count = (/1, Nsamples, Ndimensions/)

   end subroutine

   subroutine write_to_netcdf(output_variables, trajectoryIndex, sampleIndex, netcdfID)
      !Real(kind=DBprec) :: times(:), configurations(:,:,:), gradients(:,:,:), cofms(:,:)
      type(calculated_variables), intent(in) :: output_variables
      Integer, intent(in) ::  trajectoryIndex, netcdfID, sampleIndex

      time_start = (/trajectoryIndex, sampleIndex/)
      config_start = (/trajectoryIndex, sampleIndex, 1, 1/)
      cofm_start = (/trajectoryIndex, sampleIndex, 1/)

      call check(nf90_put_var(netcdfID, varid=varid_time, &
                 values=output_variables%true_times_at_sample_points, start=time_start, count=time_count))
      call check( nf90_put_var(netcdfID, varid=varid_confi, &
                 values=output_variables%chain_configuration_at_sample_points, start=config_start, count=config_count) )
      call check(nf90_put_var(netcdfID, varid=varid_grad, &
                 values=output_variables%total_force_at_sample_points, start=config_start, count=config_count))
      call check(nf90_put_var(netcdfID, varid=varid_cofm, &
                 values=output_variables%center_of_mass_cumulative, start=cofm_start, count=cofm_count))
      call check(nf90_sync(netcdfID))

   end subroutine

   subroutine get_config_and_time_from_netcdf(filePath, config, time, cofm_initial, trajectoryIndex, &
                                              sampleIndex, NDimensions, NBeads)
      Character(len=*), intent(in) :: filePath
      integer, intent(in) :: trajectoryIndex, sampleIndex, NDimensions, NBeads
      Real(kind=DBprec), intent(out) :: time, config(:,:), cofm_initial(:)
      Real(kind=DBprec) :: times(1,1), configs(1,1,NDimensions,NBeads), cofms(1,1,Ndim)
      logical :: fileExists
      Integer :: netcdfID
      integer :: cofm_status

      inquire(file=filePath, exist=fileExists)
      if (fileExists) then
         call check(nf90_open(filePath, NF90_NOWRITE, netcdfID))
         ! get the dimension IDs
         call check( nf90_inq_dimid(netcdfID, "No_of_samples", sample_dimid) )
         call check( nf90_inq_dimid(netcdfID, "Ndim", dimensions_dimid) )
         call check(nf90_inq_dimid(netcdfID, "NBeads", beads_dimid))
         call check(nf90_inq_dimid(netcdfID, "Ntraj", traj_dimid))
         call check(nf90_inq_varid(netcdfID, "Time", varid_time))
         call check( nf90_inq_varid(netcdfID, "configuration", varid_confi))
         call check(nf90_inq_varid(netcdfID, "Gradient", varid_grad))
         cofm_status = nf90_inq_varid(netcdfID, "cofm", varid_cofm)

         call check(nf90_get_var(netcdfID, varid=varid_time, values=times, &
                     start=(/trajectoryIndex, sampleIndex/), count=(/1,1/)))
         call check(nf90_get_var(netcdfID, varid=varid_confi, values=configs, &
                     start=(/trajectoryIndex, sampleIndex, 1, 1/), count=(/1,1,NDimensions,NBeads/)))
         if (cofm_status == nf90_noerr) then
            call check(nf90_get_var(netcdfID, varid=varid_cofm, values=cofms, &
                        start=(/trajectoryIndex, sampleIndex, 1/), count=(/1,1,NDimensions/)))
         else
            cofms(1,1,:) = 0.d0
         end if

         time = times(1,1)
         config = configs(1,1,:,:)
         cofm_initial = cofms(1,1,:)
      else
         print *, "NetCDF file " // trim(filePath) // " does not exist"
      end if

   end subroutine

   subroutine close_netcdf_file(netcdfID)
      Integer, intent(in) :: netcdfID

      call check(nf90_close(netcdfID))
   end subroutine

   subroutine check(status)
      integer, intent ( in) :: status

      if(status /= nf90_noerr) then
         print *, trim(nf90_strerror(status))
         stop 2
      end if
   end subroutine check

End Module
