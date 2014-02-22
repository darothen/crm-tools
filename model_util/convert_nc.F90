!
! Convert CRM71 binary output files into NetCDF
! $gfortran -o convert_nc -I/usr/local/include convert_nc.F90 -lnetcdff -lnetcdf

module nc_helper 

    use netcdf
    implicit none

    contains 

        subroutine check(status)
            integer, intent(in) :: status

            if (status /= nf90_noerr) then
                print *, trim(nf90_strerror(status))
                stop "Stopped"
            end if

        end subroutine check

end module nc_helper


program convert_nc

    use netcdf
    use nc_helper
    implicit none

!----------------------------------------------------------------------

    ! File name to write to
    character(len=*), parameter :: FILE_NAME = "test.nc"

    ! Set-up data dimension grid
    integer, parameter :: NDIMS = 2
    integer, parameter :: NX = 6, NY = 12

    ! NetCDF id containers - one for files, variabiles, dimensions
    integer :: ncid, varid, dimids(NDIMS)
    integer :: x_dimid, y_dimid

    ! Working data array to write to disk
    integer :: data_out(NY, NX)

    ! Loop indexes; error handling
    integer :: x, y

    ! Create some imgainary data (model output)
    do x = 1, NX
        do y = 1, NY
            data_out(y, x) = (x - 1)*NY + (y - 1)
        end do
    end do

    ! Always check the return code of every netCDF function call. In
    ! this example program, wrapping netCDF calls with "call check()"
    ! makes sure that any return which is not equal to nf90_noerr (0)
    ! will print a netCDF error message and exit.

    ! Create the netCDF file. The nf90_clobber parameter tells netCDF to
    ! overwrite this file, if it already exists.
    call check( nf90_create(FILE_NAME, NF90_CLOBBER, ncid) )

    ! Define the dimensions. NetCDF will hand back an ID for each one.
    call check( nf90_def_dim(ncid, "x", NX, x_dimid) )
    call check( nf90_def_dim(ncid, "y", NY, y_dimid) )

    ! The dimids array is used to pass the IDs of the dimensions of
    ! the variables. Note that in fortran arrays are stored in
    ! column-major format.
    dimids = (/ y_dimid, x_dimid /)

    ! Define the variable, using 4-byte integers
    call check( nf90_def_var(ncid, "data", NF90_INT, dimids, varid) )

    ! End the definitions mode; NetCDF now knows wer're done defining metadata
    call check( nf90_enddef(ncid) )

    ! Write the data to the file in one single operation
    call check( nf90_put_var(ncid, varid, data_out) )

    ! Close the file; this frees up any internal resources associated with
    ! it and flushes any buffers
    call check( nf90_close(ncid) )

    print *, "Finished writing " // FILE_NAME

end program convert_nc