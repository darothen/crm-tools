! crm_io.F90 --
!
! Collection of subroutines useful for reading and writing output and input files
! for the cloud resolving model (Wang and Chang, 1993). Although not explicited wrapped
! here as a module, these subroutines are intended to be compiled with the Numpy 'f2py' 
! utility and accessed through Python. 
!
! To compile, execute the following command:
!   bash$ f2py -c -m crm_io crm_io.F90

subroutine read_2aq(filename, nt, nx, ny, nz, aq_c, aq_r)

    implicit none

!-- Input Variables
    character(*), intent(in) :: filename
    integer, intent(in)      :: nt, nx, ny, nz
!-- Output Variables
    real(4), dimension(nt,nx,ny,nz), intent(out) :: aq_c, aq_r

!-- Local Variables
    real(4), dimension(nx,ny,nz) :: c_step, r_step
    integer :: t
    integer, dimension(3) :: shape

!-- F2PY VARIABLE BINDINGS
    !  f2py intent(in) :: filename, nt, nx, ny, nz
    !  f2py intent(out) :: aq_c, aq_r
    !  f2py intent(hide) :: shape, c_step, fr_step, t

    shape = (/ nx, ny, nz /)
    open (unit=1, file=filename, form="unformatted", access="sequential")
    do t = 1, nt
        read (1) c_step
        aq_c(t,:,:,:) = reshape(c_step, shape)
        read (1) r_step
        aq_r(t,:,:,:) = reshape(r_step, shape)
    end do
    close (unit=1)

end subroutine

subroutine read_diag(filename, nz, spmd, counter)

    implicit none

!-- Input Variables
    character(*), intent(in) :: filename
    integer, intent(in) :: nz
    logical, intent(in) :: spmd

!-- Output Variables
    integer, intent(out) :: counter

!-- Local Variables
    integer :: status, nn, ntime, nhour, nmin, nsec
    real(4), dimension(nz, 46) :: tdiag

!-- F2PY VARIABLE BINDINGS
    ! f2py intent(in) :: filename, nz, spmd
    ! f2py intent(out) :: counter

!-- Main routine
    open (unit=90, file=filename, form="unformatted", access="sequential")
101 format (2x,i2,":",i0.2,":",i0.2,"     |     (",i0.8," seconds)")

    counter = 0
    if (spmd) then
        read_file : do 
            read (90,iostat=status) ntime, nhour, nmin, nsec, (tdiag(:,nn), nn=1,46)
            if (status > 0) then
                print *, "Error reading " // filename 
                exit
            else if (status < 0) then
                print *, "Finished reading " // filename
                exit
            else ! successfully read from file
                write (*,101) nhour, nmin, nsec, ntime
                counter = counter + 1
                write (*,*) counter, tdiag(16,1)
            end if
        end do read_file
    else
        print *, "Not implemented yet"
    endif

    !rewind (90)
    print *, "End of reading program"
    close (90)

end subroutine read_diag

subroutine save_diag(filename, nz, nt, spmd, all_time, all_tdiag)

    implicit none

!-- Input Variables
    character(*), intent(in) :: filename
    integer, intent(in) :: nz
    integer, intent(in) :: nt
    logical, intent(in) :: spmd

!-- Output Variables
    integer, dimension(nt), intent(out) :: all_time
    real(8), dimension(nt,nz,46), intent(out) :: all_tdiag

!-- Local Variables
    integer :: status, nn, ntime, nhour, nmin, nsec, t
    real(8), dimension(nz, 46) :: tdiag

!-- F2PY VARIABLE BINDINGS
    ! f2py intent(in) :: filename, nz, nt, spmd
    ! f2py intent(hide) :: tdiag, status, nn, ntime, nhour, nmin, nsec, t
    ! f2py intent(out) :: all_tdiag, all_time

!-- Main routine
    open (unit=90, file=filename, form="unformatted", access="sequential")

    if (spmd) then
        read_file : do t = 1, nt
            read (90,iostat=status) ntime, nhour, nmin, nsec, (tdiag(:,nn), nn=1,46)
            if (status > 0) then
                print *, "Error reading " // filename 
                exit
            else if (status < 0) then
                print *, "Finished reading " // filename
                exit
            else ! successfully read from file
                print *, t, nhour, nmin, nsec, ntime, tdiag(10,36)
                all_time(t) = ntime
                all_tdiag(t,:,:) = tdiag(:,:)
            end if
        end do read_file
    else
        print *, "Not implemented yet"
    endif

    print *, "End of reading program"
    close (90)

end subroutine save_diag


subroutine read(filename, nt, nx, ny, nz, data)

    implicit none
    
!-- Input Variables
    integer, intent(in)   :: nt, nx, ny, nz
    character(*), intent(in) :: filename
    
!-- Output Variables
    real(4), dimension(nt,nx,ny,nz), intent(out) :: data

!-- Local Variables
    real(4), dimension(nx*ny*nz) :: data_record
    integer, dimension(3) :: shape
    integer :: t, status
    
!-- F2PY VARIABLE BINDINGS
    !  f2py intent(in) :: filename, nt, nx, ny, nz
    !  f2py intent(out) :: data
    !  f2py intent(hide) :: data_record, shape, t

    shape = (/ nx, ny, nz /)
    open (unit=1, file=filename, form="unformatted", access="sequential")
    
    read_loop : do t = 1, nt
        !print *, t
        read (1,iostat=status) data_record
        if (status > 0) then
            print *, "Error reading " // filename
            exit
        else if (status < 0) then
            print *, "Finished reading " // filename
            exit
        else ! the read was successful, so do some work
            write (*,'(i3,"  MAX: ",e10.3," | MIN: ",e10.3)') t, maxval(data_record), minval(data_record)
            data(t,:,:,:) = reshape(data_record, shape)
        end if
    end do read_loop
    print *, "closing"
    close (unit=1)
    print *, "done"

end subroutine read

subroutine read_basic_state(filename, nz, p0, pt0, T0, qv0)

    implicit none

!-- Input Variables
    integer, intent(in)   :: nz
    character(*), intent(in) :: filename
    
!-- Output Variables
    real (kind=8), dimension(nz), intent(out) :: p0, pt0, T0, qv0

    integer :: k

!-- F2PY VARIABLE BINDINGS
    !  f2py intent(in) :: filename, nz
    !  f2py intent(out) :: p0, pt0, T0, qv0

    open (unit=1, file=filename, form="unformatted", access="sequential")
    read (1) p0, pt0, T0, qv0
    close (unit=1)

    do k = 1, nz
        print *, p0(k), pt0(k), T0(k), qv0(k)
    end do


end subroutine read_basic_state

subroutine write_initial_profile(nk, filename, p0, height0, t0, u00, v00, qv0)

    implicit none

!-- Input variables
    character(*), intent(in) :: filename
    real, dimension(nk), intent(in) :: p0, height0, t0, u00, v00, qv0

!-- Output variables
!   - none

!-- Local variables
    integer ::  iuseless, nk, k
    real :: useless
    real(4) :: p, height, t, u, v, qv

!-- F2PY VARIABLE BINDINGS
    !  f2py intent(in) :: filename, p0, height0, t0, u00, v00, qv0
    !  f2py intent(hide) :: k, p, height, t, u, v, qv, useless, iuseless, nk
    
    open(unit=10, file=filename, form="unformatted", access="sequential")

    useless = 1.0
    iuseless = 1

    do k = 1, nk
       print *, k
       p = p0(k)
       height = height0(k)
       t = t0(k)
       u = u00(k)
       v = v00(k)
       qv = qv0(k)
       iuseless = k
       write(10) p, height, t, useless, u, v, qv, iuseless
    end do

    close(10)

end subroutine

subroutine read_initial_profile(filename, nk, annika, profiles)

    implicit none 

!-- Input variables
    character(*), intent(in) :: filename
    integer, intent(in) :: nk
    logical, intent(in) :: annika

!-- Output variables
    real, dimension(nk, 8), intent(out) :: profiles

!-- Local variables
    real(4) :: p, h, t, u, v, qv
    real :: useless
    real, dimension(nk) :: var_profile
    integer :: iuseless, k, nvar

!-- F2PY VARIABLE BINDINGS
    !  f2py intent(in) :: filename, nk, annika
    !  f2py intent(hide) :: p, h, t, u, v, qv, useless, iuseless, k, nvar, var_profile
    !  f2py intent(out) :: profiles

    open(unit=10, file=filename, form="unformatted", access="sequential")

    if (annika) then
        print *, "reading annika format"
        nvar = 5
        do k=1,nvar
           read(10) var_profile
           print *, var_profile

           profiles(:, k) = var_profile
        end do
    else
        print *, "reading standard format"
        print *, "P   H   T   useless   u   v   qv  iuseless"
        do k=1,nk
           read(10) p, h, t, useless, u, v, qv, iuseless
           print *, p, h, t, useless, u, v, qv, iuseless

           profiles(k, :) = (/ p, h, t, real(useless), u, v, qv, real(iuseless) /)
        end do
    endif 

    close(10)

end subroutine

subroutine read_initial_chem(filename, nk, gas)

    implicit none

    character(*), intent(in) :: filename
    integer, intent(in) :: nk
    real(4), dimension(nk, 17), intent(out) :: gas    

    real(4), dimension(nk) :: o3, co, zco2, ho, ho2, h2o2, xno, xno2, xno3
    real(4), dimension(nk) :: xn2o5, hno3, ch4, ch2o, ch3o2h, so2, h2so4, dms
   
!-- F2PY VARIABLE BINDINGS
    !  f2py intent(in) :: filename, nk
    !  f2py intent(out) :: gas
    !  f2py intent(hide) :: o3, co, zco2, ho, ho2, h2o2, xno, xno2, xno3, xn2o5, hno3, ch4, ch2o, ch3o2h, so2, h2so4, dms

    open(unit=10, file=filename, form="unformatted", access="sequential")
    read(10) o3, co, zco2, ho, ho2, h2o2, xno, xno2, xno3, xn2o5, hno3, ch4, ch2o, ch3o2h, so2, h2so4, dms
    close(10)

    gas(:, 1) = o3
    gas(:, 2) = co
    gas(:, 3) = zco2
    gas(:, 4) = ho
    gas(:, 5) = ho2
    gas(:, 6) = h2o2
    gas(:, 7) = xno
    gas(:, 8) = xno2
    gas(:, 9) = xno3
    gas(:,10) = xn2o5
    gas(:,11) = hno3
    gas(:,12) = ch4
    gas(:,13) = ch2o
    gas(:,14) = ch3o2h
    gas(:,15) = so2
    gas(:,16) = h2so4
    gas(:,17) = dms

    return

end subroutine

subroutine write_initial_chem(filename, nk, gas)

    implicit none

!-- Input variables
    character(*), intent(in) :: filename
    real(4), dimension(0:nk-1, 17) :: gas

!-- Local variables
    integer :: nk
    real(4), dimension(nk) :: o3, co, zco2, ho, ho2, h2o2, xno, xno2, xno3
    real(4), dimension(nk) :: xn2o5, hno3, ch4, ch2o, ch3o2h, so2, h2so4, dms

!-- F2PY VARIABLE BINDINGS
    !  f2py intent(in) :: filename, gas
    !  f2py intent(hide) :: nk, o3, co, zco2, ho, ho2, h2o2, xno, xno2, xno3, xn2o5, hno3, ch4, ch2o, ch3o2h, so2, h2so4, dms 

    o3     = gas(:, 1)
    co     = gas(:, 2)
    zco2   = gas(:, 3)
    ho     = gas(:, 4)
    ho2    = gas(:, 5)
    h2o2   = gas(:, 6)
    xno    = gas(:, 7)
    xno2   = gas(:, 8)
    xno3   = gas(:, 9)
    xn2o5  = gas(:,10)
    hno3   = gas(:,11)
    ch4    = gas(:,12)
    ch2o   = gas(:,13)
    ch3o2h = gas(:,14)
    so2    = gas(:,15)
    h2so4  = gas(:,16)
    dms    = gas(:,17)

    open(unit=10, file=filename, form="unformatted", access="sequential")
    write(10) o3, co, zco2, ho, ho2, h2o2, xno, xno2, xno3, xn2o5, hno3, ch4, ch2o, ch3o2h, so2, h2so4, dms
    close(10)

end subroutine

subroutine calc_perturb(data_arr, nt, nx, ny, nz, means, perts)

    implicit none
    
!-- Input Variables
    real(4), dimension(nt, nx, ny, nz), intent(in) :: data_arr
    
!-- Output Variables
    real(4), dimension(nt, nx, ny, nz), intent(out) :: means, perts

!-- Local Variables
    real(4), dimension(nz) :: vertical_mean
    integer :: nt, nx, ny, nz
    integer :: t, x, y
    
!-- F2PY VARIABLE BINDINGS
    !  f2py intent(in) :: data_arr
    !  f2py intent(out) :: means, perts
    !  f2py intent(hide) :: nt, nx, ny, nz, t, x, y, vertical_mean

    do t = 1, nt
        do x = 1, nx
            do y = 1, ny
                vertical_mean = sum(data_arr(t, x, y, :))/nz 
                means(t, x, y, :) = vertical_mean   
                perts(t, x, y, :) = data_arr(t, x, y, :)-vertical_mean
            end do
        end do
    end do

end subroutine
