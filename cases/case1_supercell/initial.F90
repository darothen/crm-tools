
#include "ctrparam.h"

! ================================================================
!
!  INITIAL.F90                   
!
!  Purpose:
!	A package containing initialization subroutines.			  
!
!  Author
!	Chien Wang
!	MIT Joint Program on Science and Policy of Global Change
!
!  Revision:
!	Date	By		Brief Description
!	----	--		-----------------
!	021898	Chien Wang	F90/F95 version	
!	070698	Chien Wang	cpp
!	070698	Chien Wang	revised
!	110698	Chien Wang	aqueous chemistry
!	061099	Chien Wang	revised printouts
!	020100	Chien Wang	add avden0
!	062200	Chien Wang	spmd
!	071900	Chien Wang	deleted water and wind and
!				  renamed initial to prt_conf	
!				  renamed status to getstatus
!	072000	Chien Wang	deleted dummies in several subroutines
!	091400	Chien Wang	use hydrometeor
!	092100	Chien Wang	rearrange readdata
!	101800	Chien Wang	use shared data modules
!	121900	Chien Wang	move init of 2/3d arrays to setrun.F
!	011101	Chien Wang	distributed memory init
!	011801	Chien Wang	state, fk, aerosol, & rad distmem
!	013101	Chien Wang	wind distmem
!	020901	Chien Wang	new bc
!	032101	Chien Wang	add multiple bubble option
!	032201	Chien Wang	add Coriolis
!	040301	Chien Wang	write info for multi-bulb cases
!	042401	Chien Wang	a small fix on printout lines
!	042701	Chien Wang	add u0_sin_sin & u0_cos_cos
!	041102	Chien Wang	do loop in f90/f95 format
!	041902	Chien Wang	add heterogeneous chemistry
!	103102	Chien Wang	add print out for aerosol scheme
!	021903	Chien Wang	rev. for INPUT_32B
!	021904	Chien Wang	add XYBULB
!	032304	Chien Wang	add time mean data description
!	060104	Chien Wang	add initial cn for AERO_ENABLE
!	061604	Chien Wang	change float() to real()
!	062204	Chien Wang	calculating fcoriolis
!	120304	Chien Wang	change xnc0_c to nsait_0 for init aitken
!	060705	Chien Wang	add lightning description
!	102906	Chien Wang	add pt0_g & qv0_g
!	081607	Chien Wang	f90 free format & changes in write()
!	110907	Chien Wang	TKE
!	060210	Chien Wang	add if_firststep in init_disturb
!
! ================================================================
! $Id: initial.F90 24 2013-10-07 19:59:00Z wangc $
! ================================================================

	subroutine prt_conf

!    =================================================
!                                                             
!    Subroutine for launching initial perturbation
!                 of initial values of u,v,w,ptil,& qv
!
!		F90/F95 version
!               ---------------------------------    
!    Author:   Chien Wang                           
!              MIT Joint Program on Sciece & Policy of
!                  Global Change                     
!                                                   
!    Last Revision:	July 6, 1998
!			July 19, 2000 renamed from
!				initial to prt_conf 	   
!                                                 
!    =================================================

	USE gridno
	USE shared_data
	
	character (len = 12) :: yes_no
	character (len = 12) :: yes_no2
	character (len = 4 ) :: cwarm
	
	if(nwarm.eq.1)then
	  cwarm='warm'
	else
	  cwarm='cold'
	endif
         
	write(6,'(1x/6x,"CM4.0 ",i2,"-dimensional ",a4," version  Case: ",a20,/)')ndims,cwarm,casename
	write(6,'(8x,"sin(d)sin(lat) & cos(d)cos(lat) = ",2f9.5/)')u0_sin_sin,u0_cos_cos

#ifdef MEANDATA
          yes_no = 'Time Means'
#else
          yes_no = 'Instant Data'
#endif
	write(6,'(8x,"Output Fields are Derived based on ",a12)')yes_no
	write(6,'(8x,48("=")/)')

#ifdef XPLBC_ENABLE
	  yes_no = 'X-Periodical'
#else
	  yes_no = 'X-Radiation'
#endif

#ifdef YPLBC_ENABLE
	  yes_no2 = 'Y-Periodical'
#else
	  yes_no2 = 'Y-Radiation'
#endif

	write(6,'(8x,"Basic time step = ",f8.2," second")'                    )dt
	write(6,'("        Small time step = ",f8.2," second")'               )dtau
	write(6,'("        dx & dy         = ",2f7.1," meter")'               )dx,dy
	write(6,'("        Gridpoint number= ",3i7," in x, y, & z direction")')nx,ny,nz
	write(6,'("        Lateral boundary: ",a12," type condition")'	      )yes_no
	write(6,'("                          ",a12," type condition")'	      )yes_no2

#ifdef RLVBC_ENABLE
          yes_no = 'Rigid-lid'
#else
          yes_no = 'Natural'
#endif
	write(6,'(8x,"Vertical boundary: ",a12," type condition")')yes_no
	write(6,'(8x,"  Top damp layer is:",i3," z-grid thick")')kdamp
	if(if_lbcdamp.ne.0)then
	  write(6,'(8x,"Apply lbc-damp with tau_lbcdamp = ",i9," times of dtau")')if_lbcdamp
	else
	  write(6,'(8x,"No lbc-damp applied")')
	endif
	write(6,'(8x,"  lbc relexation coefficient = ",f6.2)')eps

#ifdef FILTER_ENABLE
        if(if_filter.eq.0)then
	  write(6,'(8x,"No filter applied")')
	else
	  xxx = real(if_filter)*dt/60.0
          write(6,'(8x,"Apply filter every",f7.2," minutes")')xxx
        endif
#else
	write(6,'(8x,"No filter applied")')
#endif

	write(6,*)

#if ( defined SPMD )
	write(6,'(8x,"SPMD run using: ",i6," processors")')nproc
#else
	write(6,'(8x,"Single CPU run")')
#endif

#if ( ! defined INIT_3D )
	write(6,'(8x,"Horizontally Homogeneous Initial Fields")')
	write(6,'("          Initial data file = ",a50)')file_init
#else
	write(6,'(8x,"3D Initial Fields")')
	write(6,'("          Initial data file = ",a50)')file_init3d
#endif
	write(6,'(8x,"  nsmooth(1/0 call smooth) start with ",i2)')nsmooth
	write(6,'(8x,"  u_shift = ",f6.2," v_shift = ",f6.2," m/s")')u0shift,v0shift

#if ( defined CORIOLIS )
	write(6,'(8x,"Coriolis coefficient (1,1) and (nx,ny) = ",2f12.5)')fcor(1,1),fcor(nx,ny)
#endif

	write(6,*)

#ifdef TKE
	yes_no   = 'TKE'
#else
	yes_no   = '1st Order'
#endif
	write(6,'(8x,"Subgrid-scale closure scheme: ",a12)')yes_no

	write(6,*)

#ifdef AERO_ENABLE
	write(6,'(8x,"Initial Aitken mode CN = ",f10.2," 1/cc")')nsait_0*1.e-6
	write(6,'(8x,"Initial Accuml-mode CN = ",f10.2," 1/cc; IN = ",f10.2," 1/L")')xn_ccn0*1.e-6,xn_in0*1.e-3
#else
	write(6,'(8x,"In calculating CN nucleation rate Nc = Cs^k")')
	write(6,'(8x,"  C = ",f10.2," 1/cc;  k = ",f10.2)')xnc0_c*1.e-6,xnc0_k

	write(6,'(8x,"Reference concentrations are ")')
	write(6,'(8x,"  CCN = ",f10.2," 1/cc; IN = ",f10.2," 1/L")')xn_ccn0*1.e-6,xn_in0*1.e-3
#endif

#ifdef BERRY
	write(6,*)"       Berry-type autoconversion"
#endif
#ifdef KESSLER
	write(6,*)"       Kessler-type autoconversion"
#endif

	write(6,'(8x)')

#ifdef RAD_ENABLE
	write(6,'(8x,"Radiation is calculated at every ",f7.1," minute")')iradx*dt/60.0
#else
	write(6,'(8x,"Run without Radiation ")')
#endif

	write(6,'(8x)')

#ifdef CHEM_ENABLE
	write(6,'(8x,"Chemistry initial profile (1D) data file = ",a50)')file_cheminit

#ifdef CHEM_REACT
	write(6,'(8x,"Chemical Reactions are calculated at every ",f7.1," minute")')i_creact*dt/60.0
	write(6,'("          Photochemical tabulated data file = ",a50)')file_rktbl
#else
	write(6,'(8x,"Transport-only Chemistry ")')
#endif
#else
	write(6,'(8x,"Run without Chemistry ")')
#endif

#ifdef AQCHEM_ENABLE
	write(6,'(8x,"Aqueous Reactions are calculated at every ",f7.1," minute")')i_creact*dt/60.0
#else
	write(6,'(8x,"Run without Aqueous Chemistry ")')
#endif

#ifdef SOLIDCHEM_ENABLE
	write(6,'(8x,"Heterogeneous Reactions are calculated at every ",f7.1," minute")')i_creact*dt/60.0
#else
	write(6,'(8x,"Run without Heterogeneous Chemistry ")')
#endif

#ifdef AERO_ENABLE
	write(6,'(8x,"Run with multi-mode and multi-moment aerosol scheme ")')
#else
	write(6,'(8x,"Run with single mode and moment aerosol scheme ")')
#endif

	write(6,'(8x)')

#ifdef LIGHTNING
	write(6,'(8x,"Run with parameterized lightning scheme ")' )
	write(6,'(10x,"NO prodrate of cg & ic:",2e12.3," molecules/flash")')PROD_CG,PROD_IC
#endif

	write(6,'(8x,60("=")/)')

	return
	 end


!	===============================================
	subroutine prt_profile (u00, v00, qv_p, ptil_p)
!	===============================================
!
! 	PRINT INITIAL VALUES:
!

	USE gridno
	USE shared_data

	IMPLICIT NONE
	
	integer :: mmx, k
	real	:: z, qq

#ifdef INPUT_32B
	real(4), dimension(nz) :: u00, v00
#else
	real, dimension(nz) :: u00, v00
#endif
	real, dimension(nz) :: qv_p, ptil_p
	
! -----------------------------------------------------------------------

	write(6,'("        Max & min dz    = ",2f7.1," meter"/)'               )maxval(deltaz),minval(deltaz)

	write(6,'(8x,60("=")/)')
	write(6,'(9x,"Z(M)",7x,"PAI0    T0(K)     PT0(K)  QV0(G/KG) U0(m/s)  V0(m/s)   P(mb)"/)')

	mmx = nx/2
	z   = sum(deltaz(2:nz))
	do k = nz,1,-1
          write(6,'(2x,i3,2x,f8.1,2x,f8.2,2x,f8.2,2x,f8.2,2x,f8.4,2x,f6.2,2x,f6.2,2x,f8.2)')	&
	  	k,z,p0(k),t0(k),ptil_p(k),qv_p(k),u00(k),v00(k),avp(k)        
          z  = z - deltaz(k)                         
	end do

        write(6,'(2x,a3,32x,f8.2,2x,f8.4)')"0",pt0_g(mmx,1),qv0_g(mmx,1)*1000.

	write(76)u00,v00

      return
       end

!	===============================
	Subroutine getstatus (u00, v00)
!	===============================

	USE gridno
	USE shared_data

	IMPLICIT NONE

	integer :: i, k, kmin, iuseless
	
	real	:: qv0min, useless
	real, dimension(nz) :: height
	
#ifdef INPUT_32B
	real(4) :: r_p, r_h, r_t, r_u, r_v, r_q
!        real(4), dimension (nz) :: r_p, r_h, r_t, r_u, r_v, r_q
	real(4), dimension (nz) :: u00,   v00
#else
	real, dimension (nz) :: u00,   v00
#endif

! ===============================================

  file_read: select case (casename)

	!
	! === ISDAC
	!
	case ("ISDAC")

	  do k = 1,nz
#if ( defined INPUT_32B )
	  	read(10) r_t,r_p,r_q,r_h,r_u,r_v
     	  	p0    (k) = r_p
	  	height(k) = r_h
	  	T0    (k) = r_t
	  	u00   (k) = r_u
	  	v00   (k) = r_v
	  	qv0   (k) = r_q
#else
	  	read(10) T0(k),p0(k),qv0(k),height(k),u00(k),v00(k)
#endif
		! --- unit = K, hPa, g/kg, m, m/s, m/s
	  end do

	  deltaz(1)    = height(1)
	  deltaz(2:nz) = height(2:nz) - height(1:nz-1)
	  deltaz(nwz)  = 0.5*deltaz(nz)
	  height(2:nz) = height(2:nz) - height(1)
	  height(1)    = 0
	  qv0  = qv0*1.e-3	! g/kg to kg/kg
  
	!
	! === CEPEX M8
	!
	case ("CEPEXM0308")
	  do k=1,nz

#ifdef INPUT_32B
        	read(10)		&
!        	read(10,106)		&
        	       r_p,	&
        	       r_h,	&
        	       r_t,	&
        	       r_u,	&
        	       r_u,	&
        	       r_v,	&
        	       r_q,	&
        	       iuseless
     		p0 (k)    = r_p
		height(k) = r_h
		T0 (k)    = r_t
		u00(k)    = r_u
		v00(k)    = r_v
		qv0(k)    = r_q
#else
        	read(10)			&
!        	read(10,106)			&
        	       p0(k),		&
        	       height(k),  	&
        	       T0(k),		&
        	       useless, 	&
        	       u00(k),  	&
        	       v00(k),  	&
        	       qv0(k),  	&
        	       iuseless
#endif
	  end do

!106	format(1x,f8.3,1x,f8.1,4(1x,f7.3),1x,e10.3,1x,i3)

	  do k=2,nz
		deltaz(k) = height(k) - height(k-1)
	  end do
		deltaz(1)   = 0.5*deltaz(2)
		deltaz(nwz) = 0.5*deltaz(nz)
  	
	!
	! === Others
	!
	case default
	  do k=1,nz

#ifdef INPUT_32B
        	read(10)		&
!        	read(10,106)		&
        	       r_p,	&
        	       r_h,	&
        	       r_t,	&
        	       r_u,	&
        	       r_u,	&
        	       r_v,	&
        	       r_q,	&
        	       iuseless
     		p0 (k)    = r_p
		height(k) = r_h
		T0 (k)    = r_t
		u00(k)    = r_u
		v00(k)    = r_v
		qv0(k)    = r_q
#else
        	read(10)			&
!        	read(10,106)			&
        	       p0(k),		&
        	       height(k),  	&
        	       T0(k),		&
        	       useless, 	&
        	       u00(k),  	&
        	       v00(k),  	&
        	       qv0(k),  	&
        	       iuseless
#endif
	  end do

	  do k=2,nz
		deltaz(k) = height(k) - height(k-1)
	  end do
		deltaz(1)   = 0.5*deltaz(2)
		deltaz(nwz) = 0.5*deltaz(nz)

  end select file_read

	kmin   = 1
	qv0min = 100.0
	
	avp = p0
	p0  = p0*100.0

!	T0(1) = T0(1) - 2.0
    !!! Hardcoded for supercell!
	PT0 = T0*(1.e5/p0)**(0.286)

    pt0_g(1:nx, 1:ny) = 314.094
    qv0_g(1:nx, 1:ny) = 13.0800/1000.

	!pt0_g(1:nx,1:ny) = T0(1)*(1.e5/p0(1))**(0.286)
	!qv0_g(1:nx,1:ny) = qv0(1)
    !!! - end

	do k=1,nz
	  if(qv0(k).le.qv0min)then
	  	kmin   = k
		qv0min = qv0(k)
	  endif
	end do
      
	if(kmin.lt.nz)then
	  do k=kmin,nz
		qv0(k) = max(0.0,qv0min)
	  end do
	endif

    !!! Modification for 2D squall case
	p0   = (p0*1.e-5)**(2.87e2/1.005e3)
	!ps   = 1.e5/2.87e2*p0**(1.005e3/2.87e2-1.)
    ps   = 978.0
	!!! - end
    pp0  = 1.e5*p0**(1.005e3/2.87e2)
        den0 = pp0/(2.87e2*pt0*p0)

	avden0(1:nzi) = (den0(1:nzi) + den0(2:nz))*0.5
	avden0(nz)    =  den0(nz)

        pxx (1:nz) = (1.e5/pp0(1:nz))**(.286)
        pxxi(1:nz) = (1.e5/pp0(1:nz))**0.3

	write(75)p0,pt0,T0,qv0

	return
	 end
                                                              
!	==================================================
	subroutine init_disturb (if_firststep, qvd, ptild)
!	==================================================

!	=============================================
!
!	Subroutine for launching initial disturbance
!	-------------------------------------------- 
!	Author:		Chien Wang 
!			MIT Joint Program for Science 
!			 and Policy of Global Change                       
!                                                         
!	Most Recent Revision:  February 23, 1998
!			       July 19, 2000	add dummy       
!                                                     
!	=============================================

	USE gridno
	USE shared_data

	IMPLICIT NONE

#if ( defined MULTIPLE_XBULB )||( defined MULTIPLE_YBULB )
	integer, parameter   :: b_start = BULB_START
	integer, parameter   :: b_int   = BULB_INT
#elif ( defined MULTIPLE_XYBULB )
	integer, parameter   :: b_start = BULB_START	! this is for x_start
	integer, parameter   :: b_starty= BULB_STARTY
	integer, parameter   :: b_int   = BULB_INT
#endif

	integer		     :: if_firststep	! 1/0 = yes/no
	integer		     :: m, i, j, k, nb, nby1, nby2, nbx1, nbx2
	real		     :: qvs
	real		     :: rx2, ry2, rz2
	real		     :: xx1, xx2, x, y, zz, xxx
	character (len = 12) :: yes_no

	real, dimension(nx,ny,nz) :: qvd, ptild
    real :: height_above_sfc, slope, incpt
		
! -----------------------------------------------------------------------

	if ( if_firststep .eq. 1 ) then ! only print this once
#if   ( defined MULTIPLE_XBULB )
	write(6,'(8x,"Convection was initiated by",i3," warm bulbs along ",a11)')nbulb,'x-direction'
	write(6,'(10x,"Start from grid",i5," with intervals of",i5," grids"/)')b_start+b_int,b_int
#elif ( defined MULTIPLE_YBULB ) 
	write(6,'(8x,"Convection was initiated by",i3," warm bulbs along ",a11)')nbulb,'y-direction'
	write(6,'(10x,"Start from grid",i5," with intervals of",i5," grids"/)')b_start+b_int,b_int
#elif ( defined MULTIPLE_XYBULB )
	write(6,'(8x,"Convection was initiated by",i3," warm bulbs along ")')nbulb
	write(6,'(10x,a11,"Start from grid",i5," with an interval of",i5," grids")')'x-direction',b_start,b_int
	write(6,'(10x,a11,"Start from grid",i5," with an interval of",i5," grids"/)')'y-direction',b_starty,b_int
#else
	write(6,'(8x,"Convection was initiated by one bulb(1) or random(0): ",i2/)')nbulb
#endif

      	write(6,'(8x,"position of bulb (center, left, right, if symmetric):")')
      	write(6,'("	    mx, mx1, mx2, ioe: ",4i6)')mx,mx1,mx2,ioe 
      	write(6,'("	    my, my1, my2, joe: ",4i6)')my,my1,my2,joe
      	write(6,'("	    mz, mz1, mz2, koe: ",4i6)')mz,mz1,mz2,koe 
      	write(6,'("                 delta Ptil: ",f6.2)')dptil 
      	write(6,'("                 delta Qv  : ",f6.2/)')dqv 

      	endif 

!----------------------------------------------
!  Initial perturbation for temperature and qv:
!

	rx2 = ((real(mx2 + mx1)*0.5)*dx)**2
#ifdef MODEL_3D
        ry2 = ((real(my2 + my1)*0.5)*dy)**2
#else
        ry2 = rx2
#endif

#if ( ! defined MULTIPLE_XBULB ) && ( ! defined MULTIPLE_XYBULB )
	mx1 = mx - mx1
	mx2 = mx + mx2
#endif

#if ( ! defined MULTIPLE_YBULB ) && ( ! defined MULTIPLE_XYBULB )
#ifdef MODEL_3D
        my1 = my - my1
        my2 = my + my2
#else
        my  = 1
        my1 = 1
        my2 = 1
#endif
#endif

	mz2 = mz + mz2
	mz1 = mz - mz1

	rz2 = ( 0.5*(sum(deltaz(1:mz2)) - sum(deltaz(1:mz1))) )**2   

  pert_method: select case (nbulb)

    !
    ! === Warm bulb method:
    !
    case (1)

#if ( defined MULTIPLE_XBULB )
	mx = b_start

	do 2012 nb = 1, nbulb

	mx   = mx + b_int
        nbx1 = mx - mx1
        nbx2 = mx + mx2

      do k=mz1,mz2
      do j=my1,my2
      do i=nbx1,nbx2

#elif ( defined MULTIPLE_YBULB ) && ( defined MODEL_3D )
	my = b_start

	do 2012 nb = 1, nbulb

	my   = my + b_int 
        nby1 = my - my1
        nby2 = my + my2

      do k=mz1,mz2
      do j=nby1,nby2
      do i=mx1,mx2
#elif ( defined MULTIPLE_XYBULB ) && ( defined MODEL_3D )
	!
	! === note the bulbs located starting from x:west to east
	!		and y:north to south
	!
	mx = b_start  - b_int
	my = b_starty + b_int	

	do 2012 nb = 1, nbulb

	mx   = mx + b_int
        nbx1 = mx - mx1
        nbx2 = mx + mx2

	my   = my - b_int 
        nby1 = my - my1
        nby2 = my + my2

      do k=mz1,mz2
      do j=nby1,nby2
      do i=nbx1,nbx2
#else
      do k=mz1,mz2
      do j=my1,my2
      do i=mx1,mx2
#endif
                                                            
	if(ioe.eq.1)then
          if(i.le.mx)then   
          	x = (abs(real(i-mx))+0.5)*dx  
          else                          
          	x = (real(i-mx-1)+0.5)*dx   
          endif                       
        else
          x = abs(real(i-mx))*dx
        endif

#ifdef MODEL_3D
	if(joe.eq.1)then
	  if(j.le.my)then          
		y = (abs(real(j-my))+0.5)*dy   
	  else                           
		y = (real(j-my-1)+0.5)*dy    
	  endif                        
	else
		y = abs(real(j-my))*dy
	endif
#else
		y = 0.0
#endif
                                 
        if(koe.eq.1)then
          if(k.le.mz)then
            zz = abs( sum(deltaz(1:k)) - sum(deltaz(1:mz)) ) + 0.5*deltaz(mz+1)
          else
            zz = sum(deltaz(1:k)) - sum(deltaz(1:mz)) - 0.5*deltaz(mz+1)
          endif
        else
          zz = abs( sum(deltaz(1:k)) - sum(deltaz(1:mz)) )
        endif

        x = x**2
        y = y**2
        zz= zz**2

	! === to make sure a round area:
	if ( (x + y + zz).le.(rx2 + ry2 + rz2) ) then
          X = sqrt( min(1.0, x/rx2 + y/ry2 + zz/rz2) )

          ptild(i,j,k) = ptild(i,j,k) + dptil*(cos(x*3.1415926*.5))**2


! 091696
	  xx1 = 1.e5*p0(k)**3.5017422
	  xx2 = pt0(k)*p0(k)
          qvs = 380./xx1*exp(17.269388*(1.0-237./(xx2-36.))) - qv0(k)

#ifdef MODEL_3D
	  qvd(i,j,k) = dqv*(1.-(X/sqrt(6.))**2)*qvs + qv0(k)
#else
	  qvd(i,j,k) = dqv*(1.-(X*.5)**2)*qvs + qv0(k)
	  !qvd(i,j,k) = dqv*(1.-(X/0.8)**2)*qvs + qv0(k)
	  !qvd(i,j,k) = dqv*qvs*(cos(x*3.1415926*.5))**2 + qv0(k)
#endif                         
	end if
            
      end do
      end do
      end do

2012	continue

    !
    ! === random perturbation:
    !
    case (0)
      
      m=783637
      do j=my1,my2
      do i=mx1,mx2

	!===========================================
	!   A random number generator from
	!     Phil Rasch     Feb. 1994
	!===========================================

        m   = 125*m
        m   = m - (m/2796203)*2796203
        xxx = real(m)/2796202.		!from 0-1

	xxx = max(0.0,2.0*(xxx - 0.5))	!only positive number allowed
	!===========================================
	
      	do k=mz1,mz2
! 091696
        	xx1 = 1.e5*p0(k)**3.5017422
        	xx2 = pt0(k)*p0(k)
        	qvs = 380./xx1*exp(17.269388*(1.0-237./(xx2-36.))) - qv0(k)

        	ptild(i,j,k) = ptild(i,j,k) + dptil*xxx
        	qvd  (i,j,k) = qv0(k)	    + dqv*qvs*xxx

	 end do

      write(7,*)xxx

      end do
      end do

	close (7)

    ! -- Cold Pool - 2D Squall Case
    case (2)
        write(6,'(8x,"Placing Cold Pool")')
        do i = 1, 200
            j = 1
            ! apply a  perturbation to the initial theta of -5 K at a height of 0 km, 
            ! linearly decreasing to 0 K at a height of 4.5 km
            do k = 1, nz
                height_above_sfc = sum(deltaz(1:k))
                slope = -5./4500.
                incpt = 5.
                if (height_above_sfc .lt. 4500.) then
                    ptild(i,j,k) = ptild(i,j,k) - (slope*height_above_sfc + incpt)
                    if (i .eq. 1) write(6,*) k, height_above_sfc, (slope*height_above_sfc + incpt)
                endif
            end do
        end do


  end select pert_method
  
      return
       end

!	===========================================================
	Subroutine readdata ( u2, u1, u0, v2, v1, v0, w2, w1, 		&    
        		      state1, state2, hydromtr1, hydromtr2,	&
        		      xn_ccn, xn_in,  gas, aqc, aqr, solidi )
!	===========================================================

!	=====================================================
!  	
!	Subroutine for initialization using saved output data	
!	-----------------------------------------------------	
!	Author:     Chien Wang 	
!		    MIT Joint Program for Science and Policy
!		        of Global Change
!	
!	Most Recent Revision:  November 6, 1998 	
! 	
!	=====================================================

	USE gridno
	USE typedef_state
	USE typedef_hydrometeor
	USE typedef_gas
	USE typedef_aq
	USE typedef_solid
	USE shared_data

	character (len = 12)	:: yes_no
	character (len = 12)	:: yes_no2
	real,dimension(nx,ny,nz):: den

  ! --- define wind components
#if ( defined SPMD )
	real,dimension (ipu_start:ipu_end,jp_start:jp_end,1:nz) :: u2, u1, u0
	real,dimension (ip_start:ip_end,jpv_start:jpv_end,1:nz) :: v2, v1, v0
	real,dimension (ip_start:ip_end,jp_start:jp_end, 1:nwz) :: w2, w1
#else
	real, dimension(nux,ny,nz) :: u2, u1, u0	! u wind        
	real, dimension(nx,nvy,nz) :: v2, v1, v0	! v wind
	real, dimension(nx,ny,nwz) :: w2, w1		! w wind
#endif

	type (atm_state), dimension(nx,ny,nz) :: state1, state2
	type (hydrometeor), dimension(nx,ny,nz) :: hydromtr1, hydromtr2
	real, dimension(nx,ny,nz) :: xn_ccn, xn_in
		
#ifdef CHEM_ENABLE
        type (gas_chemical), dimension(nx,ny,nz) :: gas
#else
        type (gas_chemical) :: gas
#endif

#ifdef AQCHEM_ENABLE
        type(aq_chemical), dimension(nx,ny,nz) :: aqc, aqr
#else
        type(aq_chemical) :: aqc, aqr
#endif

#ifdef SOLIDCHEM_ENABLE
        type(solid_chemical), dimension(nx,ny,nz) :: solidi
#else
        type(solid_chemical) :: solidi
#endif

!	--------------------------------------------------

	read(103)p0,pt0,T0,qv0
	read(104)u0,v0

	read(105)ntime,nhour,nmin,nsec,nwarm,dtp,dtaup
     
	read(105)u1
	read(105)u2
	read(105)v1
	read(105)v2
	read(105)w1
	read(105)w2
	read(105)ub
	read(105)vb
	read(105)wb
	read(105)state1%p
	read(105)state2%p
	read(105)state1%T
	read(105)state2%T
	read(105)state1%ptil
	read(105)state2%ptil
	read(105)state1%qv
	read(105)state2%qv
	read(105)hydromtr1
	read(105)hydromtr2
	read(105)xn_ccn
	read(105)xn_in
	read(105)dtnet
	read(105)fluxsu
	read(105)fluxsd
	read(105)gas%o3
	read(105)gas%co
	read(105)gas%ho
	read(105)gas%ho2
	read(105)gas%h2o2
	read(105)gas%xno
	read(105)gas%xno2
	read(105)gas%hno3
	read(105)gas%ch4
	read(105)gas%ch2o
	read(105)gas%ch3o2h
	read(105)gas%so2
	read(105)gas%h2so4
	read(105)gas%dms
	read(105)aqc%o3
	read(105)aqr%o3
	read(105)aqc%civ
	read(105)aqr%civ
	read(105)aqc%h2o2
	read(105)aqr%h2o2
	read(105)aqc%xnv
	read(105)aqr%xnv
	read(105)aqc%siv
	read(105)aqr%siv
	read(105)aqc%svi
	read(105)aqr%svi
	read(105)aqc%ch2o
	read(105)aqr%ch2o
	read(105)aqc%ch3o2h
	read(105)aqr%ch3o2h
	read(105)aqc%hplus
	read(105)solidi%o3
	read(105)solidi%h2o2
	read(105)solidi%xnv
	read(105)solidi%ch2o
	read(105)solidi%ch3o2h
	read(105)solidi%siv
	read(105)solidi%svi
     
	do k=1,nz
	do j=1,ny
	do i=1,nx
#ifdef CHEM_ENABLE
	  if (gas(i,j,k)%o3    .lt.0.0) gas(i,j,k)%o3	  = 0.0
	  if (gas(i,j,k)%co    .lt.0.0) gas(i,j,k)%co	  = 0.0
	  if (gas(i,j,k)%zco2  .lt.0.0) gas(i,j,k)%zco2   = 0.0
	  if (gas(i,j,k)%ho    .lt.0.0) gas(i,j,k)%ho	  = 0.0
	  if (gas(i,j,k)%ho2   .lt.0.0) gas(i,j,k)%ho2    = 0.0
	  if (gas(i,j,k)%h2o2  .lt.0.0) gas(i,j,k)%h2o2   = 0.0
	  if (gas(i,j,k)%xno   .lt.0.0) gas(i,j,k)%xno    = 0.0
	  if (gas(i,j,k)%xno2  .lt.0.0) gas(i,j,k)%xno2   = 0.0
	  if (gas(i,j,k)%xno3  .lt.0.0) gas(i,j,k)%xno3   = 0.0
	  if (gas(i,j,k)%xn2o5 .lt.0.0) gas(i,j,k)%xn2o5  = 0.0
	  if (gas(i,j,k)%hno3  .lt.0.0) gas(i,j,k)%hno3   = 0.0
	  if (gas(i,j,k)%ch4   .lt.0.0) gas(i,j,k)%ch4    = 0.0
	  if (gas(i,j,k)%ch2o  .lt.0.0) gas(i,j,k)%ch2o   = 0.0
	  if (gas(i,j,k)%ch3o2h.lt.0.0) gas(i,j,k)%ch3o2h = 0.0
	  if (gas(i,j,k)%so2   .lt.0.0) gas(i,j,k)%so2    = 0.0
	  if (gas(i,j,k)%h2so4 .lt.0.0) gas(i,j,k)%h2so4  = 0.0
	  if (gas(i,j,k)%dms   .lt.0.0) gas(i,j,k)%dms    = 0.0
#endif

#ifdef AQCHEM_ENABLE
	  if (aqc(i,j,k)%o3    .lt.0.0) aqc(i,j,k)%o3    = 0.0
	  if (aqc(i,j,k)%civ   .lt.0.0) aqc(i,j,k)%civ   = 0.0
	  if (aqc(i,j,k)%h2o2  .lt.0.0) aqc(i,j,k)%h2o2  = 0.0
	  if (aqc(i,j,k)%xnv   .lt.0.0) aqc(i,j,k)%xnv   = 0.0
	  if (aqc(i,j,k)%siv   .lt.0.0) aqc(i,j,k)%siv   = 0.0
	  if (aqc(i,j,k)%svi   .lt.0.0) aqc(i,j,k)%svi   = 0.0
	  if (aqc(i,j,k)%ch2o  .lt.0.0) aqc(i,j,k)%ch2o  = 0.0
	  if (aqc(i,j,k)%ch3o2h.lt.0.0) aqc(i,j,k)%ch3o2h= 0.0
	  if (aqr(i,j,k)%o3    .lt.0.0) aqr(i,j,k)%o3    = 0.0
	  if (aqr(i,j,k)%civ   .lt.0.0) aqr(i,j,k)%civ   = 0.0
	  if (aqr(i,j,k)%h2o2  .lt.0.0) aqr(i,j,k)%h2o2  = 0.0
	  if (aqr(i,j,k)%xnv   .lt.0.0) aqr(i,j,k)%xnv   = 0.0
	  if (aqr(i,j,k)%siv   .lt.0.0) aqr(i,j,k)%siv   = 0.0
	  if (aqr(i,j,k)%svi   .lt.0.0) aqr(i,j,k)%svi   = 0.0
	  if (aqr(i,j,k)%ch2o  .lt.0.0) aqr(i,j,k)%ch2o  = 0.0
	  if (aqr(i,j,k)%ch3o2h.lt.0.0) aqr(i,j,k)%ch3o2h= 0.0
#endif

#ifdef SOLIDCHEM_ENABLE
	  if (solidi(i,j,k)%o3    .lt.0.0) solidi(i,j,k)%o3    = 0.0
	  if (solidi(i,j,k)%h2o2  .lt.0.0) solidi(i,j,k)%h2o2  = 0.0
	  if (solidi(i,j,k)%xnv   .lt.0.0) solidi(i,j,k)%xnv   = 0.0
	  if (solidi(i,j,k)%ch2o  .lt.0.0) solidi(i,j,k)%ch2o  = 0.0
	  if (solidi(i,j,k)%ch3o2h.lt.0.0) solidi(i,j,k)%ch3o2h= 0.0
	  if (solidi(i,j,k)%siv   .lt.0.0) solidi(i,j,k)%siv   = 0.0
	  if (solidi(i,j,k)%svi   .lt.0.0) solidi(i,j,k)%svi   = 0.0
#endif

	  if (xn_ccn(i,j,k).lt.0.0) xn_ccn(i,j,k) = 0.0
	  if (xn_in (i,j,k).lt.0.0) xn_in (i,j,k) = 0.0

	end do
	end do
	end do
	
!	data_max  = maxval (xn_ccn)
!	index_max = maxloc (xn_ccn)
!	data_min  = minval (xn_ccn)
!	index_min = minloc (xn_ccn)
	
!	write(6,*)"Max * min CCN =",data_max,index_max,data_min,index_min

!	data_max  = maxval (xn_in)
!	index_max = maxloc (xn_in)
!	data_min  = minval (xn_in)
!	index_min = minloc (xn_in)

!	write(6,*)"Max * min IN =",data_max,index_max,data_min,index_min

	close(103)
	close(104)
	close(105)
 
	avp  = p0*0.01
	ps   = 1.e5/2.87e2*p0**(1.005e3/2.87e2-1.)
        pp0  = 1.e5*p0**(1.005e3/2.87e2)
        den0 = pp0/(2.87e2*pt0*p0)

	avden0(1:nzi) = (den0(1:nzi) + den0(2:nz))*0.5
	avden0(nz)    =  den0(nz)

        pxx (1:nz) = (1.e5/pp0(1:nz))**(.286)
        pxxi(1:nz) = (1.e5/pp0(1:nz))**0.3

	do k=1,40
	  xn_ccn(1:nx,1:ny,k) = xn_ccn0*den0(1)/den0(k)
	  xn_in (1:nx,1:ny,k) = xn_in0 *den0(1)/den0(k)
!	  xn_ccn(1:nx,1:ny,k) = xn_ccn0
!	  xn_in (1:nx,1:ny,k) = xn_in0
	end do
	
	if(dt.ne.dtp.or.dtau.ne.dtaup)then
	  write(6,*)"Time step conflicated"
	  call exit(1)
	endif

	nstart = ntime

	write(6,'(1x/6x,"CM4.0 ",i2,"-dimensional ",a4," version","  Case: ",a20//)')ndims,cwarm,casename
	write(6,'(8x,"This is a run started with two-step-data (resume mode).")')
	write(6,'(10x,"The actual starting time is:",i5," hr ",i5," min ", i5," sec ")')nhour,nmin,nsec
	write(6,'(24x,"or, from step:",i8/)')ntime

#ifdef XPLBC_ENABLE
	  yes_no = 'X-Periodical'
#else
	  yes_no = 'X-Radiation'
#endif

#ifdef YPLBC_ENABLE
	  yes_no2 = 'Y-Periodical'
#else
	  yes_no2 = 'Y-Radiation'
#endif

	write(6,'(8x,"Basic time step = ",f8.2," second")')dt
	write(6,'("	   Small time step = ",f8.2," second")')dtau
	write(6,'("	   dx & dy         = ",2f7.1," meter")')dx,dy
	write(6,'("	   Gridpoint number= ",3i7," in x, y, & z direction")')nx,ny,nz
	write(6,'("	   Lateral boundary: ",a12,"-type conditions")')yes_no
	write(6,'("			     ",a12,"-type conditions"/)')yes_no2

        write(6,'(8x,"In calculating CN nucleation rate Nc = Cs^k")')xnc0_c*1.e-6
        write(6,'("        C = ",f10.2," 1/cc;  k = ",f10.2/)')xnc0_k

        write(6,'(8x,"u_shift = ",f6.2," v_shift = ",f6.2," m/s")')u0shift,v0shift
        write(6,'("	   nsmooth(1/0 call smooth) start with ",i2)')nsmooth
        write(6,'("	   lbc relexation coefficient = ",f6.2/)')eps

! ==========================================
! ===        Print Info of Filter:       ===
!
        if(if_filter.eq.0)then
          write(6,'(8x,"No filter applied"/)')
        else
          xxx = real(if_filter)*dt/60.0
          write(6,'(8x,"Apply filter every", f7.2, " minutes"/)')xxx
        endif

        if(if_lbcdamp.ne.0)then
          write(6,'(8x,"Apply lbc-damp with tau_lbcdamp = ",i9," times of dtau"/)')if_lbcdamp
        else
          write(6,'(8x,"No lbc-damp applied"/)')
        endif

	write(6,'(8x,"Radiation are calculated at every ",f7.1," minute"/)')iradx*dt/60.0

#ifdef CHEM_ENABLE
        write(6,'(8x,"Chemical Reaction are calculated at every ",f7.1," minute"/)')i_creact*dt/60.0
#else
	write(6,'(8x,"Run without Chemistry "/)')
#endif

        write(6,'(8x,60("=")/)')

      return
       end
