subroutine SteepestDescent(n,m,x,l,u,scaleF,checkder,outiter,time,nfev,theta,inform)

	use globals, only: sF,nfinner,xinner,d,JF

	implicit none
	
	! SCALARS ARGUMENTS
	integer,intent(in)       :: n,m
	logical,intent(in)       :: scaleF,checkder
	integer,intent(out)      :: outiter,nfev,inform
	real(kind=8),intent(out) :: time
	real(kind=8),intent(out) :: theta
	
	! ARRAY ARGUMENTS
	real(kind=8),intent(inout) :: x(n)
	real(kind=8),intent(in)    :: l(n),u(n)
	
	! LOCAL SCALARS
	integer :: ind,j,infoIS,infoLS,infofun,allocerr,maxoutiter,&
	           nfevLS,ngev		 
	real(kind=8) :: stp,stpmin,epsopt,ftol,eta,q,norm2dSD
	real(kind=8) :: timeini,timefin
	logical :: iprint,iprintLS

	! LOCAL ARRAYS
	real(kind=8) :: g(n),lambda(m),gphi0(m),Fx(m),C(m)
	
	! EXTERNAL SUBROUTINES
	external :: evalphi
	
  ! Given a  continuous differentiable function F:R^n->R^m, this 
  ! subroutine solves the multiobjective optimization problem
  !
  ! Minimize F(x) subject to x in R^n,
  !
  ! using the Steepest Descent Method described in
  !
  ! M. L. N. Gonçalves, F. S. Lima, and L. F. Prudente, Globally 
  ! convergent Newton-type methods for multiobjective optimization, 
  ! technical report, 2022.
  ! 
  ! The user is required to code the function values and their gradients
  ! in the subroutine "myproblem".
  !
  ! The steepest descent direction is calculated in innersolver.f90 
  ! routine. This routine uses the software Algencan contained in 
  ! libalgencan.a. A step size satifying the nonmonotone Armijo 
  ! condition is computed at each iteration by means of the line search 
  ! routine armijo.f90. 
  !
  ! inform:
  !
  !  0: Optimality satisfied
  !  1: Maximum number of iterations reached
  !  2: An error occurred during the execution of the inner algorithm
  ! -1: An error occurred during the execution of the main algorithm
  !
  ! June 2022.
  !
  
  	!---------------------------------------------------------------------   
	!     Initialization
	!---------------------------------------------------------------------
		
	! Check derivatives
	
	if ( checkder ) call checkdF(n,m,x,l,u)
	
	! Start timing
		
	call cpu_time(timeini)
	
	! Set default parameters 

	epsopt = 5.0d0 * sqrt( 2.0d0 ** (-52) )
	ftol   = 1.d-04
	eta    = 0.85d0
	q      = 1.0d0
	maxoutiter = 2000
	stpmin     = 1.0d-15
		
	! Allocating the arrays
	
	allocate(d(n),xinner(n),sF(m),JF(m,n),stat=allocerr)
	if ( allocerr .ne. 0 ) then
		write(*,*) 'Allocation error in main program'
		stop
	end if
	
	! Saving the inner informations
	
	nfinner = n
	xinner  = x
	
	! Print problem information	
	
	iprint   = .true.	
	iprintLS = .false.
	
	if ( iprint ) then
		write(*,1000)
		write(*,1001)
		write(*,1010) n, m
	end if
	
	! Scale problem
	
	call scalefactor(n,m,x,scaleF)	

	! Counters
	
	outiter = 0	
	nfev    = 0
	ngev    = 0
		
	!------------------------------------------------------------------- 
	!     Main loop
	!-------------------------------------------------------------------    

	Main_loop: do
	
		!---------------------------------------------------------------
		!     Prepare the iteration
		!---------------------------------------------------------------
		
		! Compute the Jacobian JF
				
		do ind = 1,m
			call sevalg(n,x,g,ind,infofun)
			ngev = ngev + 1
			
			JF(ind,:) = g
			
			if ( infofun /= 0 ) then
				inform = -1
				deallocate(d,xinner,sF,JF)
				return
			end if
		end do
		
		! Compute the steepest descent direction
				
		call evaldSD(n,m,l,u,d,lambda,infoIS)
				
		norm2dSD = norm2(d)
		
		! Compute theta
		
		if ( infoIS /= 0 ) then
			theta = 1.0d99
		else
			theta = - norm2dSD ** 2 / 2.0d0
		end if
	
		! Print information

		if ( iprint ) then
			if ( outiter == 0 ) then
				write(*,1030) epsopt
				write(*,1040)
				write(*,1050) outiter,abs(theta),infoIS,nfev,ngev
			else
				if ( mod(outiter,10) == 0 ) write(*,1040)
				write(*,1060) outiter,abs(theta),infoLS,infoIS,nfev,ngev
			end if
		end if
	
		!-------------------------------------------------------------------
		!     Test stopping criteria
		!-------------------------------------------------------------------
		
		! Test optimality	
		
		if ( abs( theta ) <= epsopt  ) then
	
			inform = 0
			
			! Stop timing
			
			call cpu_time(timefin)
			
			time = timefin - timeini
  
			if ( iprint ) write(*,1070) nfev, ngev, time
			deallocate(d,xinner,sF,JF)
			return
			
		end if		
		
		! Test whether the number of iterations is exhausted
		
		if ( outiter >= maxoutiter ) then
				inform = 1
				
				! Stop timing
			
				call cpu_time(timefin)
			
				time = timefin - timeini
					
				if ( iprint ) write(*,1080) nfev, ngev, time

				deallocate(d,xinner,sF,JF)
				return
		end if
		
		! Test for errors in the inner solver
		
		if (  infoIS /= 0 ) then
				inform = 2
				
				! Stop timing
			
				call cpu_time(timefin)
			
				time = timefin - timeini
					
				if ( iprint ) write(*,1090) nfev, ngev, time

				deallocate(d,xinner,sF,JF)
				return
		end if

		!-------------------------------------------------------------------
		!     Iteration
		!-------------------------------------------------------------------
		
		! Increment outiter
				
		outiter = outiter + 1
						
		! Compute initial step
		
		stp = 1.0d0
				
		! Compute C
		
		if ( outiter == 1 ) then 
			do ind = 1,m
				call sevalf(n,x,C(ind),ind,infofun)
				nfev = nfev + 1   				
			end do
		else
			q = eta * q  + 1.0d0
			C = ( q - 1.0d0 ) / q * C + 1.0d0 / q * Fx
		end if
		
		! Compute gphi0 = [gphi_1(0),...,gphi_m(0)]
		
		gphi0 = matmul( JF, d )			
		
		! Compute the stepsize satisfying the Wolfe conditions 
											  
		call armijo(evalphi,stp,stpmin,m,C,gphi0,Fx,ftol,iprintLS,  &
		            nfevLS,infoLS) 
		            					
		! Check for errors in lsvecopt
		
		if ( infoLS <= - 1 ) then
			inform = - 1
			
			! Stop timing
			
			call cpu_time(timefin)
		
			time = timefin - timeini
			
			if ( iprint ) write(*,1100) nfev, ngev, time

			deallocate(d,xinner,sF,JF)
			return
		end if
				 
		nfev = nfev + nfevLS
		
		! Update x	
				
		x = x + stp * d
				
		! Saving the inner informations
		
		xinner  = x
		
	!---------------------------------------------------------------------
	!     Iterate
	!---------------------------------------------------------------------    

	end do Main_loop

	!--------------------------------------------------------------------- 
	!     End of main loop
	!---------------------------------------------------------------------
	
	! Non-executable statements
	
	1000 format(/,1X,'-------------------------------------------------------------------------------',/,1X,& 
					 '             Newton-type method for Multiobjective Optimization                ',/,1X,& 
					 '-------------------------------------------------------------------------------')
					 
	1001 format(/,1X,'Solver: Steepest descent method')
					 								 
	1010 format(/,1X,'Number of variables:',1X,I6,/,1X,'Number of functions:',1X,I6)	
	
	1030 format(/,1X,'Optimality tolerance:',1X,1P,D7.1)
	
	1040 format(/,4X,'out',2X,'|theta|',3X,'LS',1X,'IS',1X,'#evalf',2X,'#evalg')
	
	1050 format(1X,I5,3X,1P,D8.2,3X,'-',2X,I1,1X,I6,2X,I6)
	
	1060 format(1X,I5,3X,1P,D8.2,3X,I1,2X,I1,1X,I6,2X,I6)		
	
	1070 format(/,1X,'Flag of MOPsolver: solution was found',/,/,1X, &
					 'Number of functions evaluations:               ',I6,/,1X, &
					 'Number of derivatives evaluations:             ',I6/,/,1X, &
					 'Total CPU time in seconds: ',0P,F8.2)						
	
	1080 format(/,1X,'Flag of MOPsolver: maximum of iterations reached',/,/,1X,&
					 'Number of functions evaluations:               ',I6,/,1X, &
					 'Number of derivatives evaluations:             ',I6/,/,1X, &
					 'Total CPU time in seconds: ',0P,F8.2)
					 
	1090 format(/,1X,'Flag of MOPsolver: it was not possible to solve the subproblem',/,/,1X,&
					 'Number of functions evaluations:               ',I6,/,1X, &
					 'Number of derivatives evaluations:             ',I6/,/,1X, &
					 'Total CPU time in seconds: ',0P,F8.2)
					 
	1100 format(/,1X,'Flag of MOPsolver: an error occurred during the linesearch',/,/,1X,&
					 'Number of functions evaluations:               ',I6,/,1X, &
					 'Number of derivatives evaluations:             ',I6/,/,1X, &
					 'Total CPU time in seconds: ',0P,F8.2)
					 
end subroutine SteepestDescent
