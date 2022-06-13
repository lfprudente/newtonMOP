subroutine Newton(n,m,x,l,u,strconvex,scaleF,checkder,outiter,time,nfev,nhev,chol,theta,inform)

	use globals, only: sF,nfinner,xinner,d,JF,B

	implicit none
	
	! SCALARS ARGUMENTS
	integer,intent(in)       :: n,m
	logical,intent(in)       :: scaleF,checkder
	integer,intent(out)      :: outiter,nfev,nhev,chol,inform
	real(kind=8),intent(out) :: time
	real(kind=8),intent(out) :: theta
	
	! ARRAY ARGUMENTS
	real(kind=8),intent(inout) :: x(n)
	real(kind=8),intent(in)    :: l(n),u(n)
	logical,intent(in)         :: strconvex(m)
	
	! LOCAL SCALARS
	integer :: ind,j,infoIS,infoISSD,infoLS,infofun,allocerr,maxoutiter,&
	           nfevLS,infoChol,ngev				 
	real(kind=8) :: stp,stpmin,epsopt,sqrtepsopt,ftol,eta,q,mu,multmu,&
					norm2d,norm2dmu,gamma1,gamma2,minB,thetaSD
	real(kind=8) :: timeini,timefin
	logical :: iprint,iprintLS

	! LOCAL ARRAYS
	real(kind=8) :: g(n),lambda(m),dmu(n),gphi0(m),Fx(m),C(m),Btmp(n,n),&
	                dSD(n),lambdaSD(m)
	
	! EXTERNAL SUBROUTINES
	external :: evalphi
	
  ! Given a  continuous differentiable function F:R^n->R^m, this 
  ! subroutine solves the multiobjective optimization problem
  !
  ! Minimize F(x) subject to x in R^n,
  !
  ! using the Newton Method with Safeguards described in
  !
  ! M. L. N. Gonçalves, F. S. Lima, and L. F. Prudente, Globally 
  ! convergent Newton-type methods for multiobjective optimization, 
  ! technical report, 2022.
  ! 
  ! The user is required to code the function values, their gradients
  ! and Hessians in the subroutine "myproblem".
  !
  ! The Newton direction is calculated in innersolver.f90 routine. This 
  ! routine uses the software Algencan contained in libalgencan.a. A 
  ! stepsize satifying the nonmonotone Armijo condition is computed at 
  ! each iteration by means of the line search routine armijo.f90. 
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

	epsopt     = 5.0d0 * sqrt( 2.0d0 ** (-52) )
	sqrtepsopt = sqrt(epsopt)
	ftol       = 1.d-04
	eta        = 0.85d0
	q          = 1.0d0
	multmu     = 2.0d0
	gamma1     = 1.0d-6
	gamma2     = 1.0d-2
	maxoutiter = 2000
	stpmin     = 1.0d-15
		
	! Allocating the arrays
	
	allocate(d(n),xinner(n),sF(m),JF(m,n),B(n,n,m),stat=allocerr)
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
	nhev    = 0
	chol    = 0
	
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
				deallocate(d,xinner,sF,JF,B)
				return
			end if
		end do
		
		! Compute a positive definite approximations of the Hessians
		
		do ind = 1,m
			call sevalh(n,x,B(:,:,ind),ind,infofun)
			nhev = nhev + 1
			
			if ( infofun /= 0 ) then
				inform = -1
				deallocate(d,xinner,sF,JF,B)
				return
			end if
			
			if ( .not. strconvex(ind) ) then
			
				minB = B(1,1,ind)
				do j = 2,n
					minB = min( minB, B(j,j,ind) )
				end do
				
				if ( minB > 0.0d0 ) then
					mu = 0.0d0
				else
					mu = - minB + 1.0d0
				end if
								
				do	
					Btmp(:,:) = B(:,:,ind)
					forall (j = 1:n) Btmp(j,j) = B(j,j,ind) + mu
					
					! Try to compute the Cholesky factorization of Btmp
					
					call DPOTRF('L',n,Btmp,n,infoChol)
					chol = chol + 1
					
					if ( infoChol == 0 ) then
						forall (j = 1:n) B(j,j,ind) = B(j,j,ind) + mu
						exit
					else
						mu = max( multmu * mu, 1.0d0 )
					end if
				end do
			end if
		end do
					
		! Compute the search direction
		
		mu = 0.0d0
				
		do

			call evaldN(n,m,l,u,d,lambda,infoIS)
					
			! Compute theta	
			
			if ( infoIS /= 0 ) then
				theta = 1.0d99
			else
				dmu(:) = 0.0d0
				do ind = 1,m
					 dmu = dmu + lambda(ind) * JF(ind,:)
				end do

				theta = - 0.5d0 * norm2( dmu )**2
				if ( abs( theta ) <= epsopt ) exit
				
			end if
		
			norm2d   = norm2(d)			
			norm2dmu = norm2(dmu)
			
			if ( maxval( matmul( JF,d ) ) <= - gamma1 * norm2d * norm2dmu ) then
				if ( norm2d < gamma2 * norm2dmu ) d = gamma2 * norm2dmu / norm2d * d
				exit
			else
				mu = max( multmu * mu, 1.0d0 )
				forall (ind = 1:m, j = 1:n) B(j,j,ind) = B(j,j,ind) + mu
			end if
		end do
					
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
		
		if (  abs( theta ) <= epsopt  ) then

			inform = 0
			
			! Stop timing
			
			call cpu_time(timefin)
			
			time = timefin - timeini
  
			if ( iprint ) write(*,1070) nfev, ngev, time
			deallocate(d,xinner,sF,JF,B)				
			return
			
		end if		
		
		! Test whether the number of iterations is exhausted
		
		if ( outiter >= maxoutiter ) then
				inform = 1
				
				! Stop timing
			
				call cpu_time(timefin)
			
				time = timefin - timeini
					
				if ( iprint ) write(*,1080) nfev, ngev, time

				deallocate(d,xinner,sF,JF,B)
				return
		end if
		
		! Test for errors in the inner solver
		
		if (  infoIS /= 0 ) then
				inform = 2
				
				! Stop timing
			
				call cpu_time(timefin)
			
				time = timefin - timeini
					
				if ( iprint ) write(*,1090) nfev, ngev, time

				deallocate(d,xinner,sF,JF,B)
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

			deallocate(d,xinner,sF,JF,B)
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

	1001 format(/,1X,'Solver: Newton with safeguards')
					 								 
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

end subroutine Newton
