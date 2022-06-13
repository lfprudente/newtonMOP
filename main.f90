program main

	use globals, only: problem,seed
	use myproblem, only: inip
	implicit none
	
	! LOCAL SCALARS
	integer :: solver,n,m,outiter,nfev,nhev,chol,inform
	logical :: scaleF,checkder
	real(kind=8) :: time
	real(kind=8) :: theta
	
	! LOCAL ARRAYS
	real(kind=8),allocatable :: x(:),l(:),u(:),f(:)
	logical,allocatable :: strconvex(:)
	
	problem = 'AP1'
	
	! Choose the solver to be used:

	! solver = 1: Newton method with safeguards
	! solver = 2: Newton-Gradient method with safeguards
	! solver = 3: Newton method without safeguards
	! solver = 4: Steepest Descent method
	! solver = 5: Scalarized Newton method
	
	solver = 1

	seed = 123456.0d0
	
	call inip(n,m,x,l,u,strconvex,scaleF,checkder)

	if ( solver == 1 ) call Newton(n,m,x,l,u,strconvex,scaleF,checkder,outiter,time,nfev,nhev,chol,theta,inform)
	if ( solver == 2 ) call NewtonGradient(n,m,x,l,u,scaleF,checkder,outiter,time,nfev,nhev,chol,theta,inform)
	if ( solver == 3 ) call NewtonWoSafeg(n,m,x,l,u,scaleF,checkder,outiter,time,nfev,theta,inform)
	if ( solver == 4 ) call SteepestDescent(n,m,x,l,u,scaleF,checkder,outiter,time,nfev,theta,inform)
	if ( solver == 5 ) call NewtonScalar(n,m,x,l,u,scaleF,checkder,outiter,time,nfev,nhev,chol,theta,inform)
	
	deallocate(x,l,u,strconvex)
	
end program main
