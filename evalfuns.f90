!***********************************************************************
!***********************************************************************

subroutine sevalf(n,x,f,ind,inform)

	use globals, only: sF
	use myproblem

	implicit none

	! SCALAR ARGUMENTS
	integer, intent(in)  :: n,ind
	integer, intent(out) :: inform
	real(kind=8), intent(in) :: x(n)	
	real(kind=8), intent(out) :: f

	! FUNCTIONS
	logical :: IsANumber

	inform = 0 

	call evalf(n,x,f,ind)

	if ( .not. IsANumber(f) ) then		
		inform = -1
		write(*,1000) ind
	end if

	f = sF(ind) * f
                
 1000	format(/,1X,'WARNING: The objective function value ',   &
                'computed by the user-supplied',/,1X,'subroutine ', &
                'evalf may be +Inf, -Inf or NaN. Function number:',1X,I4)	
                
end subroutine sevalf

!***********************************************************************
!***********************************************************************

subroutine sevalg(n,x,g,ind,inform)

	use globals, only: sF
	use myproblem
	
	implicit none
	
	! SCALAR ARGUMENTS
	integer, intent(in)  :: n,ind
	integer, intent(out) :: inform
	
	! ARRAY ARGUMENTS
	real(kind=8), intent(in)  :: x(n)
	real(kind=8), intent(out) :: g(n)
	
	! LOCAL SCALAR
	integer :: i
	
	! FUNCTIONS
	logical :: IsANumber
  
	inform = 0 
	
	call evalg(n,x,g,ind)
	
	do i = 1,n
		if ( .not. IsANumber(g(i)) ) then
			inform = -1
			write(*,1000) ind,i
		end if
	end do
	
	g = sF(ind) * g
	
 1000	format(/,1X,'WARNING: There is an element whose value may be +Inf, -Inf or NaN in the gradient',/,10X,&
					'of the objective computed by the user-supplied subroutine evalg.',/,10X,& 
					'Function number:',1X,I4,2X,'Coordenate:',1X,I5)	
end subroutine sevalg

!***********************************************************************
!***********************************************************************

subroutine sevalh(n,x,H,ind,inform)

	use globals, only: sF
	use myproblem
	
	implicit none
	
	! SCALAR ARGUMENTS
	integer, intent(in)  :: n,ind
	integer, intent(out) :: inform
	
	! ARRAY ARGUMENTS
	real(kind=8), intent(in)  :: x(n)
	real(kind=8), intent(out) :: H(n,n)
	
	! LOCAL SCALAR
	integer :: i,j
	
	! FUNCTIONS
	logical :: IsANumber
  
	inform = 0 
	
	call evalh(n,x,H,ind)
	
	do i = 1,n
		do j = 1,i
			if ( .not. IsANumber(H(i,j)) ) then
				inform = -1
				write(*,1000) ind,i,j
			end if
		end do
	end do
	
	H = sF(ind) * H
	
 1000	format(/,1X,'WARNING: There is an element whose value may be +Inf, -Inf or NaN in the Hessian',/,10X,&
					'of the objective computed by the user-supplied subroutine evalh.',/,10X,& 
					'Function number:',1X,I4,2X,'Coordenate:',1X,I5,1X,',',1X,I5)	
end subroutine sevalh


!***********************************************************************
!***********************************************************************

subroutine evalphi(stp,phi,ind)

	use globals, only: n => nfinner, x => xinner, d

	implicit none

	! SCALAR ARGUMENTS
	integer, intent(in) :: ind
	real(kind=8), intent(in) :: stp
	real(kind=8), intent(out) :: phi
	
	! LOCAL SCALAR
	integer :: infofun

	call sevalf(n,x + stp * d,phi,ind,infofun)

end subroutine evalphi

!!***********************************************************************
!!***********************************************************************

logical function IsANumber(x)

  implicit none

  ! SCALAR ARGUMENTS
  real(kind=8), intent(in) :: x
  
  ! LOCAL SCALAR
  real(kind=8) :: bignum
  
  bignum = 1.0d+99

  IsANumber = .true.
  if ( .not. abs( x ) .le. bignum ) IsANumber = .false.

end function IsANumber
