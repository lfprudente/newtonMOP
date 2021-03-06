! *****************************************************************
! *****************************************************************

subroutine tnls(n,nind,x,l,u,m,lambda,rho,equatn,linear,f,g,amax,d, &
     rbdnnz,rbdind,rbdtype,xp,fp,gp,lsinfo,inform)

  use modmachconst
  use modalgconst
  use modouttyp
  use problvlv, only: fcnt

  implicit none

  ! SCALAR ARGUMENTS
  integer,      intent(in)    :: m,n,nind,rbdnnz
  integer,      intent(inout) :: inform
  integer,      intent(out)   :: lsinfo
  real(kind=8), intent(in)    :: amax,f
  real(kind=8), intent(out)   :: fp

  ! ARRAY ARGUMENTS
  integer,      intent(in)  :: rbdind(nind)
  character,    intent(in)  :: rbdtype(nind)
  logical,      intent(in)  :: equatn(m),linear(m)
  real(kind=8), intent(in)  :: d(nind),g(nind),l(nind),lambda(m), &
                               rho(m),u(nind),x(nind)
  real(kind=8), intent(out) :: gp(n),xp(nind)

  ! (On return: gp must be the shrinked full-space gradient)

  ! This subroutine computes a line search in direction d.
  !
  ! lsinfo:
  !
  ! At the first trial:
  !
  ! 5: x + amax d is at the boundary and f(x + amax d) is smaller than f

  ! Extrapolation:
  !
  ! 2: Unbounded objective function?
  ! 4: beta-condition holds. No extrapolation is done
  ! 6: Maximum number of extrapolations reached
  ! 7: Similar consecutive projected points
  ! 8: Not-well-defined objective function
  ! 9: Functional value increases

  ! Backtracking:
  !
  ! 0: Armijo satisfied
  ! 1: Small step with functional value similar to the current one
  ! 2: Unbounded objective function?
  ! 3: Too small backtracking step. Wrong gradient?

  ! LOCAL SCALARS
  logical      :: boundary
  integer      :: i
  real(kind=8) :: alpha,dsupn,gtd,xsupn

  ! EXTERNAL SUBROUTINES
  external :: calcal,csetp

  ! ==================================================================
  ! ==================================================================

  ! Test Armijo condition and beta condition. Decide between accept
  ! the first trial, extrapolate or backtrack.

  ! ==================================================================
  ! ==================================================================

  ! ------------------------------------------------------------------
  ! Compute directional derivative, dsupn and xsupn
  ! ------------------------------------------------------------------

  gtd = 0.0d0
  dsupn = 0.0d0
  xsupn = 0.0d0
  do i = 1,nind
     gtd = gtd + g(i) * d(i)
     dsupn = max( dsupn, abs( d(i) ) )
     xsupn = max( xsupn, abs( x(i) ) )
  end do

  if ( iprintinn .ge. 6 ) then
     write(*, 100) xsupn,amax,dsupn
     write(10,100) xsupn,amax,dsupn
  end if

  ! ------------------------------------------------------------------
  ! Compute first trial (projected point)
  ! ------------------------------------------------------------------

  alpha = 1.0d0

  boundary = .false.
  do i = 1,nind
     xp(i) = x(i) + d(i)

     if ( xp(i) .lt. l(i) .or. xp(i) .gt. u(i) ) then
        boundary = .true.
        xp(i) = max( l(i), min( xp(i), u(i) ) )
     end if
  end do

  if ( amax .eq. 1.0d0 ) then
     boundary = .true.

     do i = 1,rbdnnz
        if ( rbdtype(i) .eq. 'L' ) then
           xp(rbdind(i)) = l(rbdind(i))
        else if ( rbdtype(i) .eq. 'U' ) then
           xp(rbdind(i)) = u(rbdind(i))
        end if
     end do
  end if

  call csetp(nind,xp,inform)
  if ( inform .ne. 0 ) return

  call calcal(nind,xp,m,lambda,rho,equatn,linear,fp,inform)
  if ( inform .ne. 0 ) return

  if ( .not. boundary ) then

     if ( iprintinn .ge. 6 ) then
        write(*, 110) fp,fcnt
        write(10,110) fp,fcnt
     end if

  else
     if ( iprintinn .ge. 6 ) then
        write(*, 120) fp,fcnt
        write(10,120) fp,fcnt
     end if
  end if

  ! ------------------------------------------------------------------
  ! The first trial is an interior point.
  ! ------------------------------------------------------------------

  if ( .not. boundary ) then

     if ( iprintinn .ge. 6 ) then
        write(*, 140)
        write(10,140)
     end if

     ! Armijo condition holds.

     if ( fp .le. f + alpha * gamma * gtd ) then

        if ( iprintinn .ge. 6 ) then
           write(*, 150)
           write(10,150)
        end if

        go to 1000

     end if

     ! Armijo condition does not hold. We will do backtracking.

     if ( iprintinn .ge. 6 ) then
        write(* ,180)
        write(10,180)
     end if

     go to 2000

  end if

  ! ------------------------------------------------------------------
  ! The first trial is at the boundary.
  ! ------------------------------------------------------------------

  if ( iprintinn .ge. 6 ) then
     write(*, 190)
     write(10,190)
  end if

  ! Function value is smaller than at the current point. We will
  ! extrapolate.

  if ( fp .lt. f ) then

     if ( iprintinn .ge. 6 ) then
        write(*, 200)
        write(10,200)
     end if

     go to 1000

  end if

  ! Discard the projected point and consider x + amax d

  if ( iprintinn .ge. 6 ) then
     write(*, 210)
     write(10,210)
  end if

  alpha = amax

  do i = 1,nind
     xp(i) = x(i) + alpha * d(i)
  end do

  do i = 1,rbdnnz
     if ( rbdtype(i) .eq. 'L' ) then
        xp(rbdind(i)) = l(rbdind(i))
     else if ( rbdtype(i) .eq. 'U' ) then
        xp(rbdind(i)) = u(rbdind(i))
     end if
  end do

  call csetp(nind,xp,inform)
  if ( inform .ne. 0 ) return

  call calcal(nind,xp,m,lambda,rho,equatn,linear,fp,inform)
  if ( inform .ne. 0 ) return

  if ( iprintinn .ge. 6 ) then
     write(*, 130) alpha,fp,fcnt
     write(10,130) alpha,fp,fcnt
  end if

  ! Function value is smaller than or equal to (or even just a little
  ! bit greater than) at the current point. Line search is over.

  if ( fp .le. f + macheps23 * abs(f) ) then

     if ( iprintinn .ge. 6 ) then
        write(*, 220)
        write(10,220)
     end if

     call calcnal(nind,xp,m,lambda,rho,equatn,linear,gp,inform)
     if ( inform .ne. 0 ) return

     lsinfo = 5

     if ( iprintinn .ge. 6 ) then
        write(*, 900)
        write(10,900)
     end if

     return

  end if

  ! Function value is greater than at the current point. We will
  ! do backtracking.

  if ( iprintinn .ge. 6 ) then
     write(*, 230)
     write(10,230)
  end if

  go to 2000

  ! ==================================================================
  ! ==================================================================

  ! Extrapolation

  ! ==================================================================
  ! ==================================================================

 1000 continue

  call extrapolation(n,nind,x,l,u,m,lambda,rho,equatn,linear,g,xp,fp,    &
       gp,d,alpha,amax,rbdnnz,rbdind,rbdtype,fmin,beta,etaext,maxextrap, &
       lsinfo,inform)

  return

  ! ==================================================================
  ! ==================================================================
  
  ! Backtracking
  
  ! ==================================================================
  ! ==================================================================
  
 2000 continue

  call backtracking(nind,x,m,lambda,rho,equatn,linear,f,d,gtd,alpha, &
       fp,xp,calcal,csetp,lsinfo,inform)
  if ( inform .ne. 0 ) return

  call calcnal(nind,xp,m,lambda,rho,equatn,linear,gp,inform)
  if ( inform .ne. 0 ) return

  ! ==================================================================
  ! End of backtracking
  ! ==================================================================

  ! NON-EXECUTABLE STATEMENTS

 100  format(/,5X,'TN Line search (xsupn = ',1P,D7.1,', amax = ',        &
                   1P,D7.1,', dsupn = ',1P,D7.1,')')
 110  format(  5X,'Unitary step    F = ',1P,D24.16,' FE = ',I7)
 120  format(  5X,'Projected point F = ',1P,D24.16,' FE = ',I7)
 130  format(  5X,'Alpha = ',1P,D7.1,' F = ',1P,D24.16,' FE = ',I7)
 140  format(  5X,'The first trial is an interior point.')
 150  format(  5X,'Armijo condition holds.')
 180  format(  5X,'Armijo condition does not hold. We will backtrack.')
 190  format(  5X,'The first trial is at the boundary.')
 200  format(  5X,'Function value at the boundary is smaller than at ',  &
                  'the current point.',/,5X,'We will extrapolate.')
 210  format(  5X,'Discarding projected point. We will now consider x ', &
                  '+ amax d.')
 220  format(  5X,'Function value at the boundary is smaller than or ',  &
                  'equal to than at the',/,5X,'current point. Line ',    &
                  'search is over.')
 230  format(  5X,'Function value at the boundary is greater than at ',  &
                  'the current point.',/,5X,'We will backtrack.')
 900  format(  5X,'Flag of TN Line search: First trial accepted.')

end subroutine tnls
