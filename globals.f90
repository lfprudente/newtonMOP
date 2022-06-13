module globals

implicit none

character(len=10) :: problem
integer :: nfinner
real(kind=8) :: seed,xbi
real(kind=8),allocatable :: xinner(:),d(:),sF(:),JF(:,:),B(:,:,:)

end module globals
