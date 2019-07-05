! fortran3 . F90
! using modules
program main


  ! MODULES and subroutines to be used
  use llt_module, only : induced_vel
  use llt_module, only : get_AIC_matrix

  ! DECLARATIONS
  implicit none
  integer :: n_vort, ii,jj
  real , dimension(3) :: Uinf, X1, X2, Xp, gamma_ind
  real , dimension(:,:), allocatable :: X
  real , dimension(:,:,:), allocatable :: AIC

  ! EXECUTION induced vel vector
  Uinf = (/1.0, 0.0, 0.0/)
  X1 = (/0.0, 2.0, 0.0/)
  X2 = (/0.0, 3.0, 0.0/)
  Xp = (/0.0, 1.5, 0.0/)
  call induced_vel(Uinf, X1, X2, Xp, gamma_ind)
  print *,'gamma iduz: ', gamma_ind

  ! EXECUTION AIC matrix

  ! Number of vortex to be used
  n_vort = 3

  ! Allovtion of X and AIC matrix
  allocate(X(3,n_vort+1),AIC(3,n_vort,n_vort))
  ! X = reshape((/0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 2.0, 0.0, 0.0, 3.0, 0.0/),(/3,n_vort+1/))
  X = reshape((/0.0, 0.0, 0.0, 0.0, 1.0, 0.1, 0.0, 2.0, 0.2, 0.0, 3.0, 0.3/),(/3,n_vort+1/))

  ! AIC matrix calculation
  call get_AIC_matrix(n_vort, X, Uinf, AIC)

  ! Print AIC Matrixx
  print *,'AIC matrix'
  do ii = 1,3
      print *,AIC(1,ii,:)
  end do
  deallocate(X,AIC)
end program main
