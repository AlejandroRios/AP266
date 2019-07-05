! fortran3 . F90
! using modules
program main


  ! MODULES
  use llt_module, only : induced_vel
  use llt_module, only : get_AIC_matrix
  use llt_module, only : get_residuals
  use llt_module, only : get_functions

  ! DECLARATIONS
  implicit none
  integer :: n_vort, ii,jj
  real , dimension(3) :: Uinf, X1, X2, Xp, gamma_ind, Gama, alpha0,chords,cl0,cla,res
  real , dimension(:,:), allocatable :: X
  real , dimension(:,:,:), allocatable :: AIC
  real  :: rho, Sref, CL, CD, L , D

  ! EXECUTION induced vel vector
  Uinf = (/1.0, 0.0, 0.0/)
  X1 = (/0.0, 2.0, 0.0/)
  X2 = (/0.0, 3.0, 0.0/)
  Xp = (/0.0, 1.5, 0.0/)
  call induced_vel(Uinf, X1, X2, Xp, gamma_ind)
  ! print *,'gamma iduz: ', gamma_ind

  ! EXECUTION AIC matrix
  ! Number of vortex
  n_vort = 3

  ! Allocation of allocatable matix X and AIC
  allocate(X(3,n_vort+1),AIC(3,n_vort,n_vort))
  ! X = reshape((/0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 2.0, 0.0, 0.0, 3.0, 0.0/),(/3,n_vort+1/))
  X = reshape((/0.0, 0.0, 0.2, 0.0, 1.0, 0.0, 0.0, 2.0, 0.0, 0.0, 3.0, 0.2/),(/3,n_vort+1/))

  ! get AIC matix
  call get_AIC_matrix(n_vort, X, Uinf, AIC)


  ! EXECUTION residuals
  Gama = (/1.0, 2.0, 1.0/)
  alpha0 = (/0.0, 0.1, 0.0/)
  chords = (/0.3, 0.6, 0.3/)
  cl0 = (/0.0, 0.0, 0.0/)
  cla = (/6.283, 6.283, 6.283/)
  rho = 1.0
  call get_residuals(n_vort, X, Gama, alpha0, chords, cl0, cla, Uinf, rho, res)
  print *,'Residuals : ', res

  ! EXECUTION get get_functions
  call get_functions(n_vort, X, Gama, alpha0, chords, cl0, cla, Uinf, rho, Sref, CL, CD, L, D)
  print *,'Sref : ', Sref
  print *,'CL : ', CL
  print *,'CD : ', CD
  print *,'L : ', L
  print *,'D : ', D


  deallocate(X,AIC)
end program main
