! fortran3 . F90
! using modules
program main


  ! MODULES
  use llt_module, only : induced_vel
  use llt_module, only : get_AIC_matrix
  use llt_module, only : get_residuals
  use llt_module, only : get_functions
  use llt_module_d, only : GET_RESIDUALS_D
  use llt_module_d, only : GET_FUNCTIONS_D

  ! DECLARATIONS
  implicit none
  integer :: n_vort, ii,jj
  real , dimension(3) :: Uinf, X1, X2, Xp, gamma_ind, Gama, alpha0,chords,cl0,cla,res
  real , dimension(3) :: Gamad, alpha0d,chordsd, resd, dresdx_FD, dresdx_AD, resAD, resdAD, difference_res
  real , dimension(:,:), allocatable :: X, Xd
  real , dimension(:,:,:), allocatable :: AIC
  real  :: rho, Sref, CL, CD, L , D, h, Srefd, CLd, CDd, Ld, Dd, CDAD, CLAD, LAD, DAD, SrefAD
  real  :: dSrefdx_FD, dCLdx_FD, dCDdx_FD, dLdx_FD, dDdx_FD, dSrefdx_AD, dCLdx_AD, dCDdx_AD, dLdx_AD, dDdx_AD
  real :: SrefdAD, CLdAD, CDdAD, LdAD, DdAD
  real ::  difference_Sref,  difference_CL, difference_CD, difference_L, difference_D


  ! Execution of iduced vel subroutine to obtain gama_ind
  Uinf = (/1.0, 0.0, 0.0/)

  X1 = (/0.0, 2.0, 0.0/)
  X2 = (/0.0, 3.0, 0.0/)
  Xp = (/0.0, 1.5, 0.0/)
  call induced_vel(Uinf, X1, X2, Xp, gamma_ind)

  ! Execution get residualds and get get_functions

  ! number of vortex
  n_vort = 3

  ! Matrix allocation X and Xd
  allocate(X(3,n_vort+1),Xd(3,n_vort+1))
  X = reshape((/0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 2.0, 0.0, 0.0, 3.0, 0.0/),(/3,n_vort+1/))
  Xd = reshape((/-0.1, 1.0, 0.3, 0.5, -0.7, 0.0, 0.2, 0.1, -0.6, 0.0, -0.9, 0.2/),(/3,n_vort+1/))

  ! Input variables DECLARATIONS

  Gama = (/2.0, 2.0, 2.0/)
  alpha0 = (/0.0, 0.0, 0.0/)
  chords = (/0.3, 0.3, 0.3/)
  cl0 = (/0.0, 0.0, 0.0/)
  cla = (/6.283, 6.283, 6.283/)
  rho = 1.0

  ! Calling get residuals and get functions t
  call get_residuals(n_vort, X, Gama, alpha0, chords, cl0, cla, Uinf, rho, res)
  call get_functions(n_vort, X, Gama, alpha0, chords, cl0, cla, Uinf, rho, Sref, CL, CD, L, D)


  ! value of h
  h = 1e-12

  ! Perturbation variables delcaration
  Gamad = (/-1.0, 0.5, 0.3/)
  alpha0d = (/0.0, 0.1, -0.1/)
  chordsd = (/1.0, -0.2, 0.7/)
  cl0 = (/0.0, 0.0, 0.0/)
  cla = (/6.283, 6.283, 6.283/)
  rho = 1.0

  ! Resuduals and functions from FD
  call get_residuals(n_vort, X + h*Xd, Gama + h*Gamad, alpha0 + h*alpha0d, chords + h*chordsd, cl0, cla, Uinf, rho, resd)
  call get_functions(n_vort, X + h*Xd, Gama + h*Gamad, alpha0 + h*alpha0d&
&   , chords + h*chordsd, cl0, cla, Uinf&
&   , rho, Srefd, CLd, CDd, Ld, Dd)

  !directional derivative FD
  dresdx_FD = (1/h)*[resd - res]
  dSrefdx_FD = (1/h)*(Srefd - Sref)
  dCLdx_FD = (1/h)*(CLd - CL)
  dCDdx_FD = (1/h)*(CDd - CD)
  dLdx_FD = (1/h)*(Ld - L)
  dDdx_FD = (1/h)*(Dd - D)


  ! Resuduals and functions from AD
  call GET_RESIDUALS_D(n_vort, x, xd, gama, gamad, alpha0, alpha0d, chords, chordsd, cl0, cla, Uinf, rho, resAD, resdAD)
  call GET_FUNCTIONS_D(n_vort, x, xd, gama&
&   , gamad, alpha0, alpha0d, chords&
&   , chordsd, cl0, cla, Uinf&
&   , rho, SrefAD, SrefdAD, CLAD&
&   , CLdAD, CDAD, CDdAD, LAD, LdAD, DAD, DdAD)
  !directional derivative AD
  dresdx_AD = resdAD
  dSrefdx_AD = SrefdAD
  dCLdx_AD = CLdAD
  dCDdx_AD = CDdAD
  dLdx_AD = LdAD
  dDdx_AD = DdAD
  print *,'Directional derivative AD: ',dDdx_AD

  ! Differences calculation
  difference_res = abs(dresdx_FD - dresdx_AD)
  difference_Sref = abs(dSrefdx_FD - dSrefdx_AD)
  difference_CL = abs(dCLdx_FD - dCLdx_AD)
  difference_CD = abs(dCDdx_FD - dCDdx_AD)
  difference_L = abs(dLdx_FD - dLdx_AD)
  difference_D = abs(dDdx_FD - dDdx_AD)

  print *,'Diference FD - AD: ', difference_D





  deallocate(X,Xd)
end program main
