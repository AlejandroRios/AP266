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
  use llt_module_b, only : GET_RESIDUALS_B
  use llt_module_b, only : GET_FUNCTIONS_B

  ! DECLARATIONS
  implicit none
  integer :: n_vort, ii,jj
  real , dimension(3) :: Uinf, X1, X2, Xp, gamma_ind, Gama, alpha0,chords,cl0,cla,res
  real , dimension(3) :: Gamad, alpha0d,chordsd, resd, dresdx_FD, dresdx_AD, resAD, resdAD, difference_res
  real , dimension(:,:), allocatable :: X, Xd, Xb
  ! real , dimension(:,:,:), allocatable :: AIC
  real  :: rho, Sref, CL, CD, L , D, h, Srefd, CLd, CDd, Ld, Dd, CDAD, CLAD, LAD, DAD, SrefAD
  real  :: dSrefdx_FD, dCLdx_FD, dCDdx_FD, dLdx_FD, dDdx_FD, dSrefdx_AD, dCLdx_AD, dCDdx_AD, dLdx_AD, dDdx_AD
  real :: SrefdAD, CLdAD, CDdAD, LdAD, DdAD
  real ::  difference_Sref,  difference_CL, difference_CD, difference_L, difference_D, dotprod
  real , dimension(3) :: chordsb, alpha0b, Gamab, resb, resbb
  real :: Srefb, CLb, CDb, Lb, Db
  real :: Srefbb, CLbb, CDbb, Lbb, Dbb


  ! EXECUTION induced vel vector
  Uinf = (/1.0, 0.0, 0.0/)

  X1 = (/0.0, 2.0, 0.0/)
  X2 = (/0.0, 3.0, 0.0/)
  Xp = (/0.0, 1.5, 0.0/)
  call induced_vel(Uinf, X1, X2, Xp, gamma_ind)
  ! print *,'gamma iduz: ', gamma_ind

  ! EXECUTION
  n_vort = 3
  allocate(X(3,n_vort+1),Xd(3,n_vort+1),Xb(3,n_vort+1))
  X = reshape((/0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 2.0, 0.0, 0.0, 3.0, 0.0/),(/3,n_vort+1/))
  Xd = reshape((/-0.1, 1.0, 0.3, 0.5, -0.7, 0.0, 0.2, 0.1, -0.6, 0.0, -0.9, 0.2/),(/3,n_vort+1/))
  gamab = 0
  alpha0b = 0
  chordsb = 0

  ! Input variables DECLARATIONS
  Gama = (/2.0, 2.0, 2.0/)
  alpha0 = (/0.0, 0.0, 0.0/)
  chords = (/0.3, 0.3, 0.3/)
  cl0 = (/0.0, 0.0, 0.0/)
  cla = (/6.283, 6.283, 6.283/)
  rho = 1.0

  ! Set of arbitrary derivative seeds DECLARATION
  gamad = (/-1.0, 0.5, 0.3/)
  alpha0d = (/0.0, 0.1, -0.1/)
  chordsd = (/1.0, -0.2, 0.7/)
  cl0 = (/0.0, 0.0, 0.0/)
  cla = (/6.283, 6.283, 6.283/)
  rho = 1.0


  !###############################################################################
  ! ! Resuduals from forward AD
  ! call GET_RESIDUALS_D(n_vort, x, xd, gama, gamad, alpha0, alpha0d, chords, chordsd, cl0, cla, Uinf, rho, resAD, resdAD)
  ! ! Derivative seeds of the outputs:
  ! print *,'resd: ', resdAD
  !
  ! ! Resuduals from reverse AD
  ! ! Set of arbitrary derivative seeds DECLARATION
  ! resb = (/0.1, -0.2, 0.3/)
  ! resbb = (/0.1, -0.2, 0.3/)
  ! call GET_RESIDUALS_B(n_vort, x, xb, gama, gamab&
  ! & , alpha0, alpha0b, chords, chordsb&
  ! &, cl0, cla, Uinf, rho&
  ! &, res, resbb)
  !
  ! ! Derivative seeds of the inputs:
  ! print *,'Xb'
  ! do ii = 1,3
  ! print *,Xb(:,ii)
  ! end do
  ! print *,'gamab: ', gamab
  ! print *,'alpha0b: ', alpha0b
  ! print *,'chordsb: ', chordsb
  !
  !
  ! ! dot product of get_residuals
  ! dotprod =  sum(Xd*Xb) + sum(Gamad*Gamab) + sum(alpha0d*alpha0b) + sum(chordsd*chordsb) - sum(resdAD*resb)
  ! print *,'dotprod : ', dotprod

 !###############################################################################
  ! Befor running this part of the program comment the last part of the code (line 68-95)

  ! Resuduals from get functions forward AD
  call GET_FUNCTIONS_D(n_vort, x, xd, gama&
&   , gamad, alpha0, alpha0d, chords&
&   , chordsd, cl0, cla, Uinf&
&   , rho, SrefAD, SrefdAD, CLAD&
&   , CLdAD, CDAD, CDdAD, LAD, LdAD, DAD, DdAD)

  !set of arbitrary derivative seeds:
  Srefbb = 0.6
  CLbb = -0.2
  CDbb = 0.1
  Lbb = -0.3
  Dbb = 0.7

  Srefb = 0.6
  CLb = -0.2
  CDb = 0.1
  Lb = -0.3
  Db = 0.7

  ! Resuduals from get functions reverse AD
  call GET_FUNCTIONS_B(n_vort, x, xb, gama, gamab, alpha0, alpha0b&
&   , chords, chordsb, cl0, cla, Uinf, rho, Sref, Srefbb, CL, CLbb, CD, &
&   CDbb, L, Lbb, D, Dbb)

  ! Printing results
  print *,'Srefd: ', SrefdAD
  print *,'CLd: ', CLdAD
  print *,'CDd: ', CDdAD
  print *,'Ld: ', LdAd
  print *,'Dd: ', DdAD

  print *,'Xb'
  do ii = 1,3
  print *,Xb(:,ii)
  end do

  print *,'gamab: ', gamab
  print *,'alpha0b: ', alpha0b
  print *,'chordsb: ', chordsb

  ! Dot Product calculation
  dotprod = sum(Xd*Xb) + sum(Gamad*Gamab) + sum(alpha0d*alpha0b) + sum(chordsd*chordsb)&
  & - SrefdAD*Srefb - CLdAD*CLb - CDdAD*CDb - LdAD*Lb - DdAD*Db
  print *,'dotprod: ', dotprod

end program main
