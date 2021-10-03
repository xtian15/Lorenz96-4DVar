SUBROUTINE lorenzN(xin,xout,nx,force,dt)
! Lorenz N model (1 scale)
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: nx
  REAL(KIND=8),INTENT(IN) :: xin(nx),force(nx),dt
  REAL(KIND=8),INTENT(OUT) :: xout(nx)
  INTEGER :: i

  xout(1) = xin(nx) * ( xin(2) - xin(nx-1) ) - xin(1) + force(1)
  xout(2) = xin(1) * ( xin(3) - xin(nx) ) - xin(2) + force(2)
  DO i=3,nx-1
    xout(i) = xin(i-1) * ( xin(i+1) - xin(i-2) ) - xin(i) + force(i)
  END DO
  xout(nx) = xin(nx-1) * ( xin(1) - xin(nx-2) ) - xin(nx) + force(nx)

  xout(:) = dt * xout(:)

  RETURN
END SUBROUTINE lorenzN

SUBROUTINE lorenzN_tlm(xin, xin_tl, xout, xout_tl,nx,force, force_tl,dt)
   ! Lorenz N model (1 scale)
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: nx
  REAL(KIND=8),INTENT(IN) :: xin(nx),force(nx),dt
  REAL(KIND=8),INTENT(IN) :: xin_tl(nx),force_tl(nx)
  REAL(KIND=8),INTENT(OUT) :: xout(nx)
  REAL(KIND=8),INTENT(OUT) :: xout_tl(nx)
  INTEGER :: i

  xout_tl(1) = xin_tl(nx) * ( xin(2) - xin(nx-1) ) &
      + xin(nx) * ( xin_tl(2) - xin_tl(nx-1) ) &
      - xin_tl(1) + force_tl(1)
  xout(1) = xin(nx) * ( xin(2) - xin(nx-1) ) - xin(1) + force(1)

  xout_tl(2) = xin_tl(1) * ( xin(3) - xin(nx) ) &
      + xin(1) * ( xin_tl(3) - xin_tl(nx) ) &
      - xin_tl(2) + force_tl(2)
  xout(2) = xin(1) * ( xin(3) - xin(nx) ) - xin(2) + force(2)

  DO i=3,nx-1
    xout_tl(i) = xin_tl(i-1) * ( xin(i+1) - xin(i-2) ) &
      + xin(i-1) * ( xin_tl(i+1) - xin_tl(i-2) ) &
      - xin_tl(i) + force_tl(i)
    xout(i) = xin(i-1) * ( xin(i+1) - xin(i-2) ) - xin(i) + force(i)
  END DO

  xout_tl(nx) = xin_tl(nx-1) * ( xin(1) - xin(nx-2) ) &
      + xin(nx-1) * ( xin_tl(1) - xin_tl(nx-2) ) &
      - xin_tl(nx) + force_tl(nx)
  xout(nx) = xin(nx-1) * ( xin(1) - xin(nx-2) ) - xin(nx) + force(nx)

  xout_tl(:)=dt* xout_tl(:)
  xout(:) = dt * xout(:)

end subroutine

SUBROUTINE lorenzN_adj(xin, xin_ad, xout_ad, nx, force_ad, dt)
   ! Lorenz N model (1 scale)
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: nx
  REAL(KIND=8),INTENT(IN) :: xin(nx), dt
  REAL(KIND=8),INTENT(INOUT) :: xin_ad(nx),force_ad(nx)
  REAL(KIND=8),INTENT(INOUT) :: xout_ad(nx)
  INTEGER :: i

  !f2py intent(in) xin, nx, dt
  !f2py intent(in) xin_ad, force_ad, xout_ad
  !f2py intent(out) xin_ad, force_ad, xout_ad

  xout_ad(:)=xout_ad(:)*dt

  force_ad(nx)=force_ad(nx) + xout_ad(nx)
  xin_ad(nx)  =xin_ad(nx) -xout_ad(nx)
  xin_ad(nx-2)=xin_ad(nx-2) -xin(nx-1) *xout_ad(nx)
  xin_ad(1)   =xin_ad(1) +xin(nx-1) *xout_ad(nx)
  xin_ad(nx-1)=xin_ad(nx-1) +( xin(1) - xin(nx-2) )*xout_ad(nx)
  xout_ad(nx) =0.

  do i=nx-1, 3, -1
    force_ad(i)=force_ad(i) + xout_ad(i)
    xin_ad(i)  =xin_ad(i) -xout_ad(i)
    xin_ad(i-2)=xin_ad(i-2) -xin(i-1) *xout_ad(i)
    xin_ad(i+1)=xin_ad(i+1) +xin(i-1) *xout_ad(i)
    xin_ad(i-1)=xin_ad(i-1) +( xin(i+1) - xin(i-2) )*xout_ad(i)
    xout_ad(i) =0.
  end do

  force_ad(2)=force_ad(2) +xout_ad(2)
  xin_ad(2)  =xin_ad(2) -xout_ad(2)
  xin_ad(nx) =xin_ad(nx) -xin(1)*xout_ad(2)
  xin_ad(3)  =xin_ad(3) +xin(1)*xout_ad(2)
  xin_ad(1)  =xin_ad(1) +( xin(3) - xin(nx) )*xout_ad(2)
  xout_ad(2) =0.

  force_ad(1) =force_ad(1) +xout_ad(1)
  xin_ad(1)   =xin_ad(1) -xout_ad(1)
  xin_ad(nx-1)=xin_ad(nx-1) -xin(nx)*xout_ad(1)
  xin_ad(2)   =xin_ad(2) +xin(nx)*xout_ad(1)
  xin_ad(nx)  =xin_ad(nx) +( xin(2) - xin(nx-1) )*xout_ad(1)
  xout_ad(1)  =0.

end subroutine

subroutine lz_rk4(xin, xout, nx, force, dt)
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: nx
  REAL(KIND=8),INTENT(IN) :: xin(nx),force(nx),dt
  REAL(KIND=8),INTENT(OUT) :: xout(nx)
  INTEGER :: i

  !f2py intent(in) nx, xin, force, dt
  !f2py intent(out) xout

  real(kind=8), dimension(nx) :: xtmp, xf1, xf2, xf3, xf4

  xtmp=xin
  call lorenzN(xtmp, xf1, nx, force, dt)
  xtmp=xin + 0.5*xf1
  call lorenzN(xtmp, xf2, nx, force, dt)
  xtmp=xin + 0.5*xf2
  call lorenzN(xtmp, xf3, nx, force, dt)
  xtmp=xin + xf3
  call lorenzN(xtmp, xf4, nx, force, dt)

  xout=xin + (xf1 + 2.*xf2 + 2.*xf3 + xf4)/6.

  return
end subroutine lz_rk4

subroutine lz_rk4_tlm(xin, xin_tl, xout, xout_tl, nx, force, force_tl, dt)
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: nx
  REAL(KIND=8),INTENT(IN) :: xin(nx),force(nx),dt
  REAL(KIND=8),INTENT(IN) :: xin_tl(nx),force_tl(nx)
  REAL(KIND=8),INTENT(OUT) :: xout(nx)
  REAL(KIND=8),INTENT(OUT) :: xout_tl(nx)
  INTEGER :: i

  !f2py intent(in) nx, xin, force, dt
  !f2py intent(out) xout

  real(kind=8), dimension(nx) :: xtmp, xf1, xf2, xf3, xf4
  real(kind=8), dimension(nx) :: xtmp_tl, xf1_tl, xf2_tl, xf3_tl, xf4_tl

  xtmp_tl=xin_tl
  xtmp=xin
  call lorenzN_tlm(xtmp, xtmp_tl, xf1, xf1_tl, nx, force, force_tl, dt)

  xtmp_tl=xin_tl + 0.5*xf1_tl
  xtmp=xin + 0.5*xf1
  call lorenzN_tlm(xtmp, xtmp_tl, xf2, xf2_tl, nx, force, force_tl, dt)

  xtmp_tl=xin_tl + 0.5*xf2_tl
  xtmp=xin + 0.5*xf2
  call lorenzN_tlm(xtmp, xtmp_tl, xf3, xf3_tl, nx, force, force_tl, dt)

  xtmp_tl=xin_tl + xf3_tl
  xtmp=xin + xf3
  call lorenzN_tlm(xtmp, xtmp_tl, xf4, xf4_tl, nx, force, force_tl, dt)

  xout_tl=xin_tl + (xf1_tl + 2.*xf2_tl + 2.*xf3_tl + xf4_tl)/6.
  xout=xin + (xf1 + 2.*xf2 + 2.*xf3 + xf4)/6.

  return
end subroutine lz_rk4_tlm

subroutine lz_rk4_adj(xin, xin_ad, xout_ad, nx, force, force_ad, dt)
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: nx
  REAL(KIND=8),INTENT(IN) :: xin(nx),force(nx),dt
  REAL(KIND=8),INTENT(INOUT) :: xin_ad(nx),force_ad(nx)
  REAL(KIND=8),INTENT(INOUT) :: xout_ad(nx)
  INTEGER :: i

  !f2py intent(in) xin, force, nx, dt
  !f2py intent(in)  xin_ad, xout_ad, force_ad
  !f2py intent(out) xin_ad, xout_ad, force_ad

  real(kind=8), dimension(nx) :: xtmp, xf1, xf2, xf3, xf4
  real(kind=8), dimension(nx) :: xtmp1, xtmp2, xtmp3
  real(kind=8), dimension(nx) :: xtmp_ad, xf1_ad, xf2_ad, xf3_ad, xf4_ad

  xf4_ad= 1./6.*xout_ad
  xf3_ad= 2./6.*xout_ad
  xf2_ad= 2./6.*xout_ad
  xf1_ad= 1./6.*xout_ad
  xin_ad= xin_ad +xout_ad
  xout_ad=0.

  ! ----- forward calculation -----
  xtmp=xin
  call lorenzN(xtmp, xf1, nx, force, dt)
  xtmp=xin + 0.5*xf1
  xtmp1=xtmp
  call lorenzN(xtmp, xf2, nx, force, dt)
  xtmp=xin + 0.5*xf2
  xtmp2=xtmp
  call lorenzN(xtmp, xf3, nx, force, dt)
  xtmp=xin + xf3
  xtmp3=xtmp
  ! ----- end forward calculation -----

  xtmp_ad=0.
  call lorenzN_adj(xtmp3, xtmp_ad, xf4_ad, nx, force_ad, dt)
  xf3_ad=xf3_ad+xtmp_ad
  xin_ad=xin_ad+xtmp_ad
  xtmp_ad=0.

  call lorenzN_adj(xtmp2, xtmp_ad, xf3_ad, nx, force_ad, dt)
  xf2_ad=xf2_ad+0.5*xtmp_ad
  xin_ad=xin_ad+xtmp_ad
  xtmp_ad=0.

  call lorenzN_adj(xtmp1, xtmp_ad, xf2_ad, nx, force_ad, dt)
  xf1_ad=xf1_ad+0.5*xtmp_ad
  xin_ad=xin_ad+xtmp_ad
  xtmp_ad=0.

  call lorenzN_adj(xin, xtmp_ad, xf1_ad, nx, force_ad, dt)
  xin_ad=xin_ad+xtmp_ad
  xtmp_ad=0.

  return
end subroutine lz_rk4_adj
