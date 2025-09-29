!**********************************************************!
!
! THIS FILE HAS BEEN ADAPTED FROM I-PI (i-pi/tools/f90)
!
!**********************************************************!

subroutine UpdateAnglePDF(gTheta, theta_idx, atxyz, theta_min, theta_max, cell, invCell, nframes, natoms, nangles, nbins)
    !
    ! Calculate the distribution of internal angles for 
    !
    IMPLICIT NONE
    !
    ! Input
    !
    INTEGER, INTENT(in) :: nframes, natoms, nangles, nbins, theta_idx(3*nangles)
    REAL(8), INTENT(in) :: theta_min, theta_max
    REAL(8), INTENT(in) :: atxyz(nframes,3*natoms), cell(3,3), invCell(3,3)
    !
    ! Output
    !
    REAL(8), INTENT(inout) :: gTheta(nbins,2)
    !f2py intent(hide) :: nframes
    !f2py intent(hide) :: nangles
    !f2py intent(hide) :: natoms
    !f2py intent(hide) :: nbins
    ! Local variables
    !
    INTEGER :: i_frame, i_angle, i_bin, i_A, i_B1, i_B2
    REAL(8) :: delta_theta, rAB1(3), rAB2(3), dAB1, dAB2, cos_theta, theta, norm
    !
    ! Histogram step initialization
    !
    delta_theta = gTheta(2,1) - gTheta(1,1)
    !
    ! Start computing g(θ) from MD configurations
    !
    norm=1.D0 / DBLE(nangles)
    DO i_frame=1,nframes
      DO i_angle=1,nangles
        i_A = theta_idx(3*i_angle-2)     ! index of central atom
        i_B1 = theta_idx(3*i_angle-1)  ! index of first edge atom
        i_B2 = theta_idx(3*i_angle)  ! index of second edge atom
        ! Compute the displacements from the central atom according to the minimum-image convention 
        CALL CalcMICDisplacement( & 
          cell, invCell, &
          atxyz(i_frame,3*i_B1-2:3*i_B1), &
          atxyz(i_frame,3*i_A -2:3*i_A ), &
          rAB1, dAB1)
        CALL CalcMICDisplacement( & 
          cell, invCell, &
          atxyz(i_frame,3*i_B2-2:3*i_B2), &
          atxyz(i_frame,3*i_A -2:3*i_A ), &
          rAB2, dAB2)
        ! Compute the angle:
        cos_theta = DOT_PRODUCT(rAB1, rAB2) / (dAB1 * dAB2)
        theta = ACOS(cos_theta)
        ! Screen angles that are outside the desired range
        IF (theta.LT.theta_max.AND.theta.GT.theta_min) THEN
          i_bin=INT((theta-theta_min)/delta_theta)+1  !bin/histogram position
          gTheta(i_bin,2) = gTheta(i_bin,2) + norm
        END IF
      END DO !i_angle 
    END DO !i_frame
    ! Normalize the histogram
END SUBROUTINE UpdateAnglePDF

SUBROUTINE CalcMICDisplacement(cell,invCell,rA,rB,rAB,dAB)
    !
    IMPLICIT NONE
    !
    REAL(8), INTENT(IN) :: rA(3), rB(3), invCell(3,3), cell(3,3)
    REAL(8), INTENT(OUT) :: rAB(3), dAB
    ! Local
    REAL(8) :: sAB(3)

    !Initialization of distance
    dAB=0.0
    !
    ! Compute distance between atom A and atom B (according to the minimum
    ! image convention)...
    !
    rAB = rA-rB
    !
    ! Displacement vector in fractional coordinates, s_AB = h^-1 r_AB 
    !
    sAB = MATMUL(invCell,rAB)
    !
    ! Impose the minimum-image convention
    !
    sAB = sAB - IDNINT(sAB)
    !
    ! Convert back to Cartesian coordinates
    !
    rAB = MATMUL(cell, sAB)
    dAB=DSQRT( DOT_PRODUCT(rAB, rAB) )
END SUBROUTINE CalcMICDisplacement
