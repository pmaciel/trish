LOGICAL FUNCTION LEFT(X1, Y1, Z1, X2, Y2, Z2, X0, Y0, Z0)
    REAL X1, Y1, Z1, X2, Y2, Z2, X0, Y0, Z0
!
!***********************************************************
!   This function determines whether node N0 is in the
! (closed) left hemisphere defined by the plane containing
! N1, N2, and the origin, where left is defined relative to
! an observer at N1 facing N2.
!
!
! On input:
!
!       X1,Y1,Z1 = Coordinates of N1.
!
!       X2,Y2,Z2 = Coordinates of N2.
!
!       X0,Y0,Z0 = Coordinates of N0.
!
! Input parameters are not altered by this function.
!
! On output:
!
!       LEFT = TRUE if and only if N0 is in the closed
!              left hemisphere.
!
! Modules required by LEFT:  None
!
!***********************************************************
!
! LEFT = TRUE iff <N0,N1 X N2> = det(N0,N1,N2) .GE. 0.
!
    LEFT = X0 * (Y1 * Z2 - Y2 * Z1) - Y0 * (X1 * Z2 - X2 * Z1) + Z0 * (X1 * Y2 - X2 * Y1) .GE. 0.
    RETURN
END
