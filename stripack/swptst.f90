LOGICAL FUNCTION SWPTST(N1, N2, N3, N4, X, Y, Z)
    INTEGER N1, N2, N3, N4
    REAL X(*), Y(*), Z(*)
!
!***********************************************************
!   This function decides whether or not to replace a
! diagonal arc in a quadrilateral with the other diagonal.
! The decision will be to swap (SWPTST = TRUE) if and only
! if N4 lies above the plane (in the half-space not contain-
! ing the origin) defined by (N1,N2,N3), or equivalently, if
! the projection of N4 onto this plane is interior to the
! circumcircle of (N1,N2,N3).  The decision will be for no
! swap if the quadrilateral is not strictly convex.
!
!
! On input:
!
!       N1,N2,N3,N4 = Indexes of the four nodes defining the
!                     quadrilateral with N1 adjacent to N2,
!                     and (N1,N2,N3) in counterclockwise
!                     order.  The arc connecting N1 to N2
!                     should be replaced by an arc connec-
!                     ting N3 to N4 if SWPTST = TRUE.  Refer
!                     to Subroutine SWAP.
!
!       X,Y,Z = Arrays of length N containing the Cartesian
!               coordinates of the nodes.  (X(I),Y(I),Z(I))
!               define node I for I = N1, N2, N3, and N4.
!
! Input parameters are not altered by this routine.
!
! On output:
!
!       SWPTST = TRUE if and only if the arc connecting N1
!                and N2 should be swapped for an arc con-
!                necting N3 and N4.
!
! Modules required by SWPTST:  None
!
!***********************************************************
!
    REAL DX1, DX2, DX3, DY1, DY2, DY3, DZ1, DZ2, DZ3, X4, Y4, Z4
!
! Local parameters:
!
! DX1,DY1,DZ1 = Coordinates of N4->N1
! DX2,DY2,DZ2 = Coordinates of N4->N2
! DX3,DY3,DZ3 = Coordinates of N4->N3
! X4,Y4,Z4 =    Coordinates of N4
!
    X4 = X(N4)
    Y4 = Y(N4)
    Z4 = Z(N4)
    DX1 = X(N1) - X4
    DX2 = X(N2) - X4
    DX3 = X(N3) - X4
    DY1 = Y(N1) - Y4
    DY2 = Y(N2) - Y4
    DY3 = Y(N3) - Y4
    DZ1 = Z(N1) - Z4
    DZ2 = Z(N2) - Z4
    DZ3 = Z(N3) - Z4
!
! N4 lies above the plane of (N1,N2,N3) iff N3 lies above
!   the plane of (N2,N1,N4) iff Det(N3-N4,N2-N4,N1-N4) =
!   (N3-N4,N2-N4 X N1-N4) > 0.
!
    SWPTST = DX3 * (DY2 * DZ1 - DY1 * DZ2) - DY3 * (DX2 * DZ1 - DX1 * DZ2) + DZ3 * (DX2 * DY1 - DX1 * DY2) .GT. 0.
    RETURN
END
