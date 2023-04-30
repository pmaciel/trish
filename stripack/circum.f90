SUBROUTINE CIRCUM(V1, V2, V3, C, IER)
    INTEGER IER
    REAL V1(3), V2(3), V3(3), C(3)
!
!***********************************************************
!   This subroutine returns the circumcenter of a spherical
! triangle on the unit sphere:  the point on the sphere sur-
! face that is equally distant from the three triangle
! vertices and lies in the same hemisphere, where distance
! is taken to be arc-length on the sphere surface.
!
!
! On input:
!
!       V1,V2,V3 = Arrays of length 3 containing the Carte-
!                  sian coordinates of the three triangle
!                  vertices (unit vectors) in CCW order.
!
! The above parameters are not altered by this routine.
!
!       C = Array of length 3.
!
! On output:
!
!       C = Cartesian coordinates of the circumcenter unless
!           IER > 0, in which case C is not defined.  C =
!           (V2-V1) X (V3-V1) normalized to a unit vector.
!
!       IER = Error indicator:
!             IER = 0 if no errors were encountered.
!             IER = 1 if V1, V2, and V3 lie on a common
!                     line:  (V2-V1) X (V3-V1) = 0.
!             (The vertices are not tested for validity.)
!
!***********************************************************
!
    INTEGER I
    REAL CNORM, CU(3), E1(3), E2(3)
!
! Local parameters:
!
! CNORM = Norm of CU:  used to compute C
! CU =    Scalar multiple of C:  E1 X E2
! E1,E2 = Edges of the underlying planar triangle:
!           V2-V1 and V3-V1, respectively
! I =     DO-loop index
!
    DO 1 I = 1, 3
        E1(I) = V2(I) - V1(I)
        E2(I) = V3(I) - V1(I)
1       CONTINUE
!
! Compute CU = E1 X E2 and CNORM**2.
!
        CU(1) = E1(2) * E2(3) - E1(3) * E2(2)
        CU(2) = E1(3) * E2(1) - E1(1) * E2(3)
        CU(3) = E1(1) * E2(2) - E1(2) * E2(1)
        CNORM = CU(1) * CU(1) + CU(2) * CU(2) + CU(3) * CU(3)
!
! The vertices lie on a common line if and only if CU is
!   the zero vector.
!
        IF (CNORM .NE. 0.) THEN
!
!   No error:  compute C.
!
            CNORM = SQRT(CNORM)
            DO 2 I = 1, 3
                C(I) = CU(I) / CNORM
2               CONTINUE
                IER = 0
                ELSE
!
!   CU = 0.
!
                IER = 1
                END IF
                RETURN
            END
