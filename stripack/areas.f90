REAL FUNCTION AREAS(V1, V2, V3)
    REAL V1(3), V2(3), V3(3)
!
!***********************************************************
!   This function returns the area of a spherical triangle
! on the unit sphere.
!
!
! On input:
!
!       V1,V2,V3 = Arrays of length 3 containing the Carte-
!                  sian coordinates of unit vectors (the
!                  three triangle vertices in any order).
!                  These vectors, if nonzero, are implicitly
!                  scaled to have length 1.
!
! Input parameters are not altered by this function.
!
! On output:
!
!       AREAS = Area of the spherical triangle defined by
!               V1, V2, and V3 in the range 0 to 2*PI (the
!               area of a hemisphere).  AREAS = 0 (or 2*PI)
!               if and only if V1, V2, and V3 lie in (or
!               close to) a plane containing the origin.
!
!***********************************************************
!
    DOUBLE PRECISION A1, A2, A3, CA1, CA2, CA3, DV1(3), DV2(3), DV3(3), S12, S23, S31, U12(3), U23(3), U31(3)
    INTEGER I
!
! Local parameters:
!
! A1,A2,A3 =    Interior angles of the spherical triangle
! CA1,CA2,CA3 = cos(A1), cos(A2), and cos(A3), respectively
! DV1,DV2,DV3 = Double Precision copies of V1, V2, and V3
! I =           DO-loop index and index for Uij
! S12,S23,S31 = Sum of squared components of U12, U23, U31
! U12,U23,U31 = Unit normal vectors to the planes defined by
!                 pairs of triangle vertices
!
    DO 1 I = 1, 3
        DV1(I) = DBLE(V1(I))
        DV2(I) = DBLE(V2(I))
        DV3(I) = DBLE(V3(I))
1       CONTINUE
!
! Compute cross products Uij = Vi X Vj.
!
        U12(1) = DV1(2) * DV2(3) - DV1(3) * DV2(2)
        U12(2) = DV1(3) * DV2(1) - DV1(1) * DV2(3)
        U12(3) = DV1(1) * DV2(2) - DV1(2) * DV2(1)
!
        U23(1) = DV2(2) * DV3(3) - DV2(3) * DV3(2)
        U23(2) = DV2(3) * DV3(1) - DV2(1) * DV3(3)
        U23(3) = DV2(1) * DV3(2) - DV2(2) * DV3(1)
!
        U31(1) = DV3(2) * DV1(3) - DV3(3) * DV1(2)
        U31(2) = DV3(3) * DV1(1) - DV3(1) * DV1(3)
        U31(3) = DV3(1) * DV1(2) - DV3(2) * DV1(1)
!
! Normalize Uij to unit vectors.
!
        S12 = 0.D0
        S23 = 0.D0
        S31 = 0.D0
        DO 2 I = 1, 3
            S12 = S12 + U12(I) * U12(I)
            S23 = S23 + U23(I) * U23(I)
            S31 = S31 + U31(I) * U31(I)
2           CONTINUE
!
! Test for a degenerate triangle associated with collinear
!   vertices.
!
            IF (S12 .EQ. 0.D0 .OR. S23 .EQ. 0.D0 .OR. S31 .EQ. 0.D0) THEN
                AREAS = 0.
                RETURN
            END IF
            S12 = SQRT(S12)
            S23 = SQRT(S23)
            S31 = SQRT(S31)
            DO 3 I = 1, 3
                U12(I) = U12(I) / S12
                U23(I) = U23(I) / S23
                U31(I) = U31(I) / S31
3               CONTINUE
!
! Compute interior angles Ai as the dihedral angles between
!   planes:
!           CA1 = cos(A1) = -<U12,U31>
!           CA2 = cos(A2) = -<U23,U12>
!           CA3 = cos(A3) = -<U31,U23>
!
                CA1 = -U12(1) * U31(1) - U12(2) * U31(2) - U12(3) * U31(3)
                CA2 = -U23(1) * U12(1) - U23(2) * U12(2) - U23(3) * U12(3)
                CA3 = -U31(1) * U23(1) - U31(2) * U23(2) - U31(3) * U23(3)
                IF (CA1 .LT. -1.D0) CA1 = -1.D0
                IF (CA1 .GT. 1.D0) CA1 = 1.D0
                IF (CA2 .LT. -1.D0) CA2 = -1.D0
                IF (CA2 .GT. 1.D0) CA2 = 1.D0
                IF (CA3 .LT. -1.D0) CA3 = -1.D0
                IF (CA3 .GT. 1.D0) CA3 = 1.D0
                A1 = ACOS(CA1)
                A2 = ACOS(CA2)
                A3 = ACOS(CA3)
!
! Compute AREAS = A1 + A2 + A3 - PI.
!
                AREAS = REAL(A1 + A2 + A3 - ACOS(-1.D0))
                IF (AREAS .LT. 0.) AREAS = 0.
                RETURN
            END
