SUBROUTINE INTRSC(P1, P2, CN, P, IER)
    INTEGER IER
    REAL P1(3), P2(3), CN(3), P(3)
!
!***********************************************************
!   Given a great circle C and points P1 and P2 defining an
! arc A on the surface of the unit sphere, where A is the
! shorter of the two portions of the great circle C12 assoc-
! iated with P1 and P2, this subroutine returns the point
! of intersection P between C and C12 that is closer to A.
! Thus, if P1 and P2 lie in opposite hemispheres defined by
! C, P is the point of intersection of C with A.
!
!
! On input:
!
!       P1,P2 = Arrays of length 3 containing the Cartesian
!               coordinates of unit vectors.
!
!       CN = Array of length 3 containing the Cartesian
!            coordinates of a nonzero vector which defines C
!            as the intersection of the plane whose normal
!            is CN with the unit sphere.  Thus, if C is to
!            be the great circle defined by P and Q, CN
!            should be P X Q.
!
! The above parameters are not altered by this routine.
!
!       P = Array of length 3.
!
! On output:
!
!       P = Point of intersection defined above unless IER
!           .NE. 0, in which case P is not altered.
!
!       IER = Error indicator.
!             IER = 0 if no errors were encountered.
!             IER = 1 if <CN,P1> = <CN,P2>.  This occurs
!                     iff P1 = P2 or CN = 0 or there are
!                     two intersection points at the same
!                     distance from A.
!             IER = 2 if P2 = -P1 and the definition of A is
!                     therefore ambiguous.
!
!***********************************************************
!
    INTEGER I
    REAL D1, D2, PP(3), PPN, T
!
! Local parameters:
!
! D1 =  <CN,P1>
! D2 =  <CN,P2>
! I =   DO-loop index
! PP =  P1 + T*(P2-P1) = Parametric representation of the
!         line defined by P1 and P2
! PPN = Norm of PP
! T =   D1/(D1-D2) = Parameter value chosen so that PP lies
!         in the plane of C
!
    D1 = CN(1) * P1(1) + CN(2) * P1(2) + CN(3) * P1(3)
    D2 = CN(1) * P2(1) + CN(2) * P2(2) + CN(3) * P2(3)
!
    IF (D1 .EQ. D2) THEN
        IER = 1
        RETURN
    END IF
!
! Solve for T such that <PP,CN> = 0 and compute PP and PPN.
!
    T = D1 / (D1 - D2)
    PPN = 0.
    DO 1 I = 1, 3
        PP(I) = P1(I) + T * (P2(I) - P1(I))
        PPN = PPN + PP(I) * PP(I)
1       CONTINUE
!
! PPN = 0 iff PP = 0 iff P2 = -P1 (and T = .5).
!
        IF (PPN .EQ. 0.) THEN
            IER = 2
            RETURN
        END IF
        PPN = SQRT(PPN)
!
! Compute P = PP/PPN.
!
        DO 2 I = 1, 3
            P(I) = PP(I) / PPN
2           CONTINUE
            IER = 0
            RETURN
        END
