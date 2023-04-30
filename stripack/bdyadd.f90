SUBROUTINE BDYADD(KK, I1, I2, LIST, LPTR, LEND, LNEW)
    INTEGER KK, I1, I2, LIST(*), LPTR(*), LEND(*), LNEW
!
!***********************************************************
!   This subroutine adds a boundary node to a triangulation
! of a set of KK-1 points on the unit sphere.  The data
! structure is updated with the insertion of node KK, but no
! optimization is performed.
!
!   This routine is identical to the similarly named routine
! in TRIPACK.
!
!
! On input:
!
!       KK = Index of a node to be connected to the sequence
!            of all visible boundary nodes.  KK .GE. 1 and
!            KK must not be equal to I1 or I2.
!
!       I1 = First (rightmost as viewed from KK) boundary
!            node in the triangulation that is visible from
!            node KK (the line segment KK-I1 intersects no
!            arcs.
!
!       I2 = Last (leftmost) boundary node that is visible
!            from node KK.  I1 and I2 may be determined by
!            Subroutine TRFIND.
!
! The above parameters are not altered by this routine.
!
!       LIST,LPTR,LEND,LNEW = Triangulation data structure
!                             created by Subroutine TRMESH.
!                             Nodes I1 and I2 must be in-
!                             cluded in the triangulation.
!
! On output:
!
!       LIST,LPTR,LEND,LNEW = Data structure updated with
!                             the addition of node KK.  Node
!                             KK is connected to I1, I2, and
!                             all boundary nodes in between.
!
! Module required by BDYADD:  INSERT
!
!***********************************************************
!
    INTEGER K, LP, LSAV, N1, N2, NEXT, NSAV
!
! Local parameters:
!
! K =     Local copy of KK
! LP =    LIST pointer
! LSAV =  LIST pointer
! N1,N2 = Local copies of I1 and I2, respectively
! NEXT =  Boundary node visible from K
! NSAV =  Boundary node visible from K
!
    K = KK
    N1 = I1
    N2 = I2
!
! Add K as the last neighbor of N1.
!
    LP = LEND(N1)
    LSAV = LPTR(LP)
    LPTR(LP) = LNEW
    LIST(LNEW) = -K
    LPTR(LNEW) = LSAV
    LEND(N1) = LNEW
    LNEW = LNEW + 1
    NEXT = -LIST(LP)
    LIST(LP) = NEXT
    NSAV = NEXT
!
! Loop on the remaining boundary nodes between N1 and N2,
!   adding K as the first neighbor.
!
1   LP = LEND(NEXT)
    CALL INSERT(K, LP, LIST, LPTR, LNEW)
    IF (NEXT .EQ. N2) GO TO 2
    NEXT = -LIST(LP)
    LIST(LP) = NEXT
    GO TO 1
!
! Add the boundary nodes between N1 and N2 as neighbors
!   of node K.
!
2   LSAV = LNEW
    LIST(LNEW) = N1
    LPTR(LNEW) = LNEW + 1
    LNEW = LNEW + 1
    NEXT = NSAV
!
3   IF (NEXT .EQ. N2) GO TO 4
    LIST(LNEW) = NEXT
    LPTR(LNEW) = LNEW + 1
    LNEW = LNEW + 1
    LP = LEND(NEXT)
    NEXT = LIST(LP)
    GO TO 3
!
4   LIST(LNEW) = -N2
    LPTR(LNEW) = LSAV
    LEND(K) = LNEW
    LNEW = LNEW + 1
    RETURN
END
