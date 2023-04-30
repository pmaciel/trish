SUBROUTINE COVSPH(KK, N0, LIST, LPTR, LEND, LNEW)
    INTEGER KK, N0, LIST(*), LPTR(*), LEND(*), LNEW
!
!***********************************************************
!   This subroutine connects an exterior node KK to all
! boundary nodes of a triangulation of KK-1 points on the
! unit sphere, producing a triangulation that covers the
! sphere.  The data structure is updated with the addition
! of node KK, but no optimization is performed.  All boun-
! dary nodes must be visible from node KK.
!
!
! On input:
!
!       KK = Index of the node to be connected to the set of
!            all boundary nodes.  KK .GE. 4.
!
!       N0 = Index of a boundary node (in the range 1 to
!            KK-1).  N0 may be determined by Subroutine
!            TRFIND.
!
! The above parameters are not altered by this routine.
!
!       LIST,LPTR,LEND,LNEW = Triangulation data structure
!                             created by Subroutine TRMESH.
!                             Node N0 must be included in
!                             the triangulation.
!
! On output:
!
!       LIST,LPTR,LEND,LNEW = Data structure updated with
!                             the addition of node KK as the
!                             last entry.  The updated
!                             triangulation contains no
!                             boundary nodes.
!
! Module required by COVSPH:  INSERT
!
!***********************************************************
!
    INTEGER K, LP, LSAV, NEXT, NST
!
! Local parameters:
!
! K =     Local copy of KK
! LP =    LIST pointer
! LSAV =  LIST pointer
! NEXT =  Boundary node visible from K
! NST =   Local copy of N0
!
    K = KK
    NST = N0
!
! Traverse the boundary in clockwise order, inserting K as
!   the first neighbor of each boundary node, and converting
!   the boundary node to an interior node.
!
    NEXT = NST
1   LP = LEND(NEXT)
    CALL INSERT(K, LP, LIST, LPTR, LNEW)
    NEXT = -LIST(LP)
    LIST(LP) = NEXT
    IF (NEXT .NE. NST) GO TO 1
!
! Traverse the boundary again, adding each node to K's
!   adjacency list.
!
    LSAV = LNEW
2   LP = LEND(NEXT)
    LIST(LNEW) = NEXT
    LPTR(LNEW) = LNEW + 1
    LNEW = LNEW + 1
    NEXT = LIST(LP)
    IF (NEXT .NE. NST) GO TO 2
!
    LPTR(LNEW - 1) = LSAV
    LEND(K) = LNEW - 1
    RETURN
END
