SUBROUTINE INSERT(K, LP, LIST, LPTR, LNEW)
    INTEGER K, LP, LIST(*), LPTR(*), LNEW
!
!***********************************************************
!   This subroutine inserts K as a neighbor of N1 following
! N2, where LP is the LIST pointer of N2 as a neighbor of
! N1.  Note that, if N2 is the last neighbor of N1, K will
! become the first neighbor (even if N1 is a boundary node).
!
!   This routine is identical to the similarly named routine
! in TRIPACK.
!
!
! On input:
!
!       K = Index of the node to be inserted.
!
!       LP = LIST pointer of N2 as a neighbor of N1.
!
! The above parameters are not altered by this routine.
!
!       LIST,LPTR,LNEW = Data structure defining the trian-
!                        gulation.  Refer to Subroutine
!                        TRMESH.
!
! On output:
!
!       LIST,LPTR,LNEW = Data structure updated with the
!                        addition of node K.
!
! Modules required by INSERT:  None
!
!***********************************************************
!
    INTEGER LSAV
!
    LSAV = LPTR(LP)
    LPTR(LP) = LNEW
    LIST(LNEW) = K
    LPTR(LNEW) = LSAV
    LNEW = LNEW + 1
    RETURN
END
