INTEGER FUNCTION NBCNT(LPL, LPTR)
    INTEGER LPL, LPTR(*)
!
!***********************************************************
!   This function returns the number of neighbors of a node
! N0 in a triangulation created by Subroutine TRMESH.
!
!   This function is identical to the similarly named
! function in TRIPACK.
!
!
! On input:
!
!       LPL = LIST pointer to the last neighbor of N0 --
!             LPL = LEND(N0).
!
!       LPTR = Array of pointers associated with LIST.
!
! Input parameters are not altered by this function.
!
! On output:
!
!       NBCNT = Number of neighbors of N0.
!
! Modules required by NBCNT:  None
!
!***********************************************************
!
    INTEGER K, LP
!
! Local parameters:
!
! K =  Counter for computing the number of neighbors
! LP = LIST pointer
!
    LP = LPL
    K = 1
!
1   LP = LPTR(LP)
    IF (LP .EQ. LPL) GO TO 2
    K = K + 1
    GO TO 1
!
2   NBCNT = K
    RETURN
END
