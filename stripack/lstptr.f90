INTEGER FUNCTION LSTPTR(LPL, NB, LIST, LPTR)
    INTEGER LPL, NB, LIST(*), LPTR(*)
!
!***********************************************************
!   This function returns the index (LIST pointer) of NB in
! the adjacency list for N0, where LPL = LEND(N0).
!
!   This function is identical to the similarly named
! function in TRIPACK.
!
!
! On input:
!
!       LPL = LEND(N0)
!
!       NB = Index of the node whose pointer is to be re-
!            turned.  NB must be connected to N0.
!
!       LIST,LPTR = Data structure defining the triangula-
!                   tion.  Refer to Subroutine TRMESH.
!
! Input parameters are not altered by this function.
!
! On output:
!
!       LSTPTR = Pointer such that LIST(LSTPTR) = NB or
!                LIST(LSTPTR) = -NB, unless NB is not a
!                neighbor of N0, in which case LSTPTR = LPL.
!
! Modules required by LSTPTR:  None
!
!***********************************************************
!
    INTEGER LP, ND
!
! Local parameters:
!
! LP = LIST pointer
! ND = Nodal index
!
    LP = LPTR(LPL)
1   ND = LIST(LP)
    IF (ND .EQ. NB) GO TO 2
    LP = LPTR(LP)
    IF (LP .NE. LPL) GO TO 1
!
2   LSTPTR = LP
    RETURN
END
