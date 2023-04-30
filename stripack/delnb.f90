SUBROUTINE DELNB(N0, NB, N, LIST, LPTR, LEND, LNEW, LPH)
    INTEGER N0, NB, N, LIST(*), LPTR(*), LEND(N), LNEW, LPH
!
!***********************************************************
!   This subroutine deletes a neighbor NB from the adjacency
! list of node N0 (but N0 is not deleted from the adjacency
! list of NB) and, if NB is a boundary node, makes N0 a
! boundary node.  For pointer (LIST index) LPH to NB as a
! neighbor of N0, the empty LIST,LPTR location LPH is filled
! in with the values at LNEW-1, pointer LNEW-1 (in LPTR and
! possibly in LEND) is changed to LPH, and LNEW is decremen-
! ted.  This requires a search of LEND and LPTR entailing an
! expected operation count of O(N).
!
!   This routine is identical to the similarly named routine
! in TRIPACK.
!
!
! On input:
!
!       N0,NB = Indexes, in the range 1 to N, of a pair of
!               nodes such that NB is a neighbor of N0.
!               (N0 need not be a neighbor of NB.)
!
!       N = Number of nodes in the triangulation.  N .GE. 3.
!
! The above parameters are not altered by this routine.
!
!       LIST,LPTR,LEND,LNEW = Data structure defining the
!                             triangulation.
!
! On output:
!
!       LIST,LPTR,LEND,LNEW = Data structure updated with
!                             the removal of NB from the ad-
!                             jacency list of N0 unless
!                             LPH < 0.
!
!       LPH = List pointer to the hole (NB as a neighbor of
!             N0) filled in by the values at LNEW-1 or error
!             indicator:
!             LPH > 0 if no errors were encountered.
!             LPH = -1 if N0, NB, or N is outside its valid
!                      range.
!             LPH = -2 if NB is not a neighbor of N0.
!
!***********************************************************
!
    INTEGER I, LNW, LP, LPB, LPL, LPP, NN
!
! Local parameters:
!
! I =   DO-loop index
! LNW = LNEW-1 (output value of LNEW)
! LP =  LIST pointer of the last neighbor of NB
! LPB = Pointer to NB as a neighbor of N0
! LPL = Pointer to the last neighbor of N0
! LPP = Pointer to the neighbor of N0 that precedes NB
! NN =  Local copy of N
!
    NN = N
!
! Test for error 1.
!
    IF (N0 .LT. 1 .OR. N0 .GT. NN .OR. NB .LT. 1 .OR. NB .GT. NN .OR. NN .LT. 3) THEN
        LPH = -1
        RETURN
    END IF
!
!   Find pointers to neighbors of N0:
!
!     LPL points to the last neighbor,
!     LPP points to the neighbor NP preceding NB, and
!     LPB points to NB.
!
    LPL = LEND(N0)
    LPP = LPL
    LPB = LPTR(LPP)
1   IF (LIST(LPB) .EQ. NB) GO TO 2
    LPP = LPB
    LPB = LPTR(LPP)
    IF (LPB .NE. LPL) GO TO 1
!
!   Test for error 2 (NB not found).
!
    IF (ABS(LIST(LPB)) .NE. NB) THEN
        LPH = -2
        RETURN
    END IF
!
!   NB is the last neighbor of N0.  Make NP the new last
!     neighbor and, if NB is a boundary node, then make N0
!     a boundary node.
!
    LEND(N0) = LPP
    LP = LEND(NB)
    IF (LIST(LP) .LT. 0) LIST(LPP) = -LIST(LPP)
    GO TO 3
!
!   NB is not the last neighbor of N0.  If NB is a boundary
!     node and N0 is not, then make N0 a boundary node with
!     last neighbor NP.
!
2   LP = LEND(NB)
    IF (LIST(LP) .LT. 0 .AND. LIST(LPL) .GT. 0) THEN
        LEND(N0) = LPP
        LIST(LPP) = -LIST(LPP)
    END IF
!
!   Update LPTR so that the neighbor following NB now fol-
!     lows NP, and fill in the hole at location LPB.
!
3   LPTR(LPP) = LPTR(LPB)
    LNW = LNEW - 1
    LIST(LPB) = LIST(LNW)
    LPTR(LPB) = LPTR(LNW)
    DO 4 I = NN, 1, -1
    IF (LEND(I) .EQ. LNW) THEN
        LEND(I) = LPB
        GO TO 5
    END IF
4   CONTINUE
!
5   DO 6 I = 1, LNW - 1
        IF (LPTR(I) .EQ. LNW) THEN
            LPTR(I) = LPB
        END IF
6       CONTINUE
!
! No errors encountered.
!
        LNEW = LNW
        LPH = LPB
        RETURN
    END
