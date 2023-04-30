SUBROUTINE GETNP(X, Y, Z, LIST, LPTR, LEND, L, NPTS, DF, IER)
    INTEGER LIST(*), LPTR(*), LEND(*), L, NPTS(L), IER
    REAL X(*), Y(*), Z(*), DF
!
!***********************************************************
!   Given a Delaunay triangulation of N nodes on the unit
! sphere and an array NPTS containing the indexes of L-1
! nodes ordered by angular distance from NPTS(1), this sub-
! routine sets NPTS(L) to the index of the next node in the
! sequence -- the node, other than NPTS(1),...,NPTS(L-1),
! that is closest to NPTS(1).  Thus, the ordered sequence
! of K closest nodes to N1 (including N1) may be determined
! by K-1 calls to GETNP with NPTS(1) = N1 and L = 2,3,...,K
! for K .GE. 2.
!
!   The algorithm uses the property of a Delaunay triangula-
! tion that the K-th closest node to N1 is a neighbor of one
! of the K-1 closest nodes to N1.
!
!
! On input:
!
!       X,Y,Z = Arrays of length N containing the Cartesian
!               coordinates of the nodes.
!
!       LIST,LPTR,LEND = Triangulation data structure.  Re-
!                        fer to Subroutine TRMESH.
!
!       L = Number of nodes in the sequence on output.  2
!           .LE. L .LE. N.
!
! The above parameters are not altered by this routine.
!
!       NPTS = Array of length .GE. L containing the indexes
!              of the L-1 closest nodes to NPTS(1) in the
!              first L-1 locations.
!
! On output:
!
!       NPTS = Array updated with the index of the L-th
!              closest node to NPTS(1) in position L unless
!              IER = 1.
!
!       DF = Value of an increasing function (negative cos-
!            ine) of the angular distance between NPTS(1)
!            and NPTS(L) unless IER = 1.
!
!       IER = Error indicator:
!             IER = 0 if no errors were encountered.
!             IER = 1 if L < 2.
!
!***********************************************************
!
    INTEGER I, LM1, LP, LPL, N1, NB, NI, NP
    REAL DNB, DNP, X1, Y1, Z1
!
! Local parameters:
!
! DNB,DNP =  Negative cosines of the angular distances from
!              N1 to NB and to NP, respectively
! I =        NPTS index and DO-loop index
! LM1 =      L-1
! LP =       LIST pointer of a neighbor of NI
! LPL =      Pointer to the last neighbor of NI
! N1 =       NPTS(1)
! NB =       Neighbor of NI and candidate for NP
! NI =       NPTS(I)
! NP =       Candidate for NPTS(L)
! X1,Y1,Z1 = Coordinates of N1
!
    LM1 = L - 1
    IF (LM1 .LT. 1) GO TO 6
    IER = 0
!
! Store N1 = NPTS(1) and mark the elements of NPTS.
!
    N1 = NPTS(1)
    X1 = X(N1)
    Y1 = Y(N1)
    Z1 = Z(N1)
    DO 1 I = 1, LM1
        NI = NPTS(I)
        LEND(NI) = -LEND(NI)
1       CONTINUE
!
! Candidates for NP = NPTS(L) are the unmarked neighbors
!   of nodes in NPTS.  DNP is initially greater than -cos(PI)
!   (the maximum distance).
!
        DNP = 2.
!
! Loop on nodes NI in NPTS.
!
        DO 4 I = 1, LM1
            NI = NPTS(I)
            LPL = -LEND(NI)
            LP = LPL
!
! Loop on neighbors NB of NI.
!
2           NB = ABS(LIST(LP))
            IF (LEND(NB) .LT. 0) GO TO 3
!
! NB is an unmarked neighbor of NI.  Replace NP if NB is
!   closer to N1.
!
            DNB = -(X(NB) * X1 + Y(NB) * Y1 + Z(NB) * Z1)
            IF (DNB .GE. DNP) GO TO 3
            NP = NB
            DNP = DNB
3           LP = LPTR(LP)
            IF (LP .NE. LPL) GO TO 2
4           CONTINUE
            NPTS(L) = NP
            DF = DNP
!
! Unmark the elements of NPTS.
!
            DO 5 I = 1, LM1
                NI = NPTS(I)
                LEND(NI) = -LEND(NI)
5               CONTINUE
                RETURN
!
! L is outside its valid range.
!
6               IER = 1
                RETURN
            END
