! By default, a triangulation is created from a set of N nodes consisting of the north pole and N-1 points uniformly distributed
! around the 60-degree parallel (with constant longitudinal separation). However, by enabling the READ statements below (C# in the
! first two columns), testing may be performed on an arbitrary set of up to NMAX nodes (refer to the PARAMETER statement below). A
! data set consists of the following sequence of records:
!
! N = Number of nodes (3 to NMAX)
! (RLAT(I),RLON(I), I = 1,N) = Nodal coordinates in degrees latitude (-90 to 90) and degrees longitude (-180 to 180)

PROGRAM MAIN
    INTEGER IER, IFLAG, K, LWK, LNEW, N, NCOL, NMAX, NN, NROW, NT, NT6, NTMX
    REAL DEL, SC, PHI, THETA

    PARAMETER(NMAX=100, NTMX=2 * NMAX, NT6=6 * NMAX, LWK=2 * NMAX, NCOL=NMAX, NROW=9)
    INTEGER, PARAMETER :: LOUT=7

    ! Array storage for the triangulation, work space, and nodal coordinates
    INTEGER LIST(NT6), LPTR(NT6), LEND(NMAX), IWK(LWK)
    REAL DS(NMAX), RLAT(NMAX), RLON(NMAX), X(NMAX), Y(NMAX), Z(NMAX)

    ! Array storage for the triangle list
    INTEGER LTRI(NROW, NTMX)

    ! Generate the default set of nodes as latitudinal and longitudinal coordinates. DEL is the separation in degrees between the nodes on the 60-degree line of latitude.
    N = 9
    RLAT(1) = 90.
    RLON(1) = 0.
    DEL = 360./REAL(N - 1)
    DO K = 2, N
        RLAT(K) = 60.
        RLON(K) = REAL(K - 2) * DEL
    END DO

    IF (N .LT. 3 .OR. N .GT. NMAX) THEN
        PRINT *, 'N is outside its valid range:  N =' , N
        STOP
    END IF

    ! Set X,Y to the values of RLON,RLAT [radian]
    SC = ATAN(1.) / 45.
    DO K = 1, N
        PHI = SC * RLAT(K)
        THETA = SC * RLON(K)
        X(K) = COS(PHI) * COS(THETA)
        Y(K) = COS(PHI) * SIN(THETA)
        Z(K) = SIN(PHI)
    END DO

    ! Transform spherical coordinates X and Y to Cartesian coordinates (X,Y,Z) on the unit sphere
    !CALL TRANS(N, Y, X, X, Y, Z)

    ! Create the triangulation and test the error flag
    CALL TRMESH(N, X, Y, Z, LIST, LPTR, LEND, LNEW, IWK, IWK(N + 1), DS, IER)
    IF (IER .EQ. -2) THEN
        PRINT *, 'Error: the first three nodes are collinear'
        STOP
    ELSEIF (IER .GT. 0) THEN
        PRINT *, 'Error: duplicate nodes encountered'
        STOP
    END IF

    ! Test TRLIST and TRLPRT by creating and printing a triangle list.
    CALL TRLIST(N, LIST, LPTR, LEND, NROW, NT, LTRI, IER)
    CALL TRLPRT(N, RLON, RLAT, Z, IFLAG, NROW, NT, LTRI, LOUT)
END
