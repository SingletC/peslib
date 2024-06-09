!===========================================================================
! *def DSYEVJ3(A, Q, W)
!  subroutine of diagonalizing a 3 by 3 matrix
!  A is the input matirx
!  Q is the eigenvector matrix
!  W is the eigenvalue
!===========================================================================
SUBROUTINE DSYEVJ3(A, Q, W)
      DOUBLE PRECISION A(3,3)
      DOUBLE PRECISION Q(3,3)
      DOUBLE PRECISION W(3)

      INTEGER          N
      PARAMETER        ( N = 3 )
    
      DOUBLE PRECISION SD, SO
      DOUBLE PRECISION S, C, T
      DOUBLE PRECISION G, H, Z, THETA
      DOUBLE PRECISION THRESH
      INTEGER          I, X, Y, R

      DO 10 X = 1, N
        Q(X,X) = 1.0D0
        DO 11, Y = 1, X-1
          Q(X, Y) = 0.0D0
          Q(Y, X) = 0.0D0
   11   CONTINUE
   10 CONTINUE

      DO 20 X = 1, N
        W(X) = A(X, X)
   20 CONTINUE

      SD = 0.0D0
      DO 30 X = 1, N
        SD = SD + ABS(W(X))
   30 CONTINUE
      SD = SD**2
 
      DO 40 I = 1, 50
        SO = 0.0D0
        DO 50 X = 1, N
          DO 51 Y = X+1, N
            SO = SO + ABS(A(X, Y))
   51     CONTINUE
   50   CONTINUE
        IF (SO .EQ. 0.0D0) THEN
          RETURN
        END IF

        IF (I .LT. 4) THEN
          THRESH = 0.2D0 * SO / N**2
        ELSE
          THRESH = 0.0D0
        END IF

        DO 60 X = 1, N
          DO 61 Y = X+1, N
            G = 100.0D0 * ( ABS(A(X, Y)) )
            IF ( I .GT. 4 .AND. ABS(W(X)) + G .EQ. ABS(W(X)) .AND. ABS(W(Y)) + G .EQ. ABS(W(Y)) ) THEN
              A(X, Y) = 0.0D0
            ELSE IF (ABS(A(X, Y)) .GT. THRESH) THEN
              H = W(Y) - W(X)
              IF ( ABS(H) + G .EQ. ABS(H) ) THEN
                T = A(X, Y) / H
              ELSE
                THETA = 0.5D0 * H / A(X, Y)
                IF (THETA .LT. 0.0D0) THEN
                  T = -1.0D0 / (SQRT(1.0D0 + THETA**2) - THETA)
                ELSE
                  T = 1.0D0 / (SQRT(1.0D0 + THETA**2) + THETA)
                END IF
              END IF

              C = 1.0D0 / SQRT( 1.0D0 + T**2 )
              S = T * C
              Z = T * A(X, Y)
              
              A(X, Y) = 0.0D0
              W(X)    = W(X) - Z
              W(Y)    = W(Y) + Z
              DO 70 R = 1, X-1
                T       = A(R, X)
                A(R, X) = C * T - S * A(R, Y)
                A(R, Y) = S * T + C * A(R, Y)
   70         CONTINUE
              DO 80, R = X+1, Y-1
                T       = A(X, R)
                A(X, R) = C * T - S * A(R, Y)
                A(R, Y) = S * T + C * A(R, Y)
   80         CONTINUE
              DO 90, R = Y+1, N
                T       = A(X, R)
                A(X, R) = C * T - S * A(Y, R)
                A(Y, R) = S * T + C * A(Y, R)
   90         CONTINUE

              DO 100, R = 1, N
                T       = Q(R, X)
                Q(R, X) = C * T - S * Q(R, Y)
                Q(R, Y) = S * T + C * Q(R, Y)
  100         CONTINUE
            END IF
   61     CONTINUE
   60   CONTINUE
   40 CONTINUE

      PRINT *, "DSYEVJ3: No convergence."
            
END SUBROUTINE

!===========================================================================
! *def order(a,rv)
!  subroutine of ordering three eigenstates
!  a is the eigenvalue
!  rv is the eigenvector
!===========================================================================
subroutine order(a,rv)
double precision a(3)
double precision rv(3,3)
double precision tmpv(3,3)
double precision s0, s1, s2
integer i, j

do i=1,3
do j=1,3
tmpv(i,j)=rv(i,j)
enddo
enddo

    if (a(1) < a(2)) then
      if (a(1) < a(3)) then
        if (a(2) < a(3)) then
          s0=a(1)
          s1=a(2)
          s2=a(3)
          do i=1,3
            rv(i,1)=tmpv(i,1)
          enddo
          do i=1,3
            rv(i,2)=tmpv(i,2)
          enddo
          do i=1,3
            rv(i,3)=tmpv(i,3)
          enddo
        else
          s0=a(1)
          s1=a(3)
          s2=a(2)
          do i=1,3
            rv(i,1)=tmpv(i,1)
          enddo
          do i=1,3
            rv(i,2)=tmpv(i,3)
          enddo
          do i=1,3
            rv(i,3)=tmpv(i,2)
          enddo
        end if
      else
          s0=a(3)
          s1=a(1)
          s2=a(2)
          do i=1,3
            rv(i,1)=tmpv(i,3)
          enddo
          do i=1,3
            rv(i,2)=tmpv(i,1)
          enddo
          do i=1,3
            rv(i,3)=tmpv(i,2)
          enddo
      end if
    else
      if (a(2) < a(3)) then
        if (a(1) < a(3)) then
          s0=a(2)
          s1=a(1)
          s2=a(3)
          do i=1,3
            rv(i,1)=tmpv(i,2)
          enddo
          do i=1,3
            rv(i,2)=tmpv(i,1)
          enddo
          do i=1,3
            rv(i,3)=tmpv(i,3)
          enddo
        else
          s0=a(2)
          s1=a(3)
          s2=a(1)
          do i=1,3
            rv(i,1)=tmpv(i,2)
          enddo
          do i=1,3
            rv(i,2)=tmpv(i,3)
          enddo
          do i=1,3
            rv(i,3)=tmpv(i,1)
          enddo
        end if
      else
          s0=a(3)
          s1=a(2)
          s2=a(1)
          do i=1,3
            rv(i,1)=tmpv(i,3)
          enddo
          do i=1,3
            rv(i,2)=tmpv(i,2)
          enddo
          do i=1,3
            rv(i,3)=tmpv(i,1)
          enddo
      end if
    end if

a(1)=s0
a(2)=s1
a(3)=s2
end subroutine
