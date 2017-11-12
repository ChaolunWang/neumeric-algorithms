!File: lab2.f90
!Content: this file contain the code of lab2, which is about using lapack routine to solve linear system
!Author: Chaolun Wang
!Date: 01/15/2016


PROGRAM main
    IMPLICIT NONE
    !data declearation for question 1a
    INTEGER :: i, j, k
    INTEGER,PARAMETER :: n=100
    DOUBLE PRECISION :: A(n, n), x(n), b(n),  norm
    
    !data declearation for question 1b
    DOUBLE PRECISION :: A_21(2, 2),A_22(2, 2), b_2(2), c_2(2)
    
    !data declearation for question 2a
    INTEGER,PARAMETER :: m=100, p=2
    DOUBLE PRECISION :: AB(3*p+1, m), b_3(m), x_3(m)
    
    !data declearation for question 2b
    INTEGER,PARAMETER :: X4=5                                  !X4 can be defined as certain iteration of the linear system, here it is 5 iterations
    DOUBLE PRECISION :: AB_2(3*p+1, m),  b_4(m)
    
    !Question 1
    !part a.
    
    PRINT*, 'Question 1'
    PRINT*, 'Part a.'
    !initialization of A, x 
	DO i = 1, n
	    DO j = 1, n
	       IF((i.EQ.j ).OR. (j.EQ.n)) THEN
	           A(i, j)=1
	       ELSE IF(j.GT.i) THEN
	           A(i, j)=-10
	       ELSE
	           A(i, j)=0
	       END IF   
	    END DO
	    x(i)=DBLE(i)/DBLE(n)
	END DO    
    
    !calculate b by A*x
	DO i = 1, n
        b(i)=0
	    DO j = 1, n
            b(i)=b(i)+x(j)*A(i, j)   
	    END DO	   
	END DO     
    
    CALL solveByDgesv(A, b, n)                           !using the dgesv packages to solve the linear system
    
    DO i = 1, n
        b(i)=b(i)-x(i)   
	END DO	   
    
    CALL lInfiniteNorm(b, n, norm)                       ! call the routine to calculat the infinite norm of the error vector
    PRINT*, 'The l-infinite error norm is:', norm  
	      
    !Question 1
    !part b.
    PRINT*, 'Question 1'
    PRINT*, 'Part b.'
    !initialization
    
    DATA A_21 / 4.5, 1.6, 3.1, 1.1 /
    DATA A_22 / 4.5, 1.6, 3.1, 1.1 /
    DATA b_2 / 19.249, 6.843 /
    DATA c_2 / 19.25, 6.84 /

    CALL solveByDgesv(A_21, b_2, 2)                         !using the dgesv packages to solve the linear system                                
    CALL solveByDgesv(A_22, c_2, 2)    
    PRINT*, 'The solution to Ax=b is(in the form x1 x2):'  !prin out the result on screen
    PRINT*, b_2(1), b_2(2)
    PRINT*, 'The solution to Ay=c is(in the form y1 y2):' 
    PRINT*, c_2(1), c_2(2) 

    !Question 2
    !part a.
    PRINT*, 'Question 2'
    PRINT*, 'Part a.'    
    !initialization of AB , x_3, b_3
    
    DO i = 1, 3*p+1                                                     
	    DO j = 1, m
	       IF(i.EQ.p+1 .AND. j.GT.p) THEN
	           AB(i, j)=2
	       ELSE IF(i.EQ.p+2 .AND. j.GT.p-1) THEN
	           AB(i, j)=-1
	       ELSE IF(i.EQ.p+3 ) THEN
	           AB(i, j)=8
	       ELSE IF(i.EQ.p+4 .AND. j.LT.m-p+2) THEN
	           AB(i, j)=1
	       ELSE IF(i.EQ.p+5 .AND. j.LT.m-p+1) THEN
	           AB(i, j)=3
	       ELSE
	           AB(i, j)=0
	       END IF   
	    END DO
	END DO    
    
    
    b_3(1)=9
    b_3(2)=10
    b_3(m-1)=11
    b_3(m)=12   
    DO i=3, m-2
        b_3(i)=13
    END DO
    

    
    DO i=1, m
        x_3(i)=1
    END DO
    
    !call band matrix solver routine
    CALL solveBandMatrix(AB, b_3, m, p)
    !calculate the l2 norm of error vector
    DO i = 1, m
        b_3(i)=b_3(i)-x_3(i)   
	END DO	   
    
    CALL l2Norm(b_3, m, norm)                           !print out the l2 error norm
    PRINT*, 'The l2 error norm is:', norm      
 
    !Question 2
    !part b.
    PRINT*, 'Question 2'
    PRINT*, 'Part b.'
    !initialization of AB_2 , b_4 
    DO i = 1, 3*p+1
	    DO j = 1, m
	       IF(i.EQ.p+1 .AND. j.GT.p) THEN
	           AB_2(i, j)=2
	       ELSE IF(i.EQ.p+2 .AND. j.GT.p-1) THEN
	           AB_2(i, j)=-1
	       ELSE IF(i.EQ.p+3 ) THEN
	           AB_2(i, j)=8
	       ELSE IF(i.EQ.p+4 .AND. j.LT.m-p+2) THEN
	           AB_2(i, j)=1
	       ELSE IF(i.EQ.p+5 .AND. j.LT.m-p+1) THEN
	           AB_2(i, j)=3
	       ELSE
	           AB_2(i, j)=0
	       END IF   
	    END DO
	END DO    
    
    
    b_4(1)=9
    b_4(2)=10
    b_4(m-1)=11
    b_4(m)=12   
    DO i=3, m-2
        b_4(i)=13
    END DO    
    
    !call function to solve multiple system based on requirement
    CALL solveMultipleSystem(AB_2, b_4, m, p, X4)
    !the value in b_4 will be the bn, n is the number of X4
    PRINT*, 'System have been solved'
    DO i=1, m                                                       !print out the final result of bn, n=X4, which can be changed from 1 to 5                                                         
        PRINT*, b_4(i)
    END DO        
  
END PROGRAM main

!end of main function

SUBROUTINE solveByDgesv(A, b, n)                                    !routine to solve the system by dgesv routine
   INTEGER :: n, ipiv(n), info
   DOUBLE PRECISION :: A(n, n), b(n)
   CALL dgesv(n, 1, A, n, ipiv, b, n, info)
   CALL information(info)
END SUBROUTINE solveByDgesv

SUBROUTINE solveBandMatrix(AB, b_3, m, p)                           !routine used to solve the band matrix
    INTEGER :: m, p, ipiv(m), info
    DOUBLE PRECISION :: AB(3*p+1, m), b_3(m)
    CALL dgbtrf(m, m, p, p, AB, 3*p+1, ipiv, info)                  !dgbtrf was used for factorization
    CALL information(info)
    CALL dgbtrs('N', m, p, p, 1,  AB, 3*p+1, ipiv, b_3, m, info)    !dgbtrs was used to solve the factorized system
    CALL information(info)
END SUBROUTINE solveBandMatrix



SUBROUTINE solveMultipleSystem(AB_2, b_4, m, p, k)                   !routine used to solve multiple systems
    INTEGER :: m, p, k, i, j, ipiv(m), info
    DOUBLE PRECISION :: AB_2(3*p+1, m), b_4(m), x_4(m)
    DO i=1, m
        x_4(i)=b_4(i)
    END DO 
    CALL dgbtrf(m, m, p, p, AB_2, 3*p+1, ipiv, info)                 !dgbtrf was used for factorization
    DO i=1 , k-1
        CALL dgbtrs('N', m, p, p, 1,  AB_2, 3*p+1, ipiv, b_4, m, info)   !dgbtrs was used in a loop inorder to solve a series of systems
        DO j = 1, m
           b_4(j)=b_4(j)+x_4(j)   
	    END DO 
    END DO
END SUBROUTINE solveMultipleSystem




SUBROUTINE lInfiniteNorm(b, n, norm)                                  !calculate the infinite norm
   INTEGER :: n
   DOUBLE PRECISION :: b(n), norm
   norm=0
   DO i= 1, n
      IF(ABS(b(i)) .GT. ABS(norm)) THEN
          norm=ABS(b(i))
      END IF
   END DO
END SUBROUTINE lInfiniteNorm

SUBROUTINE l2Norm(b, n, norm)                                         !calculate the l2 norm
   INTEGER :: n
   DOUBLE PRECISION :: b(n), norm
   norm=0
   DO i= 1, n
     norm=norm+b(i)*b(i)
   END DO
   norm=sqrt(norm)
END SUBROUTINE l2Norm

SUBROUTINE information(info)                                          !show the information (if the matrix solver working well)
    INTEGER :: info
    PRINT*, 'info :  ', info
    IF(info .EQ. 0) THEN
      PRINT*, 'matrix solving is successful'
    ELSE IF(info .LT. 0) THEN
      PRINT*, 'illegal argument'
    ELSE 
      PRINT*, 'factor U is singular'
    END IF
END SUBROUTINE information

!end of program
