SUBROUTINE GDR(N,NN1,KI,NN3,DENS,BOXL,CX,CY)
!*********************************************************************************************************************
!   1. AUTHORS AND E-MAIL CONTACT:
!   -Isaac Neri Gomez Sarmiento (isaacneri.gs@gmail.com)
!   -César Omar Ramírez Álvarez (cesaromarramirezalvarez@gmail.com)
!   -Jonas Valenzuela Terán     (24jonass@gmail.com)
!   -We thank the support and advice of our molecular simulations teacher Dr. Laura Lorenia Yeomans Reyna
!*********************************************************************************************************************
!   2. DESCRIPTION:
!       THIS SUBROUTINE ALLOWS US TO COMPUTE THE RADIAL DISTRIBUTION FUNCTION G(R).
!       G(R) GIVES US AN INSIGHT ON THE MEAN DISTANCE BETWEEN PARTICLES.
!       HIGH PEAKS MEAN A GREATER PROBABILITY TO FIND PARTICLES AROUND A REFERENCE PARTICLE IN THAT REGION.
!*********************************************************************************************************************

!********************
!VARIABLE DECLARATION
!********************
IMPLICIT NONE
INTEGER, INTENT(IN):: N !Number of particles
INTEGER, INTENT(IN):: NN1 !Maximum row size for particles' properties matrices: NN1=N+10
INTEGER, INTENT(IN):: KI !Counter for the number of saved simulation steps for position matrices with periodic boundary conditions
INTEGER, INTENT(IN):: NN3 !Number of histogram bins to compute g(r).
REAL, DIMENSION(NN1,KI), INTENT(IN):: CX, CY !Two dimensional arrays. Rows=particles, Columns:position at saved steps (periodic boundary conditions)
REAL, INTENT(IN):: BOXL !Size of the simulation square
REAL, INTENT(IN):: DENS ! Population density (reduced units)
REAL, INTENT(IN):: PI !Variable that holds the pi value
INTEGER :: i, j !Counters
INTEGER :: MAXBIN !Maximum number of bins for the histogram.
REAL :: DELTAR !We divide the space around a particle in MAXBIN layers of thickness DELTAR. From this we create a histogram of the number of particles in each layer averaged over the simulation time.
INTEGER :: l, m !Counters that iterate over all particles
INTEGER :: NBIN !Counter that iterates over each bin of the histogram.
REAL :: C1 !
REAL :: RL, RU !RL: Distance of a particle from its center to the lower part of a layer. !RU: Distance of a particle from its center to the upper part of a layer
REAL :: RT !Distance of a particle from its center to the half part of a layer
REAL :: GDRTA !Radial distribution function g(r)
REAL :: XL0, YL0 !These variables hold the value of the X, Y coordinates for a particle l
REAL :: XLT, YLT !These variables hold the value of the X, Y coordinates for a particle m
REAL :: XL0T, YL0T !Difference between X, Y coordinates of particle l and m.
REAL :: R0T !Distance between two particles, after imposing periodic boundary conditions.
INTEGER, DIMENSION(NN3) :: NHIST !Histogram


!Initializing histogram to 0's
DO i = 1, NN3
    NHIST(i) = 0
END DO


DELTAR = 0.01E0

MAXBIN = INT(BOXL/2.E0/DELTAR)

PI = 4.E0*ATAN(1.E0)

!Iterating over all particles and all simulation steps to obtain histogram.
DO l = 1,N
    DO m = 1,N
        IF (m == l) CYCLE !This is to avoid computing the distance of a particle with itself
            DO j = 1, KI

                XL0 = CX(l,j)
                XLT = CX(m,j)
                XL0T = XL0 - XLT

                YL0 = CY(l,j)
                YLT = CY(m,j)
                YL0T = YL0 - YLT


                XL0T = XL0T - BOXL*ANINT(XL0T/BOXL)
                YL0T = YL0T - BOXL*ANINT(YL0T/BOXL)

                R0T = SQRT(XL0T**2 + YL0T**2)

                NBIN = INT(R0T/DELTAR) + 1

                IF(NBIN <= MAXBIN) THEN

                    NHIST(NBIN) = NHIST(NBIN) + 1

                END IF
            END DO
    END DO
END DO


!Creating file where we'll store G(r)
OPEN(9, FILE='GDR.dat', STATUS='UNKNOWN')

!Computing G(R) using the histogram we created
DO NBIN = 1,MAXBIN

    RL = REAL(NBIN - 1)*DELTAR

    RU = RL + DELTAR

    RT = RL + DELTAR/2.E0

    C1 = PI*DENS*(RU**2 - RL**2)

    GDRTA = REAL(NHIST(NBIN))/REAL(KI)/REAL(N)/C1 !Computing G(r)

    !Saving r and G(r) in the created file
    WRITE(9,*) SNGL(RT),SNGL(GDRTA)

END DO
CLOSE(9)
END SUBROUTINE
