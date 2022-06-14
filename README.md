# mtgso12

## Compare results of four implementations of Gram-Schmidt Orthogonalization

* Two fortran subroutines mtgso1 and mtgso2 implements modified
Gram-Schmidt orthogonalization, the second one is optimized for 2D
surveying networks and swap colums corresponding to point coordinates
x and y to improve numerical stability.

* C++ template class GSO from the project GNU Gama
https://www.gnu.org/software/gama/. Clone Gama project from git
server: ``git clone https://git.savannah.gnu.org/git/gama.git``

* C++ class ICGS from this directory implementing "Iterated Classical
Gram-Schmidt algorithm" (ICGS)

### Parameters of fortran subroutines

`````
      SUBROUTINE MTGSO1(A,ADIM,M1,N1,N2,N4,MODE,SC,EPS,TOL,
     /                  NX,DF,PROC,IER,INDX,INDS,CP)
      IMPLICIT INTEGER (A-Z)
      INTEGER  ADIM
      INTEGER  M1,N1,N2,N4,MODE,NX,DF,IER,INDX(1),INDS(1)
      LOGICAL  SC
      REAL     A(ADIM),EPS,TOL
      DOUBLE PRECISION CP(1)
      EXTERNAL PROC

Input/output paramaters: A, DF, INDS
Output parameter: IER
Input parameters: ADIM, M1, N1, N2, N4, MODE, SC, EPS, TOL NX, INDX
Other parameters: CP
``````
