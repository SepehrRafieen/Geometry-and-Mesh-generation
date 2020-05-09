  program mesh_generation

implicit none

integer::I,J,IMAX,JMAX,N,Z,kSI1,KSII,ETA1,ETAJ
REAL::DKSI,DETA,DEX,DEY,dx,dy,EX,EY,EPSILON,RESU,AI,CI &
      ,ASTARI,CSTARI,BI,DI,DTETA1,DTETA2,DR1,DR2,DR3,PIE,DTETA,ASTAR
real, allocatable ::X(:,:),Y(:,:),Xx(:,:),Yy(:,:),XOLD(:,:),YOLD(:,:),P(:,:)&
      ,Q(:,:),Alpha(:,:),Beta(:,:),GAMA(:,:),JACO(:,:)
real, allocatable ::AIX(:),BIX(:),CIX(:),DIX(:),AIY(:),BIY(:),CIY(:),DIY(:)&
      ,SIX(:),SIY(:),AJX(:),BJX(:),CJX(:),DJX(:),AJY(:),BJY(:),CJY(:),DJY(:) &
      ,SJX(:),SJY(:),TETA1(:),TETA2(:),TETA(:)
REAL ADX,ADY,ANGLE1,ANGLE2
real, allocatable ::ORTHOGONAL(:,:),SKEWNESS(:,:)
imax=25
jmax=25
DEX=1.0/(IMAX-1)
DEY=1.0/(JMAX-1)

allocate (X(IMAX,JMAX),Y(IMAX,JMAX),Xx(IMAX,JMAX),Yy(IMAX,JMAX) &
         ,XOLD(IMAX,JMAX),YOLD(IMAX,JMAX),P(IMAX,JMAX),Q(IMAX,JMAX) &
         ,Alpha(IMAX,JMAX),BETA(IMAX,JMAX),GAMA(IMAX,JMAX),JACO(IMAX,JMAX))

allocate(AIX(JMAX),BIX(JMAX),CIX(JMAX),DIX(JMAX),AIY(JMAX),BIY(JMAX) &
         ,CIY(JMAX),DIY(JMAX),SIY(JMAX),SIX(JMAX),AJX(Imax),BJX(Imax) &
         ,CJX(Imax),DJX(Imax),AJY(Imax),BJY(Imax),CJY(Imax),DJY(Imax) &
         ,SJY(Imax),SJX(Imax),TETA1(IMAX),TETA2(IMAX),TETA(IMAX))
ALLOCATE(ORTHOGONAL(IMAX,JMAX),SKEWNESS(IMAX,JMAX))
call initialANDBC
!CALL SQUARE
CALL MYGEOMETRY1_1
!CALL MYGEOMETRY1_2
!CALL MYGEOMETRY2_1
!CALL MYGEOMETRY2_2
call geometry1

do 

 XOLD(:,:)=X(:,:)
 YOLD(:,:)=Y(:,:)
call geometry2

call IXsweep
call IYsweep
call JXsweep
call JYsweep
call resual
print *, EPSILON

if(EPSILON<0.0001)EXIT

end do

call  ope
call outp
call CALSKEWNESS
100 continue

print*,"Done"    
pause
contains
!-----------------------------
subroutine initialANDBC
 x(:,:)=0.0
 y(:,:)=0.0
  DO I=1,IMAX
  DO J=1,JMAX
  X(I,J)=I
  Y(I,J)=J
  END DO 
  END DO
end subroutine  initialANDBC 
!--------------------------------
subroutine geometry1
P(:,:)=0.0
Q(:,:)=0.0
DKSI=1.0
DETA=1.0
KSI1=30
ETA1=44
AI=0.0
ASTAR=0.0
CI=0.0
BI=0.0
DI=0.0

DO J=1,JMAX
DO I=1,IMAX
KSII=I-KSI1
ETAJ=J-ETA1

IF(ETAJ>0)THEN
Q(I,J)=-ASTAR*EXP(-CI*ABS(J-ETA1))-BI*(+1)*EXP(-DI*((KSII**2+(J-ETA1)**2)**0.5))
END IF


IF(KSII>0)THEN
P(I,J)=-AI*EXP(-CI*ABS(I-KSI1))-BI*(+1)*EXP(-DI*((KSII**2+(J-ETA1)**2)**0.5))
END IF

IF (KSII==0)THEN
P(I,J)=0.0
END IF

IF (ETAJ==0)THEN
Q(I,J)=0.0
END IF


IF (KSII<0)THEN
P(I,J)=-AI*(-1)*EXP(-CI*ABS(I-KSI1))-BI*(-1)*EXP(-DI*((KSII**2+(J-ETA1)**2)**0.5))
END IF

IF (ETAJ<0)THEN
Q(I,J)=-ASTAR*(-1)*EXP(-CI*ABS(J-ETA1))-BI*(-1)*EXP(-DI*((KSII**2+(J-ETA1)**2)**0.5))
END IF

END DO
END DO

end subroutine geometry1
!----------------------------------
!---------------MYGEOMETRY1------------
subroutine MYGEOMETRY1_1
PIE=3.1415926535897932384626433832795
DX=20./(IMAX-1)
DY=10./(JMAX-1)
!LOWER WALL
DO I=1,IMAX
  X(I,1)=(I-1)*DX
  Y(I,1)=0.0
END DO
!RIGHT WALL
DO J=1,JMAX
  Y(IMAX,J)=(J-1)*DY
  X(IMAX,J)=20.0
END DO
!UPPER BOUNDARY
DX=15.0/((IMAX)-1)
DY=10.0/((IMAX+2)/3-1)
DO I=1,(IMAX+2)/3
  X(I,JMAX)=5.0+(I-1)*DX
  Y(I,JMAX)=20.0-(I-1)*DY
END DO
DO I=(IMAX+2)/3+1,IMAX
  X(I,JMAX)=X(I-1,JMAX)+DX
  Y(I,JMAX)=10.0
END DO
!LEFT BOUNDARY
DY=10.0/((JMAX+1)/2-1)
DO J=1,(JMAX+1)/2
  Y(1,J)=(J-1)*DY
  X(1,J)=0.0
END DO
DX=5.0/((JMAX+1)/2-1)
DO J=(JMAX+1)/2+1,JMAX
  Y(1,J)=Y(1,J-1)+DY
  X(1,J)=X(1,J-1)+DX
END DO
WRITE(*,*)Y(1,JMAX)
end subroutine MYGEOMETRY1_1
!--------------------------------------------
!--------------------MYGEOMETRY1-----------------
subroutine MYGEOMETRY1_2
PIE=3.1415926535897932384626433832795
DX=20.0/(IMAX-1)
DY=10.0/(JMAX-1)
!LOWER WALL
DO I=1,IMAX
  X(I,1)=(I-1)*DX
  Y(I,1)=0.0
END DO
!RIGHT WALL
DY=10.0/(JMAX-1)
DO J=1,JMAX
  Y(IMAX,J)=(J-1)*DY
  X(IMAX,J)=20.0
END DO
!UPPER BOUNDARY
DX=20.0/((IMAX)-1)
DY=10.0/((IMAX+3)/4-1)
DO I=1,(IMAX+3)/4
  X(I,JMAX)=(I-1)*DX
  Y(I,JMAX)=10.+(I-1)*DY
END DO
DO I=(IMAX+3)/4+1,(IMAX+3)/4*2-1
  X(I,JMAX)=(I-1)*DX
  Y(I,JMAX)=Y(I-1,JMAX)-DY
END DO
DO I=(IMAX+3)/4*2,IMAX
  X(I,JMAX)=(I-1)*DX
  Y(I,JMAX)=10.0
END DO
!LEFT BOUNDARY
DY=10./(JMAX-1)
DO J=1,JMAX
  Y(1,J)=(J-1)*DY
  X(1,J)=0.0
END DO
READ(*,*)
WRITE(*,*)Y(9,JMAX),Y(8,JMAX)
!$$$$$$ DX=5./((JMAX+1)/2-1)
!$$$$$$ DO J=(JMAX+1)/2+1,JMAX
!$$$$$$   Y(1,J)=Y(1,J-1)+DY
!$$$$$$   X(1,J)=X(1,J-1)+DX
!$$$$$$ END DO
end subroutine MYGEOMETRY1_2
!--------------------------------------------
!-------------------MYGEOMETRY2_1----------------
subroutine MYGEOMETRY2_1
!DOWN WALL
DX=10./(IMAX-1)
DO I=1,IMAX
  X(I,1)=(I-1)*DX
  Y(I,1)=0.0
END DO
!UP WALL
DY=10./(IMAX-1)
DO I=1,IMAX
  Y(I,JMAX)=20.-(I-1)*DY
  X(I,JMAX)=20.
END DO
!RIGHT WALL
DY=10./((JMAX+1)/2-1)
DX=10./((JMAX+1)/2-1)
DO J=1,((JMAX+1)/2)
  Y(IMAX,J)=(J-1)*DY
  X(IMAX,J)=10.0
END DO
DO J=((JMAX+1)/2)+1,JMAX
  Y(IMAX,J)=10.
  X(IMAX,J)=X(IMAX,J-1)+DX
END DO
!LEFT WALL  
DX=10./((JMAX+2)/3-1)
DY=10./((JMAX+2)/3-1)
DO J=1,(JMAX+2)/3
  X(1,J)=0.0
  Y(1,J)=(J-1)*DY
END DO
DO J=(JMAX+2)/3+1,((JMAX+2)/3)*2-1
  X(1,J)=X(1,J-1)+DX
  Y(1,J)=Y(1,J-1)+DY
END DO
DO J=((JMAX+2)/3)*2,JMAX
  X(1,J)=X(1,J-1)+DX
  Y(1,J)=20.
END DO
end subroutine MYGEOMETRY2_1
!--------------------------------------------
!---------------------MYGEOMETRY2_2--------------
subroutine MYGEOMETRY2_2
!DOWN WALL
DX=10./(IMAX-1)
DO I=1,IMAX
  X(I,1)=(I-1)*DX
  Y(I,1)=0.0
END DO
!UP WALL
DX=20./(IMAX-1)

DO I=1,(IMAX+1)/2
  Y(I,JMAX)=10.+(I-1)*DX
  X(I,JMAX)=(I-1)*DX
END DO
DO I=(IMAX+1)/2+1,IMAX
  Y(I,JMAX)=20.0
  X(I,JMAX)=(I-1)*DX
END DO
!RIGHT WALL
DY=10./((JMAX+2)/3-1)
DX=10./((JMAX+2)/3-1)
DO J=1,((JMAX+2)/3)
  Y(IMAX,J)=(J-1)*DY
  X(IMAX,J)=10.
END DO
DO J=((JMAX+2)/3)+1,((JMAX+2)/3)*2-1
  Y(IMAX,J)=10.
  X(IMAX,J)=X(IMAX,J-1)+DX
END DO
DO J=((JMAX+2)/3)*2,JMAX
  Y(IMAX,J)=Y(IMAX,J-1)+DY
  X(IMAX,J)=20.0
END DO
!LEFT WALL  
DY=10./(JMAX-1)
DO J=1,(JMAX)
  X(1,J)=0.0
  Y(1,J)=(J-1)*DY
END DO

end subroutine MYGEOMETRY2_2
!----------------------------------
subroutine geometry2
GAMA(:,:)=0.0
ALPHA(:,:)=0.0
BETA(:,:)=0.0
JACO(:,:)=0.0
DO J=2,JMAX-1
DO I=2,IMAX-1
ALPHA(I,J)=((X(I,J+1)-X(I,J-1))/(2.0*DETA))**2+((Y(I,J+1)-Y(I,J-1))/(2.0*DETA))**2
GAMA(I,J)=((X(I+1,J)-X(I-1,J))/(2.0*DKSI))**2+((Y(I+1,J)-Y(I-1,J))/(2.0*DKSI))**2
BETA(I,J)=((Y(I,J+1)-Y(I,J-1))/(2.0*DETA))*((Y(I+1,J)-Y(I-1,J))/(2.0*DKSI))&
           +((X(I,J+1)-X(I,J-1))/(2.0*DETA))*((X(I+1,J)-X(I-1,J))/(2.0*DKSI))
JACO(I,J)=((X(I+1,J)-X(I-1,J))/(2.0*DKSI))*((Y(I,J+1)-Y(I,J-1))/(2.0*DETA)) &
           -((X(I,J+1)-X(I,J-1))/(2.0*DETA))*((Y(I+1,J)-Y(I-1,J))/(2.0*DKSI))
END DO
END DO

end subroutine geometry2
!-----------------------------------
subroutine IXsweep

do i=2,imax-1
AIX(:)=0.0
BIX(:)=0.0
CIX(:)=0.0
DIX(:)=0.0
SIX(:)=0.0
do j=2,jmax-1
AIX(J)=GAMA(I,J)/DETA**2-Q(I,J)*JACO(I,J)**2/(2*DETA)
BIX(J)=-2.0*(ALPHA(I,J)/DKSI**2+GAMA(I,J)/DETA**2)
CIX(J)=GAMA(I,J)/DETA**2+Q(I,J)*JACO(I,J)**2/(2*DETA)
DIX(J)=BETA(I,J)*(X(I+1,J+1)-X(I-1,J+1)-X(I+1,J-1)+X(I-1,J-1))/(2.0*DETA*DKSI)&
       -ALPHA(I,J)*(X(I-1,J)+X(I+1,J))/DKSI**2-(P(I,J)*JACO(I,J)**2)*(X(I+1,J)-X(I-1,J))/(2.0*DKSI)


end do
 
AIX(1)=0.0
BIX(1)=1.0
CIX(1)=0.0
DIX(1)=X(I,1)

AIX(jmax)=0.0
BIX(jmax)=1.0
CIX(jmax)=0.0
DIX(jmax)=X(I,JMAX)


call TDMAIX
X(i,:)=SIX(:)

end do
end subroutine IXsweep
!---------------------------------
subroutine IYsweep

do i=2,imax-1
AIY(:)=0.0
BIY(:)=0.0
CIY(:)=0.0
DIY(:)=0.0
SIY(:)=0.0

do j=2,jmax-1
AIY(J)=GAMA(I,J)/DETA**2-Q(I,J)*JACO(I,J)**2/(2*DETA)
BIY(J)=-2.0*(ALPHA(I,J)/DKSI**2+GAMA(I,J)/DETA**2)
CIY(J)=GAMA(I,J)/DETA**2+Q(I,J)*JACO(I,J)**2/(2*DETA)
DIY(J)=BETA(I,J)*(Y(I+1,J+1)-Y(I-1,J+1)-Y(I+1,J-1)+Y(I-1,J-1))/(2.0*DETA*DKSI) &
       -ALPHA(I,J)*(Y(I-1,J)+Y(I+1,J))/DKSI**2-(P(I,J)*JACO(I,J)**2)*(Y(I+1,J)-Y(I-1,J))/(2.0*DKSI)


end do



AIY(1)=0.0
BIY(1)=1.0
CIY(1)=0.0
DIY(1)=Y(I,1)

AIY(jmax)=0.0
BIY(jmax)=1.0
CIY(jmax)=0.0
DIY(jmax)=Y(I,JMAX)


call TDMAIY
Y(i,:)=SIY(:)

end do
end subroutine IYsweep
!---------------------------------------
subroutine JXsweep

do J=2,Jmax-1
AJX(:)=0.0
BJX(:)=0.0
CJX(:)=0.0
DJX(:)=0.0
SJX(:)=0.0
do I=2,Imax-1
AJX(I)=ALPHA(I,J)/DKSI**2-P(I,J)*JACO(I,J)**2/(2.0*DKSI)
BJX(I)=-2.0*(ALPHA(I,J)/DKSI**2+GAMA(I,J)/DETA**2)
CJX(I)=ALPHA(I,J)/DKSI**2+P(I,J)*JACO(I,J)**2/(2.0*DKSI)
DJX(I)=BETA(I,J)*(X(I+1,J+1)-X(I-1,J+1)-X(I+1,J-1)+X(I-1,J-1))&
/(2.0*DETA*DKSI)-GAMA(I,J)*(X(I,J+1)+X(I,J-1))/DETA**2-(Q(I,J)*JACO(I,J)**2)*(X(I,J+1)-X(I,J-1))/(2.0*DETA)



end do

AJX(1)=0.0
BJX(1)=1.0
CJX(1)=0.0
DJX(1)=X(1,J)

AJX(Imax)=0.0
BJX(Imax)=1.0
CJX(Imax)=0.0
DJX(Imax)=X(IMAX,J)



call TDMAJX
X(:,J)=SJX(:)

end do
end subroutine JXsweep
!------------------------------------
subroutine JYsweep
do J=2,Jmax-1



AJY(:)=0.0
BJY(:)=0.0
CJY(:)=0.0
DJY(:)=0.0
SJY(:)=0.0


do I=2,Imax-1
AJY(I)=ALPHA(I,J)/DKSI**2-P(I,J)*JACO(I,J)**2/(2.0*DKSI)
BJY(I)=-2.0*(ALPHA(I,J)/DKSI**2+GAMA(I,J)/DETA**2)
CJY(I)=ALPHA(I,J)/DKSI**2+P(I,J)*JACO(I,J)**2/(2.0*DKSI)
DJY(I)=BETA(I,J)*(Y(I+1,J+1)-Y(I-1,J+1)-Y(I+1,J-1)+Y(I-1,J-1))&
/(2.0*DETA*DKSI)-GAMA(I,J)*(Y(I,J+1)+Y(I,J-1))/DETA**2-(Q(I,J)*JACO(I,J)**2)*(Y(I,J+1)-Y(I,J-1))/(2.0*DETA)

END DO


AJY(1)=0.0
BJY(1)=1.0
CJY(1)=0.0
DJY(1)=Y(1,J)

AJY(Imax)=0.0
BJY(Imax)=1.0
CJY(Imax)=0.0
DJY(Imax)=Y(IMAX,J)

call TDMAJY
Y(:,J)=SJY(:)
end do
end subroutine JYsweep

!----------------------------------

	subroutine tdmaIX
real ::f
n=0
do z=2,Jmax
f=AIX(z)/BIX(z-1)
BIX(z)=BIX(z)-CIX(z-1)*f
DIX(z)=DIX(z)-DIX(z-1)*f
end do
SIX(Jmax)=DIX(Jmax)/BIX(Jmax)
do n=Jmax-1,1,-1
SIX(n)=(DIX(n)-CIX(n)*SIX(n+1))/BIX(n)
end do
return
end subroutine tdmaIX

!-------------------------------------------
subroutine tdmaIY
real ::f
n=0
do z=2,Jmax
f=AIY(z)/BIY(z-1)
BIY(z)=BIY(z)-CIY(z-1)*f
DIY(z)=DIY(z)-DIY(z-1)*f
end do
SIY(Jmax)=DIY(Jmax)/BIY(Jmax)
do n=Jmax-1,1,-1
SIY(n)=(DIY(n)-CIY(n)*SIY(n+1))/BIY(n)
end do
return
end subroutine tdmaIY
!-------------------------------------------
	subroutine tdmaJX
real ::f
n=0
do z=2,Imax
f=AJX(z)/BJX(z-1)
BJX(z)=BJX(z)-CJX(z-1)*f
DJX(z)=DJX(z)-DJX(z-1)*f
end do
SJX(Imax)=DJX(Imax)/BJX(Imax)
do n=Imax-1,1,-1
SJX(n)=(DJX(n)-CJX(n)*SJX(n+1))/BJX(n)
end do
return
end subroutine tdmaJX

!------------------------------------------
subroutine tdmaJY
real ::f
n=0
do z=2,Imax
f=AJY(z)/BJY(z-1)
BJY(z)=BJY(z)-CJY(z-1)*f
DJY(z)=DJY(z)-DJY(z-1)*f
end do
SJY(Imax)=DJY(Imax)/BJY(Imax)
do n=Imax-1,1,-1
SJY(n)=(DJY(n)-CJY(n)*SJY(n+1))/BJY(n)
end do
return
end subroutine tdmaJY
!-----------------------------------------
subroutine resual

    Ex=0
	Ey=0
	do  i=1,imax
	do  j=1,jmax
	EX=max(Ex,abs(X(i,j)-XOLD(i,j)))
    EY=max(Ey,abs(Y(i,j)-YOLD(i,j)))
	
	END DO
END DO

	EPSILON=Ex+Ey
end subroutine resual
!-----------------------------------
subroutine ope
open(10, file='mesh.plt')
OPEN(20, FILE='ERR.PLT')
end subroutine ope
!---------------------------------
subroutine outp
write (10,*)  "ZONE ,j=",JMAX ,"I=",IMAX
do j=1,jmax
do i=1,imax
   
 write (10,*)    x(i,j),y(i,j)
  end do
end do

 
end subroutine outp
!--------------------------------------------
!-------------------------------------
subroutine CALSKEWNESS
write (20,*)  'VARIABLES="X","Y","ORTHOGONALITY","SKEWNESS"'
write (20,*)  "ZONE ,j=",JMAX ,"I=",IMAX,'F=POINT'
PIE=4.*ATAN(1.)
do j=2,jmax-1
do i=2,imax-1
   ADX= (X(I+1,J)-X(I,J))+(X(I+1,J+1)-X(I,J+1))
   ADY= (Y(I,J+1)-Y(I,J))+(Y(I+1,J+1)-Y(I+1,J))
   ANGLE1= ATAN2((Y(I,J+1)-Y(I,J-1)),(X(I,J+1)-X(I,J-1)))
   ANGLE2= ATAN2((Y(I+1,J)-Y(I-1,J)),(X(I+1,J)-X(I-1,J)))
   ORTHOGONAL(I,J)= ABS(PIE/2.-ABS(ANGLE1-ANGLE2))*180./PIE
   SKEWNESS(I,J)=ABS(ADY/ADX)
   

END DO
END DO
DO I=2,IMAX-1
ORTHOGONAL(I,JMAX)=ORTHOGONAL(I,JMAX-1)
ORTHOGONAL(I,1)=ORTHOGONAL(I,2)
SKEWNESS(I,JMAX)=SKEWNESS(I,JMAX-1)
SKEWNESS(I,1)=SKEWNESS(I,2)
END DO
DO J=1,JMAX
ORTHOGONAL(IMAX,J)=ORTHOGONAL(IMAX-1,J)
ORTHOGONAL(1,J)=ORTHOGONAL(2,J)
SKEWNESS(IMAX,J)=SKEWNESS(IMAX-1,J)
SKEWNESS(1,J)=SKEWNESS(2,J)
END DO  

do j=1,jmax
do i=1,imax
 write (20,*)    x(i,j),y(i,j),ORTHOGONAL(I,J),SKEWNESS(I,J)
 END DO
 END DO
 
end subroutine CALSKEWNESS
!----------------------------------	

    end program
    
    