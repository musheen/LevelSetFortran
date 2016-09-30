MODULE set_subs

!*************************************************************************************!
!
! Module containing subroutines for set3d.f90
!
!*************************************************************************************!

IMPLICIT NONE

CONTAINS

!*************************************************************************************!
! Read STL and Allocate
!*************************************************************************************!

SUBROUTINE stlRead(surfX,nSurfNode,surfElem,filename,nSurfElem,surfElemTag,surfOrder,nBndComp,nBndElem,bndNormal)

CHARACTER header*80
CHARACTER, INTENT(IN) :: filename*80
INTEGER*2 padding
INTEGER*4 ntri,iunit,nSurfNode,k,i,n,p,kk,share,nSurfElem
REAL*4,ALLOCATABLE,DIMENSION(:,:) :: normals,triangles,nodesT
INTEGER*4,ALLOCATABLE,DIMENSION(:,:),INTENT(OUT) :: surfElem
REAL,ALLOCATABLE,DIMENSION(:,:),INTENT(OUT) :: surfX
INTEGER,ALLOCATABLE,DIMENSION(:),INTENT(OUT) :: surfElemTag,surfOrder
REAL,ALLOCATABLE,DIMENSION(:,:),INTENT(OUT) :: bndNormal
INTEGER,INTENT(OUT) :: nBndComp,nBndElem

PRINT*,
PRINT*, " Reading in .stl Mesh "
PRINT*,

iunit=13
OPEN(unit=iunit,file=filename,status='old',access='stream',form='unformatted')

! read .stl header info 
READ(iunit) header
READ(iunit) ntri
   
ALLOCATE(normals(3,ntri))
ALLOCATE(triangles(3,ntri*3))
ALLOCATE(surfELem(ntri,3))
 
! read .stl data
k=1
DO i = 1,ntri
   READ(iunit) normals(1,i),normals(2,i),normals(3,i)
   READ(iunit) triangles(1,k),triangles(2,k),triangles(3,k)
   READ(iunit) triangles(1,k+1),triangles(2,k+1),triangles(3,k+1)
   READ(iunit) triangles(1,k+2),triangles(2,k+2),triangles(3,k+2)
   READ(iunit) padding
  k=k+3
END DO
  
CLOSE(iunit)

ALLOCATE(nodesT(3,ntri*5))
nSurfElem = ntri

! search through data and put into surfX and surfElem style arrays
DO k = 1,ntri
  nodesT(1,k) = 1000000. 
  nodesT(2,k) = 1000000. 
  nodesT(3,k) = 1000000. 
END DO

! eliminate repeated nodes and clean up arrays
i = 1
nSurfNode = 3
k = 0;
DO n = 1,ntri
   DO p = 1,3
      share = 0 
      DO kk = 1,nSurfNode
         IF ((abs(nodesT(1,kk) - triangles(1,i)) < 1.e-13) .AND. &
             (abs(nodesT(2,kk) - triangles(2,i)) < 1.e-13) .AND. &
             (abs(nodesT(3,kk) - triangles(3,i)) < 1.e-13)) THEN
            share = kk
            EXIT
         END IF
      END DO
      IF (share > 0) THEN
         surfElem(n,p) = share
      ELSE
         k             = k+1 
         nodesT(:,k)   = triangles(:,i)
         surfElem(n,p) = k !1-based
      END IF
      i = i+1
   END DO
   nSurfNode = k 
END DO

! allocate surfX
ALLOCATE(surfX(nSurfNode,3))

! fill in surface node data 
DO k = 1,nSurfNode
   surfX(k,1) = nodesT(1,k)
   surfX(k,2) = nodesT(2,k)
   surfX(k,3) = nodesT(3,k)
END DO

! deallocate unnecessary data
DEALLOCATE(nodesT)
DEALLOCATE(triangles)
DEALLOCATE(normals)

ALLOCATE(surfOrder(nSurfElem))
ALLOCATE(surfElemTag(nSurfElem))
ALLOCATE(bndNormal(nBndComp,3))

surfOrder = 1
surfElemTag = 0
nBndElem = 0
nBndComp = 1
bndNormal = 0.


ENDSUBROUTINE stlRead

!*************************************************************************************!
! Read S3D and Allocate
!*************************************************************************************!

SUBROUTINE s3dRead(surfX,nSurfNode,surfElem,filename,nSurfElem,surfElemTag,surfOrder,nBndComp,nBndElem,bndNormal)

CHARACTER, INTENT(IN) :: filename*80
INTEGER*4 ntri,iunit,nSurfNode,k,i,n,p,kk,share,nSurfElem
REAL*4,ALLOCATABLE,DIMENSION(:,:) :: normals,triangles,nodesT
INTEGER*4,ALLOCATABLE,DIMENSION(:,:),INTENT(OUT) :: surfElem
REAL,ALLOCATABLE,DIMENSION(:,:),INTENT(OUT) :: surfX
INTEGER,ALLOCATABLE,DIMENSION(:),INTENT(OUT) :: surfElemTag,surfOrder
REAL,ALLOCATABLE,DIMENSION(:,:),INTENT(OUT) :: bndNormal
INTEGER,INTENT(OUT) :: nBndComp,nBndElem

PRINT*,
PRINT*, " Reading in .s3d Mesh "
PRINT*,


PRINT*, " Working on this ... "
STOP


END SUBROUTINE s3dRead
!*************************************************************************************!
! Determine Sign
!*************************************************************************************!

SUBROUTINE phiSign(pS,sgn,dxx,gM)

REAL,INTENT(IN) :: pS,dxx,gM
REAL,INTENT(OUT) :: sgn
REAL :: eps = 1.E-13


! non smeared
!IF (pS < 0.) THEN
!   sgn =-1.
!ELSE IF (pS > 0.) THEN
!   sgn = 1.
!ELSE
!   sgn = 0.
!END IF

! smeared
sgn = pS/sqrt(pS*pS + dxx*dxx*gM)


ENDSUBROUTINE phiSign

!*************************************************************************************!
! Determine Narrow Band
!*************************************************************************************!

SUBROUTINE narrowBand(nx,ny,nz,dx,phi,phiNB,phiSB)

INTEGER,INTENT(IN) :: nx,ny,nz
INTEGER :: i,j,k
REAL,INTENT(IN) :: dx
REAL,DIMENSION(0:nx,0:ny,0:nz),INTENT(IN) :: phi
INTEGER,DIMENSION(0:nx,0:ny,0:nz),INTENT(INOUT) :: phiNB,phiSB

phiNB = 0
phiSB = 0

DO i = 0,nx
   DO j = 0,ny
      DO k = 0,nz
         
         ! narrow band
         IF (abs(phi(i,j,k)) < 4.1*dx) THEN
            phiNB(i,j,k) = 1
         END IF

         ! stencil band
         IF (abs(phi(i,j,k)) < 8.1*dx) THEN
            phiSB(i,j,k) = 1
         END IF

      END DO
   END DO
END DO

END SUBROUTINE narrowBand

!*************************************************************************************!
! Calculate First Derivative
!*************************************************************************************!

SUBROUTINE firstDeriv(i,j,k,nx,ny,nz,dx,phi,phiX,phiY,phiZ,order,gMM,gradPhi)

INTEGER,INTENT(IN) :: i,j,k,order,nx,ny,nz
REAL,INTENT(IN) :: dx
REAL,DIMENSION(0:nx,0:ny,0:nz),INTENT(IN) :: phi
REAL,DIMENSION(0:nx,0:ny,0:nz,3),INTENT(INOUT) :: gradPhi
REAL,INTENT(OUT) :: phiX,phiY,phiZ,gMM
REAL :: aa1,aa2,aa3,aa4,aa5,aa6,aa7,dxx,aa8,aa9
INTEGER :: im,jm,km,ip,jp,kp,im2,ip2,jm2,jp2,km2,kp2
INTEGER :: ip1,jp1,kp1,im1,jm1,km1,ip3,jp3,kp3,im3,jm3,km3,im4,jm4,km4,ip4,jp4,kp4

IF (order == 1) THEN

   IF (phi(i+1,j,k) > phi(i,j,k)) THEN
      phiX = (phi(i,j,k)+phi(i+1,j,k))/dx
   ELSE
      phiX = (phi(i,j,k)-phi(i-1,j,k))/dx
   END IF

   IF (phi(i,j+1,k) > phi(i,j,k)) THEN
      phiY = (phi(i,j,k)+phi(i,j+1,k))/dx
   ELSE
      phiY = (phi(i,j,k)-phi(i,j-1,k))/dx
   END IF

   IF (phi(i,j,k+1) > phi(i,j,k)) THEN
      phiZ = (phi(i,j,k)+phi(i,j,k+1))/dx
   ELSE
      phiZ = (phi(i,j,k)-phi(i,j,k-1))/dx
   END IF


ELSEIF (order == 2) THEN
  
   dxx = 1./(dx); 

   ! second order first derivative
   phiX = (-1./2.*phi(i-1,j,k) +  1./2.*phi(i+1,j,k))*dxx
   phiY = (-1./2.*phi(i,j-1,k) +  1./2.*phi(i,j+1,k))*dxx
   phiZ = (-1./2.*phi(i,j,k-1) +  1./2.*phi(i,j,k+1))*dxx

ELSEIF (order == 4) THEN 

    dxx = 1./(12.*dx)

   ! fourth order first derivatives
   phiX = (-phi(i+2,j,k)+8.*phi(i+1,j,k)-8.*phi(i-1,j,k)+phi(i-2,j,k))*dxx
   phiY = (-phi(i,j+2,k)+8.*phi(i,j+1,k)-8.*phi(i,j-1,k)+phi(i,j-2,k))*dxx
   phiZ = (-phi(i,j,k+2)+8.*phi(i,j,k+1)-8.*phi(i,j,k-1)+phi(i,j,k-2))*dxx

ELSEIF (order == 6) THEN

   im1 = i-1
   jm1 = j-1
   km1 = k-1
   im2 = i-2
   jm2 = j-2
   km2 = k-2
   im3 = i-3
   jm3 = j-3
   km3 = k-3

   ip1 = i+1
   jp1 = j+1
   kp1 = k+1
   ip2 = i+2
   jp2 = j+2
   kp2 = k+2
   ip3 = i+3
   jp3 = j+3
   kp3 = k+3

   aa1 =-1./60.;
   aa2 = 3./20.;
   aa3 =-3./4.;
   aa4 = 0.;
   aa5 = 3./4.;
   aa6 =-3./20.;
   aa7 = 1./60.;

   ! calculate derivatives
   phiX = (phi(im3,j,k)*aa1 + phi(im2,j,k)*aa2 + phi(im1,j,k)*aa3 + &
&         phi(i,j,k)*aa4 + phi(ip1,j,k)*aa5 + phi(ip2,j,k)*aa6 + phi(ip3,j,k)*aa7)/dx    

   phiY = (phi(i,jm3,k)*aa1 + phi(i,jm2,k)*aa2 + phi(i,jm1,k)*aa3 + &
&         phi(i,j,k)*aa4 + phi(i,jp1,k)*aa5 + phi(i,jp2,k)*aa6 + phi(i,jp3,k)*aa7)/dx   

   phiZ = (phi(i,j,km3)*aa1 + phi(i,j,km2)*aa2 + phi(i,j,km1)*aa3 + &
&         phi(i,j,k)*aa4 + phi(i,j,kp1)*aa5 + phi(i,j,kp2)*aa6 + phi(i,j,kp3)*aa7)/dx   

ELSEIF (order == 8) THEN

   im1 = i-1
   jm1 = j-1
   km1 = k-1
   im2 = i-2
   jm2 = j-2
   km2 = k-2
   im3 = i-3
   jm3 = j-3
   km3 = k-3
   im4 = i-4
   jm4 = j-4
   km4 = k-4

   ip1 = i+1
   jp1 = j+1
   kp1 = k+1
   ip2 = i+2
   jp2 = j+2
   kp2 = k+2
   ip3 = i+3
   jp3 = j+3
   kp3 = k+3
   ip4 = i+4
   jp4 = j+4
   kp4 = k+4


   aa1 = 1./280.;
   aa2 =-4./105.;
   aa3 = 1./5.;
   aa4 =-4./5.;
   aa5 = 0.;
   aa6 = 4./5;
   aa7 =-1./5.;
   aa8 = 4./105.
   aa9 =-1./280.

   ! calculate derivatives
   phiX = (phi(im4,j,k)*aa1 + phi(im3,j,k)*aa2 + phi(im2,j,k)*aa3 + phi(im1,j,k)*aa4 + &
&          phi(ip1,j,k)*aa6 + phi(ip2,j,k)*aa7 + phi(ip3,j,k)*aa8 + phi(ip4,j,k)*aa9)/dx    
   phiY = (phi(i,jm4,k)*aa1 + phi(i,jm3,k)*aa2 + phi(i,jm2,k)*aa3 + phi(i,jm1,k)*aa4 + &
&          phi(i,jp1,k)*aa6 + phi(i,jp1,k)*aa7 + phi(i,jp3,k)*aa8 + phi(i,jp4,k)*aa9)/dx   

   phiZ = (phi(i,j,km4)*aa1 + phi(i,j,km3)*aa2 + phi(i,j,km2)*aa3 + phi(i,j,km1)*aa4 + &
&          phi(i,j,kp1)*aa6 + phi(i,j,kp2)*aa7 + phi(i,j,kp3)*aa8 + phi(i,j,kp4)*aa9)/dx  

ELSE 
   PRINT*," This order not set yet "
   STOP
END IF

gradPhi(i,j,k,1) = phiX
gradPhi(i,j,k,2) = phiY
gradPhi(i,j,k,3) = phiZ

gMM = phiX*phiX + phiY*phiY + phiZ*phiZ

gMM = sqrt(gMM)

END SUBROUTINE firstDeriv

!*************************************************************************************!
! Calculate Second Derivative 
!*************************************************************************************!

SUBROUTINE secondDeriv(i,j,k,nx,ny,nz,dx,phi,phiXX,phiYY,phiZZ,phiXY,phiXZ,phiYZ,order)

REAL :: dxx,dx12
INTEGER,INTENT(IN) :: i,j,k,order,nx,ny,nz
REAL,INTENT(IN) :: dx
REAL,DIMENSION(0:nx,0:ny,0:nz),INTENT(IN) :: phi
REAL,INTENT(OUT) :: phiXX,phiYY,phiZZ,phiXY,phiXZ,phiYZ
REAL :: aa1,aa2,aa3,aa4,aa5,aa6,aa7
REAL :: bb1,bb2,bb3,bb4,bb5,bb6,bb7
INTEGER :: im,jm,km,ip,jp,kp,im2,ip2,jm2,jp2,km2,kp2
INTEGER :: ip1,jp1,kp1,im1,jm1,km1,ip3,jp3,kp3,im3,jm3,km3

IF (order == 2) THEN

   dxx = 1./(dx*dx) 

   ! second order second derivatives
   phiXX = (-2.*phi(i,j,k) + phi(i+1,j,k) + phi(i-1,j,k))*dxx;          
   phiYY = (-2.*phi(i,j,k) + phi(i,j+1,k) + phi(i,j-1,k))*dxx;
   phiZZ = (-2.*phi(i,j,k) + phi(i,j,k+1) + phi(i,j,k-1))*dxx;
       
   ! second order mixed derivatives
   phiXY = phi(i+1,j+1,k  )-phi(i+1,j-1,k  )-phi(i-1,j+1,k  )+phi(i-1,j-1,k  )
   phiYZ = phi(i  ,j+1,k+1)-phi(i  ,j+1,k-1)-phi(i  ,j-1,k+1)+phi(i  ,j-1,k-1)
   phiXZ = phi(i+1,j  ,k+1)-phi(i+1,j  ,k-1)-phi(i-1,j  ,k+1)+phi(i-1,j  ,k-1)

   phiXY = phiXY*dxx/4.
   phiYZ = phiYZ*dxx/4.
   phiXZ = phiXZ*dxx/4.

ELSE 

   PRINT*," This order not set yet "
   STOP

END IF   

END SUBROUTINE secondDeriv

!*************************************************************************************!
! Calculate Min/Max Flow
!*************************************************************************************!

SUBROUTINE minMax(i,j,k,nx,ny,nz,dx,phi,grad2Phi,gradMixPhi,F,gridX)

REAL :: dxx,curv,b1,b2,b3,b4,b5,b6,denom,gM,HH,pave,thresh
REAL :: phiXX,phiYY,phiZZ,phiXY,phiXZ,phiYZ,phiX,phiY,phiZ
INTEGER :: h,order2,order1
INTEGER,INTENT(IN) :: i,j,k,nx,ny,nz
REAL,INTENT(IN) :: dx
REAL,DIMENSION(0:nx,0:ny,0:nz),INTENT(INOUT) :: phi
REAL,DIMENSION(0:nx,0:ny,0:nz,3),INTENT(INOUT) :: grad2Phi,gradMixPhi,gridX
REAL,DIMENSION(0:nx,0:ny,0:nz,3) :: gradPhi
REAL,INTENT(OUT) :: F
      

! true curvature

!order1 = 2
!order2 = 2
 
!CALL firstDeriv(i,j,k,nx,ny,nz,dx,phi,phiX,phiY,phiZ,order1,gM,gradPhi)
!CALL secondDeriv(i,j,k,nx,ny,nz,dx,phi,phiXX,phiYY,phiZZ,phiXY,phiXZ,phiYZ,order2)

!b1 = grad2Phi(i,j,k,2)+grad2Phi(i,j,k,3) 
!b2 = grad2Phi(i,j,k,1)+grad2Phi(i,j,k,3) 
!b3 = grad2Phi(i,j,k,1)+grad2Phi(i,j,k,2)
!b4 = gradPhi(i,j,k,1)*gradPhi(i,j,k,2)*gradMixPhi(i,j,k,1);
!b5 = gradPhi(i,j,k,1)*gradPhi(i,j,k,3)*gradMixPhi(i,j,k,2);
!b6 = gradPhi(i,j,k,2)*gradPhi(i,j,k,3)*gradMixPhi(i,j,k,3);

!denom = gM*gM*gM

!curv = (b1*gradPhi(i,j,k,1)*gradPhi(i,j,k,1) + b2*gradPhi(i,j,k,2)*gradPhi(i,j,k,2) &
!      + b3*gradPhi(i,j,k,3)*gradPhi(i,j,k,3) - 2.*b4 - 2.*b5 - 2.*b6)/denom

!IF (denom < 1.E-13) THEN 
!   curv = 0.
!END IF


! laplacian

phiXX = grad2Phi(i,j,k,1)
phiYY = grad2Phi(i,j,k,2)
phiZZ = grad2Phi(i,j,k,3)

phiXY = gradMixPhi(i,j,k,1)
phiXZ = gradMixPhi(i,j,k,2)
phiYZ = gradMixPhi(i,j,k,3)

curv = phiXX+phiYY+phiZZ

! this allows you to set different filter sizes depending on location
!IF (gridX(i,j,k,1) > 0.) THEN
!   h = 10
!ELSE
   h = 1
!END IF

! threshold value
thresh = 0.      

pAve = phi(i,j,k)+phi(i-h,j,k)+phi(i+h,j,k)+phi(i,j+h,k)+phi(i,j-h,k)+phi(i,j,k+h)+phi(i,j,k-h)
pAve = pAve/7.

! min/max switch
IF (pAve < thresh) THEN
   F = min(curv,0.0)
ELSE
   F = max(curv,0.0)
END IF

END SUBROUTINE minMax

!*************************************************************************************!
! Calculate WENO5 derivatives
!*************************************************************************************!

SUBROUTINE weno(gM,i,j,k,nx,ny,nz,dx,phi,gradPhi,gradPhiMag)

IMPLICIT NONE

INTEGER,INTENT(IN) :: nx,ny,nz,i,j,k
REAL,INTENT(IN) :: dx
REAL,DIMENSION(0:nx,0:ny,0:nz),INTENT(INOUT) :: phi,gradPhiMag
REAL,DIMENSION(0:nx,0:ny,0:nz,3),INTENT(INOUT) :: gradPhi
REAL,INTENT(OUT):: gM
REAL :: gradX,gradY,gradZ,pa,pb,pc,pd,pe,pf,na,nb,nc,nd,ne,nf,a,b,c,d,e,f,uzm1,uzm2,uzm3
REAL :: u,ap,am,bp,bm,cp,cm,dp,dm,IS0p,IS0m,IS1p,IS1m,IS2p,IS2m,epsp,epsm,dPlus,dMinus
REAL :: a0p,a0m,a1p,a1m,a2p,a2m,w0p,w0m,w2p,w2m,PWp,PWm,p0,p1,p2,p3,p4,p5
INTEGER :: im,jm,km,ip,jp,kp,im2,ip2,jm2,jp2,km2,kp2
INTEGER :: ip1,jp1,kp1,im1,jm1,km1,ip3,jp3,kp3,im3,jm3,km3

! calculate WENO5 derivatives

IF ((i>3).AND.(i<nx-4).AND.(j>3).AND.(j<ny-4).AND.(k>3).AND.(k<nz-4)) THEN

   ! x-direction
   ap = (phi(i+3,j,k)-2.*phi(i+2,j,k)+phi(i+1,j,k))/dx
   am = (phi(i-3,j,k)-2.*phi(i-2,j,k)+phi(i-1,j,k))/dx
   bp = (phi(i+2,j,k)-2.*phi(i+1,j,k)+phi(i  ,j,k))/dx
   bm = (phi(i-2,j,k)-2.*phi(i-1,j,k)+phi(i  ,j,k))/dx
   cp = (phi(i+1,j,k)-2.*phi(i  ,j,k)+phi(i-1,j,k))/dx
   cm = cp
   dp = bm
   dm = bp

   IS0p = 13.*(ap-bp)*(ap-bp) + 3.*(ap-3.*bp)*(ap-3.*bp)
   IS0m = 13.*(am-bm)*(am-bm) + 3.*(am-3.*bm)*(am-3.*bm)
   IS1p = 13.*(bp-cp)*(bp-cp) + 3.*(bp +  cp)*(bp +  cp)
   IS1m = 13.*(bm-cm)*(bm-cm) + 3.*(bm +  cm)*(bm +  cm)
   IS2p = 13.*(cp-dp)*(cp-dp) + 3.*(3.*cp-dp)*(3.*cp-dp)
   IS2m = 13.*(cm-dm)*(cm-dm) + 3.*(3.*cm-dm)*(3.*cm-dm)

   p0 = (phi(i-2,j,k) - phi(i-3,j,k))/dx
   p1 = (phi(i-1,j,k) - phi(i-2,j,k))/dx
   p2 = (phi(i  ,j,k) - phi(i-1,j,k))/dx
   p3 = (phi(i+1,j,k) - phi(i  ,j,k))/dx
   p4 = (phi(i+2,j,k) - phi(i+1,j,k))/dx
   p5 = (phi(i+3,j,k) - phi(i+2,j,k))/dx

   ! epsilon with scaling term
   epsp = (1.E-6)*max(p1*p1,max(p2*p2,max(p3*p3,max(p4*p4,p5*p5)))) + 1.E-99
   epsm = (1.E-6)*max(p0*p0,max(p1*p1,max(p2*p2,max(p3*p3,p4*p4)))) + 1.E-99

   a0p = 1./((epsp+IS0p)*(epsp+IS0p))
   a0m = 1./((epsm+IS0m)*(epsm+IS0m))
   a1p = 6./((epsp+IS1p)*(epsp+IS1p))
   a1m = 6./((epsm+IS1m)*(epsm+IS1m))
   a2p = 3./((epsp+IS2p)*(epsp+IS2p))
   a2m = 3./((epsm+IS2m)*(epsm+IS2m))  

   w0p = a0p/(a0p+a1p+a2p)
   w0m = a0m/(a0m+a1m+a2m)
   w2p = a2p/(a0p+a1p+a2p)
   w2m = a2m/(a0m+a1m+a2m)

   PWp = 1./3.*w0p*(ap-2.*bp+cp)+1./6.*(w2p-0.5)*(bp-2.*cp+dp)
   PWm = 1./3.*w0m*(am-2.*bm+cm)+1./6.*(w2m-0.5)*(bm-2.*cm+dm)

   a = 1./12.*(-p1+7.*p2+7.*p3-p4) - PWm
   b = 1./12.*(-p1+7.*p2+7.*p3-p4) + PWp

   ! y-direction
   ap = (phi(i,j+3,k)-2.*phi(i,j+2,k)+phi(i,j+1,k))/dx
   am = (phi(i,j-3,k)-2.*phi(i,j-2,k)+phi(i,j-1,k))/dx
   bp = (phi(i,j+2,k)-2.*phi(i,j+1,k)+phi(i,j  ,k))/dx
   bm = (phi(i,j-2,k)-2.*phi(i,j-1,k)+phi(i,j  ,k))/dx
   cp = (phi(i,j+1,k)-2.*phi(i,j  ,k)+phi(i,j-1,k))/dx
   cm = cp
   dp = bm
   dm = bp
  
   IS0p = 13.*(ap-bp)*(ap-bp) + 3.*(ap-3.*bp)*(ap-3.*bp)
   IS0m = 13.*(am-bm)*(am-bm) + 3.*(am-3.*bm)*(am-3.*bm)
   IS1p = 13.*(bp-cp)*(bp-cp) + 3.*(bp +  cp)*(bp +  cp)
   IS1m = 13.*(bm-cm)*(bm-cm) + 3.*(bm +  cm)*(bm +  cm)
   IS2p = 13.*(cp-dp)*(cp-dp) + 3.*(3.*cp-dp)*(3.*cp-dp)
   IS2m = 13.*(cm-dm)*(cm-dm) + 3.*(3.*cm-dm)*(3.*cm-dm)

   p0 = (phi(i,j-2,k) - phi(i,j-3,k))/dx
   p1 = (phi(i,j-1,k) - phi(i,j-2,k))/dx
   p2 = (phi(i,j  ,k) - phi(i,j-1,k))/dx
   p3 = (phi(i,j+1,k) - phi(i,j  ,k))/dx
   p4 = (phi(i,j+2,k) - phi(i,j+1,k))/dx
   p5 = (phi(i,j+3,k) - phi(i,j+3,k))/dx

   ! epsilon with scaling term
   epsp = (1.E-6)*max(p1*p1,max(p2*p2,max(p3*p3,max(p4*p4,p5*p5)))) + 1.E-99
   epsm = (1.E-6)*max(p0*p0,max(p1*p1,max(p2*p2,max(p3*p3,p4*p4)))) + 1.E-99

   a0p = 1./((epsp+IS0p)*(epsp+IS0p))
   a0m = 1./((epsm+IS0m)*(epsm+IS0m))
   a1p = 6./((epsp+IS1p)*(epsp+IS1p))
   a1m = 6./((epsm+IS1m)*(epsm+IS1m))
   a2p = 3./((epsp+IS2p)*(epsp+IS2p))
   a2m = 3./((epsm+IS2m)*(epsm+IS2m))  

   w0p = a0p/(a0p+a1p+a2p)
   w0m = a0m/(a0m+a1m+a2m)
   w2p = a2p/(a0p+a1p+a2p)
   w2m = a2m/(a0m+a1m+a2m)

   PWp = 1./3.*w0p*(ap-2.*bp+cp)+1./6.*(w2p-0.5)*(bp-2.*cp+dp)
   PWm = 1./3.*w0m*(am-2.*bm+cm)+1./6.*(w2m-0.5)*(bm-2.*cm+dm)

   c = 1./12.*(-p1+7.*p2+7.*p3-p4) - PWm
   d = 1./12.*(-p1+7.*p2+7.*p3-p4) + PWp

   ! z-direction
   ap = (phi(i,j,k+3)-2.*phi(i,j,k+2)+phi(i,j,k+1))/dx
   am = (phi(i,j,k-3)-2.*phi(i,j,k-2)+phi(i,j,k-1))/dx
   bp = (phi(i,j,k+2)-2.*phi(i,j,k+1)+phi(i,j,k  ))/dx
   bm = (phi(i,j,k-2)-2.*phi(i,j,k-1)+phi(i,j,k  ))/dx
   cp = (phi(i,j,k+1)-2.*phi(i,j,k  )+phi(i,j,k-1))/dx
   cm = cp
   dp = bm
   dm = bp
  
   IS0p = 13.*(ap-bp)*(ap-bp) + 3.*(ap-3.*bp)*(ap-3.*bp)
   IS0m = 13.*(am-bm)*(am-bm) + 3.*(am-3.*bm)*(am-3.*bm)
   IS1p = 13.*(bp-cp)*(bp-cp) + 3.*(bp +  cp)*(bp +  cp)
   IS1m = 13.*(bm-cm)*(bm-cm) + 3.*(bm +  cm)*(bm +  cm)
   IS2p = 13.*(cp-dp)*(cp-dp) + 3.*(3.*cp-dp)*(3.*cp-dp)
   IS2m = 13.*(cm-dm)*(cm-dm) + 3.*(3.*cm-dm)*(3.*cm-dm)

   p0 = (phi(i,j,k-2) - phi(i,j,k-3))/dx
   p1 = (phi(i,j,k-1) - phi(i,j,k-2))/dx
   p2 = (phi(i,j,k  ) - phi(i,j,k-1))/dx
   p3 = (phi(i,j,k+1) - phi(i,j,k  ))/dx
   p4 = (phi(i,j,k+2) - phi(i,j,k+1))/dx
   p5 = (phi(i,j,k+3) - phi(i,j,k+2))/dx

   ! epsilon with scaling term
   epsp = (1.E-6)*max(p1*p1,max(p2*p2,max(p3*p3,max(p4*p4,p5*p5)))) + 1.E-99
   epsm = (1.E-6)*max(p0*p0,max(p1*p1,max(p2*p2,max(p3*p3,p4*p4)))) + 1.E-99

   a0p = 1./((epsp+IS0p)*(epsp+IS0p))
   a0m = 1./((epsm+IS0m)*(epsm+IS0m))
   a1p = 6./((epsp+IS1p)*(epsp+IS1p))
   a1m = 6./((epsm+IS1m)*(epsm+IS1m))
   a2p = 3./((epsp+IS2p)*(epsp+IS2p))
   a2m = 3./((epsm+IS2m)*(epsm+IS2m))  

   w0p = a0p/(a0p+a1p+a2p)
   w0m = a0m/(a0m+a1m+a2m)
   w2p = a2p/(a0p+a1p+a2p)
   w2m = a2m/(a0m+a1m+a2m)

   PWp = 1./3.*w0p*(ap-2.*bp+cp)+1./6.*(w2p-0.5)*(bp-2.*cp+dp)
   PWm = 1./3.*w0m*(am-2.*bm+cm)+1./6.*(w2m-0.5)*(bm-2.*cm+dm)

   e = 1./12.*(-p1+7.*p2+7.*p3-p4) - PWm
   f = 1./12.*(-p1+7.*p2+7.*p3-p4) + PWp

ELSE

   ! set plus and minus integers
   im = i-1
   jm = j-1
   km = k-1
   ip = i+1
   jp = j+1
   kp = k+1
            
   ! calculate derivatives
   a = (phi(i ,j ,k ) - phi(im,j ,k ))/dx
   b = (phi(ip,j ,k ) - phi(i ,j ,k ))/dx
   c = (phi(i ,j ,k ) - phi(i ,jm,k ))/dx
   d = (phi(i ,jp,k ) - phi(i ,j ,k ))/dx
   e = (phi(i ,j ,k ) - phi(i ,j ,km))/dx
   f = (phi(i ,j ,kp) - phi(i ,j ,k ))/dx
 
END IF


! check positive   
pa = max(a,0.)
pb = max(b,0.)
pc = max(c,0.)
pd = max(d,0.)
pe = max(e,0.)
pf = max(f,0.)

! check negative
na = min(a,0.)
nb = min(b,0.)
nc = min(c,0.)
nd = min(d,0.)
ne = min(e,0.)
nf = min(f,0.)

! gradient depends on which direction data is traveling
IF (phi(i,j,k) > 0.) THEN
   gradX = max(pa*pa,nb*nb)
   gradY = max(pc*pc,nd*nd)
   gradZ = max(pe*pe,nf*nf)
ELSE
   gradX = max(pb*pb,na*na)
   gradY = max(pd*pd,nc*nc)
   gradZ = max(pf*pf,ne*ne)
END IF


! fill in WENO first derivatives now
gradPhi(i,j,k,1) = gradX
gradPhi(i,j,k,2) = gradY
gradPhi(i,j,k,3) = gradZ


! magnitude of the gradient of phi
gM = sqrt(gradX+gradY+gradZ)
gradPhiMag(i,j,k) = gM

! currently not using dplus or dminus, but keeping here just incase.
dPlus = max(a,0.)**2+min(b,0.)**2+max(c,0.)**2+min(d,0.)**2+max(e,0.)**2+min(f,0.)**2
dPlus = sqrt(dPlus)
dMinus = min(a,0.)**2+max(b,0.)**2+min(c,0.)**2+max(d,0.)**2+min(e,0.)**2+max(f,0.)**2
dMinus = sqrt(dMinus)

END SUBROUTINE weno 

!*************************************************************************************!
! Reiniitialize the signed distance function
!*************************************************************************************!

SUBROUTINE reinit(phi,gradPhi,gradPhiMag,nx,ny,nz,iter,dx,h)

INTEGER,INTENT(IN) :: nx,ny,nz,iter
REAL,INTENT(IN) :: dx,h
REAL,DIMENSION(0:nx,0:ny,0:nz),INTENT(INOUT) :: phi,gradPhiMag
REAL,DIMENSION(0:nx,0:ny,0:nz,3),INTENT(INOUT) :: gradPhi
REAL :: k1,gM,phiErr
REAL,DIMENSION(0:nx,0:ny,0:nz) :: phiS,phiN
INTEGER :: i,j,k,n,raster
REAL :: sgn

raster = 0

! set the phi sign array
phiS = phi
phiN = phi

! iterate
DO n=0,iter

   !******************************* Gauss Seidel *********************************!

   ! run a raster scan for Gauss Sidel
   raster = raster + 1

   ! scan 1
   IF (raster == 1) THEN
   DO i = 1,nx-1
      DO j = 1,ny-1
         DO k = 1,nz-1
            CALL weno(gM,i,j,k,nx,ny,nz,dx,phi,gradPhi,gradPhiMag)
            CALL phiSign(phiS(i,j,k),sgn,dx,gM)
            k1  = sgn*(1.-gM) 
            phi(i,j,k) = phi(i,j,k)+h*k1
         END DO
      END DO
   END DO
   END IF

   ! scan 2
   IF (raster == 2) THEN
   DO i = 1,nx-1
      DO j = 1,ny-1
         DO k = nz-1,1,-1
            CALL weno(gM,i,j,k,nx,ny,nz,dx,phi,gradPhi,gradPhiMag)
            CALL phiSign(phiS(i,j,k),sgn,dx,gM)
            k1  = sgn*(1.-gM) 
            phi(i,j,k) = phi(i,j,k)+h*k1
         END DO
      END DO
   END DO
   END IF

   ! scan 3
   IF (raster == 3) THEN
   DO i = 1,nx-1
      DO j = ny-1,1,-1
         DO k = nz-1,1,-1
            CALL weno(gM,i,j,k,nx,ny,nz,dx,phi,gradPhi,gradPhiMag)
            CALL phiSign(phiS(i,j,k),sgn,dx,gM)
            k1  = sgn*(1.-gM) 
            phi(i,j,k) = phi(i,j,k)+h*k1
         END DO
      END DO
   END DO
   END IF

   ! scan 4
   IF (raster == 4) THEN
   DO i = nx-1,1,-1
      DO j = ny-1,1,-1
         DO k = nz-1,1,-1
            CALL weno(gM,i,j,k,nx,ny,nz,dx,phi,gradPhi,gradPhiMag)
            CALL phiSign(phiS(i,j,k),sgn,dx,gM)
            k1  = sgn*(1.-gM) 
            phi(i,j,k) = phi(i,j,k)+h*k1
         END DO
      END DO
   END DO
   END IF

   ! scan 5
   IF (raster == 5) THEN
   DO i = nx-1,1,-1
      DO j = 1,ny-1
         DO k = nz-1,1,-1
            CALL weno(gM,i,j,k,nx,ny,nz,dx,phi,gradPhi,gradPhiMag)
            CALL phiSign(phiS(i,j,k),sgn,dx,gM)
            k1  = sgn*(1.-gM) 
            phi(i,j,k) = phi(i,j,k)+h*k1
         END DO
      END DO
   END DO
   END IF

   !scan 6
   IF (raster == 6) THEN
   DO i = nx-1,1,-1
      DO j = ny-1,1,-1
         DO k = 1,nz-1
            CALL weno(gM,i,j,k,nx,ny,nz,dx,phi,gradPhi,gradPhiMag)
            CALL phiSign(phiS(i,j,k),sgn,dx,gM)
            k1  = sgn*(1.-gM) 
            phi(i,j,k) = phi(i,j,k)+h*k1
         END DO
      END DO
   END DO
   END IF

   ! scan 7
   IF (raster == 7) THEN
   DO i = nx-1,1,-1
      DO j = 1,ny-1
         DO k = 1,nz-1
            CALL weno(gM,i,j,k,nx,ny,nz,dx,phi,gradPhi,gradPhiMag)
            CALL phiSign(phiS(i,j,k),sgn,dx,gM)
            k1  = sgn*(1.-gM) 
            phi(i,j,k) = phi(i,j,k)+h*k1
         END DO
      END DO
   END DO
   END IF

   ! scan 8
   IF (raster == 8) THEN
   DO i = 1,nx-1
      DO j = ny-1,1,-1
         DO k = 1,nz-1
            CALL weno(gM,i,j,k,nx,ny,nz,dx,phi,gradPhi,gradPhiMag)
            CALL phiSign(phiS(i,j,k),sgn,dx,gM)
            k1  = sgn*(1.-gM) 
            phi(i,j,k) = phi(i,j,k)+h*k1
         END DO
      END DO
   END DO
   END IF

   ! reset counter
   IF (raster == 8) raster = 0


   ! extrapolation boundary condition
   DO i = 0,nx
      DO j = 0,ny
         DO k = 0,nz

            ! corners
            phi(0,0,0) = phi(1,1,1) + dx 
            phi(nx,0,0) = phi(nx-1,1,1) + dx
            phi(0,ny,0) = phi(1,ny-1,1) + dx 
            phi(0,0,nz) = phi(1,1,nz-1) + dx 
            phi(nx,ny,0) = phi(nx-1,ny-1,1) + dx
            phi(0,ny,nz) = phi(1,ny-1,nz-1) + dx 
            phi(nx,0,nz) = phi(nx-1,1,nz-1) + dx 
            phi(nx,ny,nz) = phi(nx-1,ny-1,nz-1) + dx 

            ! edges
            phi(i,0,0) = phi(i,1,1) + dx 
            phi(0,j,0) = phi(1,j,1) + dx 
            phi(0,0,k) = phi(1,1,k) + dx 
            phi(i,ny,nz) = phi(i,ny-1,nz-1) + dx 
            phi(nx,j,nz) = phi(nx-1,j,nz-1) + dx 
            phi(nx,ny,k) = phi(nx-1,ny-1,k) + dx 
            phi(i,0,nz) = phi(i,1,nz-1) + dx 
            phi(nx,j,0) = phi(nx-1,j,1) + dx 
            phi(nx,0,k) = phi(nx-1,1,k) + dx
            phi(i,ny,0) = phi(i,ny-1,1) + dx 
            phi(0,j,nz) = phi(1,j,nz-1) + dx 
            phi(0,ny,k) = phi(1,ny-1,k) + dx 
            
            ! faces
            phi(0,j,k) = phi(1,j,k) + dx 
            phi(i,0,k) = phi(i,1,k) + dx 
            phi(i,j,0) = phi(i,j,1) + dx  
            phi(nx,j,k) = phi(nx-1,j,k) + dx 
            phi(i,ny,k) = phi(i,ny-1,k) + dx             
            phi(i,j,nz) = phi(i,j,nz-1) + dx

         END DO
      END DO
   END DO

   
   !********************************* RMS ***************************************!

   phiErr = 0.
  
   ! calculate RMS
   DO i = 0,nx
      DO j = 0,ny
         DO k = 0,nz
            phiErr = phiErr + (phi(i,j,k)-phiN(i,j,k))*(phi(i,j,k)-phiN(i,j,k))
         END DO
      END DO
   END DO

   ! check error
   phiErr = sqrt(phiErr/(nx*ny*nz))
   IF (phiErr < 1.E-5) THEN
      PRINT*, " Distance function time integration has reached steady state "
      EXIT
   END IF

   ! set new phi 
   phiN = phi

   PRINT*, " Iteration: ",n," ", " RMS Error: ",phiErr

   ! check for NAN
   IF (isnan(phiErr)) STOP 

END DO
PRINT*,

END SUBROUTINE reinit


!*************************************************************************************!
! Trilinear interpolation for min/max speed
!*************************************************************************************!

SUBROUTINE setSurfCurv(xLo,nx,ny,nz,dx,curvSurf,curv,nSurfNode,surfX,gradPhiSurf,gradPhi)

INTEGER,INTENT(IN) :: nSurfNode,nx,ny,nz
REAL,INTENT(IN) :: dx
REAL,DIMENSION(3),INTENT(IN) :: xLo
REAL,DIMENSION(0:nx,0:ny,0:nz),INTENT(IN) :: curv
REAL,DIMENSION(0:nx,0:ny,0:nz,3),INTENT(IN) :: gradPhi
REAL,DIMENSION(nSurfNode),INTENT(INOUT) :: curvSurf
REAL,DIMENSION(nSurfNode,3),INTENT(IN) :: surfX
REAL,DIMENSION(nSurfNode,3),INTENT(OUT) :: gradPhiSurf
INTEGER :: n,i1,j1,k1,i0,j0,k0
REAL :: x,y,z,xd,yd,zd,c00,c10,c01,c11,c0,c1,x1,y1,z1,x0,y0,z0
REAL :: gradPhiMag,gradMag2,minX,minY,minZ

minX = xLo(1)
minY = xLo(2)
minZ = xLo(3)

DO n = 1,nSurfNode

   x = surfX(n,1)
   y = surfX(n,2)
   z = surfX(n,3)

   i0 = floor((x-minX)/dx)
   j0 = floor((y-minY)/dx)
   k0 = floor((z-minZ)/dx)

   x0 = i0*dx + minX
   y0 = j0*dx + minY
   z0 = k0*dx + minZ

   i1 = i0+1
   j1 = j0+1
   k1 = k0+1

   x1 = i1*dx + minX
   y1 = j1*dx + minY
   z1 = k1*dx + minZ

   xd = (x-x0)/(x1-x0)
   yd = (y-y0)/(y1-y0)
   zd = (z-z0)/(z1-z0)

   c00 = curv(i0,j0,k0)*(1.-xd) + curv(i1,j0,k0)*xd
   c10 = curv(i0,j1,k0)*(1.-xd) + curv(i1,j1,k0)*xd
   c01 = curv(i0,j0,k1)*(1.-xd) + curv(i1,j0,k1)*xd
   c11 = curv(i0,j1,k1)*(1.-xd) + curv(i1,j1,k1)*xd

   c0 = c00*(1.-yd)+c10*yd
   c1 = c01*(1.-yd)+c11*yd

   curvSurf(n) = c0*(1.-zd)+c1*zd

   c00 = gradPhi(i0,j0,k0,1)*(1.-xd) + gradPhi(i1,j0,k0,1)*xd
   c10 = gradPhi(i0,j1,k0,1)*(1.-xd) + gradPhi(i1,j1,k0,1)*xd
   c01 = gradPhi(i0,j0,k1,1)*(1.-xd) + gradPhi(i1,j0,k1,1)*xd
   c11 = gradPhi(i0,j1,k1,1)*(1.-xd) + gradPhi(i1,j1,k1,1)*xd

   c0 = c00*(1.-yd)+c10*yd
   c1 = c01*(1.-yd)+c11*yd

   IF (x > 0.) THEN
      gradPhiSurf(n,1) = -(c0*(1.-zd)+c1*zd)
   ELSE 
      gradPhiSurf(n,1) = (c0*(1.-zd)+c1*zd)
   END IF

   c00 = gradPhi(i0,j0,k0,2)*(1.-xd) + gradPhi(i1,j0,k0,2)*xd
   c10 = gradPhi(i0,j1,k0,2)*(1.-xd) + gradPhi(i1,j1,k0,2)*xd
   c01 = gradPhi(i0,j0,k1,2)*(1.-xd) + gradPhi(i1,j0,k1,2)*xd
   c11 = gradPhi(i0,j1,k1,2)*(1.-xd) + gradPhi(i1,j1,k1,2)*xd

   c0 = c00*(1.-yd)+c10*yd
   c1 = c01*(1.-yd)+c11*yd

   IF (y > 0.) THEN
      gradPhiSurf(n,2) = -(c0*(1.-zd)+c1*zd)
   ELSE 
      gradPhiSurf(n,2) = (c0*(1.-zd)+c1*zd)
   END IF

   c00 = gradPhi(i0,j0,k0,3)*(1.-xd) + gradPhi(i1,j0,k0,3)*xd
   c10 = gradPhi(i0,j1,k0,3)*(1.-xd) + gradPhi(i1,j1,k0,3)*xd
   c01 = gradPhi(i0,j0,k1,3)*(1.-xd) + gradPhi(i1,j0,k1,3)*xd
   c11 = gradPhi(i0,j1,k1,3)*(1.-xd) + gradPhi(i1,j1,k1,3)*xd

   c0 = c00*(1.-yd)+c10*yd;
   c1 = c01*(1.-yd)+c11*yd;

   IF (z > 0.) THEN
      gradPhiSurf(n,3) = -(c0*(1.-zd)+c1*zd)
   ELSE 
      gradPhiSurf(n,3) = (c0*(1.-zd)+c1*zd)
   END IF

   gradMag2 = gradPhiSurf(n,1)*gradPhiSurf(n,1)+gradPhiSurf(n,2)*gradPhiSurf(n,2)+gradPhiSurf(n,3)*gradPhiSurf(n,3)

   IF (gradMag2 < 1.E-7) THEN
      gradPhiMag = 0.
      gradPhiSurf(n,1) = 0.
      gradPhiSurf(n,2) = 0.
      gradPhiSurf(n,3) = 0.
   ELSE
      gradPhiMag = sqrt(gradMag2)
      gradPhiSurf(n,1) = gradPhiSurf(n,1)/gradPhiMag
      gradPhiSurf(n,2) = gradPhiSurf(n,2)/gradPhiMag
      gradPhiSurf(n,3) = gradPhiSurf(n,3)/gradPhiMag
   END IF

END DO

END SUBROUTINE setSurfCurv


!*************************************************************************************!
! Trilinear interpolation for phi
!*************************************************************************************!

SUBROUTINE setPhiSurf(xLo,nx,ny,nz,dx,phiSurf,phi,nSurfNode,surfX,gradPhiSurf,gradPhi)

INTEGER,INTENT(IN) :: nSurfNode,nx,ny,nz
REAL,INTENT(IN) :: dx
REAL,DIMENSION(3),INTENT(IN) :: xLo
REAL,DIMENSION(0:nx,0:ny,0:nz),INTENT(IN) :: phi
REAL,DIMENSION(0:nx,0:ny,0:nz,3),INTENT(IN) :: gradPhi
REAL,DIMENSION(nSurfNode),INTENT(INOUT) :: phiSurf
REAL,DIMENSION(nSurfNode,3),INTENT(IN) :: surfX
REAL,DIMENSION(nSurfNode,3),INTENT(OUT) :: gradPhiSurf
INTEGER :: n,i1,j1,k1,i0,j0,k0
REAL :: x,y,z,xd,yd,zd,c00,c10,c01,c11,c0,c1,x1,y1,z1,x0,y0,z0
REAL :: gradPhiMag,gradMag2,minX,minY,minZ

minX = xLo(1)
minY = xLo(2)
minZ = xLo(3)


DO n = 1,nSurfNode

   x = surfX(n,1)
   y = surfX(n,2)
   z = surfX(n,3)

   i0 = floor((x-minX)/dx)
   j0 = floor((y-minY)/dx)
   k0 = floor((z-minZ)/dx)

   x0 = i0*dx + minX
   y0 = j0*dx + minY
   z0 = k0*dx + minZ

   i1 = i0+1
   j1 = j0+1
   k1 = k0+1

   x1 = i1*dx + minX
   y1 = j1*dx + minY
   z1 = k1*dx + minZ

   xd = (x-x0)/(x1-x0)
   yd = (y-y0)/(y1-y0)
   zd = (z-z0)/(z1-z0)

   c00 = phi(i0,j0,k0)*(1.-xd) + phi(i1,j0,k0)*xd
   c10 = phi(i0,j1,k0)*(1.-xd) + phi(i1,j1,k0)*xd
   c01 = phi(i0,j0,k1)*(1.-xd) + phi(i1,j0,k1)*xd
   c11 = phi(i0,j1,k1)*(1.-xd) + phi(i1,j1,k1)*xd

   c0 = c00*(1.-yd)+c10*yd
   c1 = c01*(1.-yd)+c11*yd

   phiSurf(n) = c0*(1.-zd)+c1*zd

   c00 = gradPhi(i0,j0,k0,1)*(1.-xd) + gradPhi(i1,j0,k0,1)*xd
   c10 = gradPhi(i0,j1,k0,1)*(1.-xd) + gradPhi(i1,j1,k0,1)*xd
   c01 = gradPhi(i0,j0,k1,1)*(1.-xd) + gradPhi(i1,j0,k1,1)*xd
   c11 = gradPhi(i0,j1,k1,1)*(1.-xd) + gradPhi(i1,j1,k1,1)*xd

   c0 = c00*(1.-yd)+c10*yd
   c1 = c01*(1.-yd)+c11*yd

   !IF (x > 0.) THEN
      gradPhiSurf(n,1) = -(c0*(1.-zd)+c1*zd)
   !ELSE 
   !   gradPhiSurf(n,1) = (c0*(1.-zd)+c1*zd)
   !END IF

   c00 = gradPhi(i0,j0,k0,2)*(1.-xd) + gradPhi(i1,j0,k0,2)*xd
   c10 = gradPhi(i0,j1,k0,2)*(1.-xd) + gradPhi(i1,j1,k0,2)*xd
   c01 = gradPhi(i0,j0,k1,2)*(1.-xd) + gradPhi(i1,j0,k1,2)*xd
   c11 = gradPhi(i0,j1,k1,2)*(1.-xd) + gradPhi(i1,j1,k1,2)*xd

   c0 = c00*(1.-yd)+c10*yd
   c1 = c01*(1.-yd)+c11*yd

   !IF (y > 0.) THEN
      gradPhiSurf(n,2) = -(c0*(1.-zd)+c1*zd)
   !ELSE 
   !   gradPhiSurf(n,2) = (c0*(1.-zd)+c1*zd)
   !END IF

   c00 = gradPhi(i0,j0,k0,3)*(1.-xd) + gradPhi(i1,j0,k0,3)*xd
   c10 = gradPhi(i0,j1,k0,3)*(1.-xd) + gradPhi(i1,j1,k0,3)*xd
   c01 = gradPhi(i0,j0,k1,3)*(1.-xd) + gradPhi(i1,j0,k1,3)*xd
   c11 = gradPhi(i0,j1,k1,3)*(1.-xd) + gradPhi(i1,j1,k1,3)*xd

   c0 = c00*(1.-yd)+c10*yd;
   c1 = c01*(1.-yd)+c11*yd;

   !IF (z > 0.) THEN
      gradPhiSurf(n,3) = -(c0*(1.-zd)+c1*zd)
   !ELSE 
   !   gradPhiSurf(n,3) = (c0*(1.-zd)+c1*zd)
   !END IF

   gradMag2 = gradPhiSurf(n,1)*gradPhiSurf(n,1)+gradPhiSurf(n,2)*gradPhiSurf(n,2)+gradPhiSurf(n,3)*gradPhiSurf(n,3)

   IF (gradMag2 < 1.E-7) THEN
      gradPhiMag = 0.
      gradPhiSurf(n,1) = 0.
      gradPhiSurf(n,2) = 0.
      gradPhiSurf(n,3) = 0.
   ELSE
      gradPhiMag = sqrt(gradMag2)
      gradPhiSurf(n,1) = gradPhiSurf(n,1)/gradPhiMag
      gradPhiSurf(n,2) = gradPhiSurf(n,2)/gradPhiMag
      gradPhiSurf(n,3) = gradPhiSurf(n,3)/gradPhiMag
   END IF

END DO

END SUBROUTINE setPhiSurf

!*************************************************************************************!
!
! Module End
!
!*************************************************************************************!


END MODULE set_subs
