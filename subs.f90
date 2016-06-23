MODULE set_subs

IMPLICIT NONE

CONTAINS

!*************************************************************************************!
! Generate Strands
!*************************************************************************************!

SUBROUTINE strandGen(triangles,nTri,nSurfNode,normals,clip,strandDist,wSpace,smooth,nodeNormals)

REAL*4,DIMENSION(:,:),INTENT(IN) :: triangles,normals
INTEGER,INTENT(IN) :: nTri,clip,strandDist
INTEGER,INTENT(INOUT) :: nSurfNode
REAL,INTENT(IN) :: wSpace,smooth

INTEGER :: n,i,k,kk,p,share
INTEGER,DIMENSION(nSurfNode) :: nSurT,nNum   !! Number of Surrounding Triangles, Node Number
INTEGER*4,DIMENSION(nTri,3) :: surfElem      !! Tag for surface nodes
REAL*4,DIMENSION(3,nTri*5)  :: nodesT        !! Node Locations
REAL*4,DIMENSION(3,nTri),INTENT(OUT) :: nodeNormals   !! Cumulative area for surrounding triangles

!! Determine surface nodes, and shared nodes from triangles

nSurT = 1
nodeNormals  = 0.

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
            nSurT(kk) = nSurT(kk) + 1
            nodeNormals(:,kk) = nodeNormals(:,kk) + normals(:,kk)
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

nodeNormals(:,:) = nodeNormals(:,:)/(sqrt(nodeNormals(:,:)*nodeNormals(:,:)))

END SUBROUTINE strandGen

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
         IF (abs(phi(i,j,k)) < 3.1*dx) THEN
            phiNB(i,j,k) = 1
         END IF

         ! stencil band
         IF (abs(phi(i,j,k)) < 5.1*dx) THEN
            phiSB(i,j,k) = 1
         END IF

      END DO
   END DO
END DO

END SUBROUTINE narrowBand

!*************************************************************************************!
! Calculate First Derivative
!*************************************************************************************!

SUBROUTINE firstDeriv(i,j,k,nx,ny,nz,dx,phi,phiX,phiY,phiZ,order,gMM)

INTEGER,INTENT(IN) :: i,j,k,order,nx,ny,nz
REAL,INTENT(IN) :: dx
REAL,DIMENSION(0:nx,0:ny,0:nz),INTENT(IN) :: phi
REAL,INTENT(OUT) :: phiX,phiY,phiZ,gMM
REAL :: aa1,aa2,aa3,aa4,aa5,aa6,aa7,dxx
INTEGER :: im,jm,km,ip,jp,kp,im2,ip2,jm2,jp2,km2,kp2
INTEGER :: ip1,jp1,kp1,im1,jm1,km1,ip3,jp3,kp3,im3,jm3,km3

IF (order == 2) THEN
  
   dxx = 1./(dx); 

   ! second order first derivative
   phiX = (-1./2.*phi(i-1,j,k) +  1./2.*phi(i+1,j,k))*dxx
   phiY = (-1./2.*phi(i,j-1,k) +  1./2.*phi(i,j+1,k))*dxx
   phiZ = (-1./2.*phi(i,j,k-1) +  1./2.*phi(i,j,k+1))*dxx

ELSE 
   PRINT*," This order not set yet "
   STOP
END IF


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
REAL :: phiXX,phiYY,phiZZ,phiXY,phiXZ,phiYZ
INTEGER :: h
INTEGER,INTENT(IN) :: i,j,k,nx,ny,nz
REAL,INTENT(IN) :: dx
REAL,DIMENSION(0:nx,0:ny,0:nz),INTENT(IN) :: phi
REAL,DIMENSION(0:nx,0:ny,0:nz,3),INTENT(IN) :: grad2Phi,gradMixPhi,gridX
REAL,INTENT(OUT) :: F
       
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
REAL,DIMENSION(0:nx,0:ny,0:nz,3,2),INTENT(INOUT) :: gradPhi
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
   IS1m = 13.*(bm-cm)*(bm-cm) + 3.*(bm +  bm)*(bm +  bm)
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
   IS1m = 13.*(bm-cm)*(bm-cm) + 3.*(bm +  bm)*(bm +  bm)
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
   IS1m = 13.*(bm-cm)*(bm-cm) + 3.*(bm +  bm)*(bm +  bm)
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


! fill in WENO first derivatives now
gradPhi(i,j,k,1,1) = a
gradPhi(i,j,k,2,1) = b
gradPhi(i,j,k,3,1) = c
gradPhi(i,j,k,1,2) = d
gradPhi(i,j,k,2,2) = e
gradPhi(i,j,k,3,2) = f

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
REAL,DIMENSION(0:nx,0:ny,0:nz,3,2),INTENT(INOUT) :: gradPhi
REAL :: k1,gM,phiErr
REAL,DIMENSION(0:nx,0:ny,0:nz) :: phiS,phiN
INTEGER :: i,j,k,n
REAL :: sgn

! set the phi sign array
phiS = phi

! iterate
DO n=0,iter

   !********************* Explicit Forward Euler Scheme **************************!

   DO i = 1,nx-1
      DO j = 1,ny-1
         DO k = 1,nz-1
            CALL weno(gM,i,j,k,nx,ny,nz,dx,phi,gradPhi,gradPhiMag)
            CALL phiSign(phiS(i,j,k),sgn,dx,gM)
            k1  = sgn*(1.-gM) 
            phiN(i,j,k) = phi(i,j,k)+h*k1
         END DO
      END DO
   END DO

   ! extrapolation boundary condition
   DO i = 0,nx
      DO j = 0,ny
         DO k = 0,nz

            ! corners
            phiN(0,0,0) = phi(1,1,1) + dx 
            phiN(nx,0,0) = phi(nx-1,1,1) + dx
            phiN(0,ny,0) = phi(1,ny-1,1) + dx 
            phiN(0,0,nz) = phi(1,1,nz-1) + dx 
            phiN(nx,ny,0) = phi(nx-1,ny-1,1) + dx
            phiN(0,ny,nz) = phi(1,ny-1,nz-1) + dx 
            phiN(nx,0,nz) = phi(nx-1,1,nz-1) + dx 
            phiN(nx,ny,nz) = phi(nx-1,ny-1,nz-1) + dx 

            ! edges
            phiN(i,0,0) = phi(i,1,1) + dx 
            phiN(0,j,0) = phi(1,j,1) + dx 
            phiN(0,0,k) = phi(1,1,k) + dx 
            phiN(i,ny,nz) = phi(i,ny-1,nz-1) + dx 
            phiN(nx,j,nz) = phi(nx-1,j,nz-1) + dx 
            phiN(nx,ny,k) = phi(nx-1,ny-1,k) + dx 
            phiN(i,0,nz) = phi(i,1,nz-1) + dx 
            phiN(nx,j,0) = phi(nx-1,j,1) + dx 
            phiN(nx,0,k) = phi(nx-1,1,k) + dx
            phiN(i,ny,0) = phi(i,ny-1,1) + dx 
            phiN(0,j,nz) = phi(1,j,nz-1) + dx 
            phiN(0,ny,k) = phi(1,ny-1,k) + dx 
            
            ! faces
            phiN(0,j,k) = phi(1,j,k) + dx 
            phiN(i,0,k) = phi(i,1,k) + dx 
            phiN(i,j,0) = phi(i,j,1) + dx  
            phiN(nx,j,k) = phi(nx-1,j,k) + dx 
            phiN(i,ny,k) = phi(i,ny-1,k) + dx             
            phiN(i,j,nz) = phi(i,j,nz-1) + dx

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
   phi = phiN

   PRINT*, " Iteration: ",n," ", " RMS Error: ",phiErr

   ! check for NAN
   IF (isnan(phiErr)) STOP 

END DO
PRINT*,

END SUBROUTINE reinit

END MODULE set_subs
