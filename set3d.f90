PROGRAM set3d

!*************************************************************************************!
!
! Min/Max Flow Level Set Solver
!
! Version: 3.0  
!
! Date: 20/1/2016
!
! Type: H-J WENO 5 Serial Explicit TVD R-K Level Set Solver
!
!*************************************************************************************!


IMPLICIT NONE

!*************************************************************************************!
! Set Data
!*************************************************************************************!

INTEGER :: j,sUnit,nbytePhi,offset,nxs,n1,n2,n3,fN,im,jm,km,ip,jp,kp,dd
INTEGER :: iter,nx,ny,nz,qG,nn,counter,orderUp,order1,order2
REAL :: a0,a1,a2,a3,a4,t,xLo(3),xHi(3),xs,ys,yp,ypp,g,gp,x0,x1,dxs,dPlus,dMinus
REAL :: pX1,pX2,pX3,pY1,pY2,pY3,pZ1,pZ2,pZ3,gX,gY,gZ,dis,minD,k1,h1,k2,k3,k4,t5
REAL :: B1,B2,B3,C1,C2,C3,pSx,pSy,pSz,pS,pX,pY,pZ,sgn,ddx,ddy,ddz,gM,gMM,phiErrS
REAL :: y1,z1,maxX,maxY,maxZ,minX,minY,minZ,h,a,b,c,d,e,f,phiErr,gradMag2,x,y,z,dxx
REAL :: dx,cX,cY,cZ,t1,t2,t3,t4,bx,phiX,phiY,phiZ,phiXX,phiYY,phiZZ,phiXZ,phiXY,phiYZ
REAL :: pi,aa,bb,period,ax,ay,az,axy,ayz,axz,bxy,by,bz,byz,bxz,cxy,cyz,cxz,pp,CFL
REAL :: aaa,bbb,ccc,ax1,ax2,ay1,ay2,az1,az2,Dijk
REAL,ALLOCATABLE,DIMENSION(:,:,:) :: phi,phiS,phiN,gradPhiMag,phiE,phi1,phi2,phi3,phiO,S
INTEGER,ALLOCATABLE,DIMENSION(:,:,:) :: phiSB,phiNB
REAL,ALLOCATABLE,DIMENSION(:,:) :: centroid,surfX
CHARACTER(LEN=1) :: lf=char(10)
CHARACTER(LEN=1024) :: extent,origin,spacing,coffset
INTEGER*4,ALLOCATABLE,DIMENSION(:,:) :: surfElem
CHARACTER header*80,filename*80
INTEGER*2 padding
INTEGER*4 ntri,iunit,nSurfNode,k,i,n,p,kk,share,nSurfElem
REAL*4,ALLOCATABLE,DIMENSION(:,:) :: normals,triangles,nodesT
REAL,ALLOCATABLE,DIMENSION(:,:,:,:) :: gridX,gradPhi,grad2Phi,gradMixPhi
INTEGER iargc
REAL,ALLOCATABLE,DIMENSION(:,:) :: q1S,q2S,q3S,q4S,q5S,q6S,q7S,q8S
INTEGER,DIMENSION(8) :: qq


!*************************************************************************************!
! Import STL Data
!*************************************************************************************!

! start system time
CALL cpu_time(t1)

! import .stl data
call getarg(1,filename)

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


!*************************************************************************************!
! Determine xLo and xHi
!*************************************************************************************!

! initialize
x1 = surfX(1,1);
y1 = surfX(1,2);
z1 = surfX(1,3);

maxX = x1;
maxY = y1;
maxZ = z1;

minX = x1;
minY = y1;
minZ = z1;

! find the max and min
DO n = 2,nSurfNode 
   x1 = surfX(n,1)
   y1 = surfX(n,2)
   z1 = surfX(n,3)

   IF (x1 > maxX) THEN
      maxX = x1
   END IF
   IF (y1 > maxY) THEN
     maxY = y1
   END IF
   IF (z1 > maxZ) THEN
      maxZ = z1
   END IF

   IF (x1 < minX) THEN
      minX = x1
   END IF
   IF (y1 < minY) THEN
     minY = y1
   END IF
   IF (z1 < minZ) THEN
      minZ = z1
   END IF
   
END DO


!*************************************************************************************!
! Define Cartesiang Grid Size and Allocate Phi
!*************************************************************************************!

! find the characteristic size of the object to add around
ddx = maxX-minX
ddy = maxY-minY
ddz = maxZ-minZ

! set dx
dx = 0.05!*100.
!dx = 0.003125*100.

! define Cartesian grid
nx = ceiling((maxX-minX)/dx)+1;
ny = ceiling((maxY-minY)/dx)+1;
nz = ceiling((maxZ-minZ)/dx)+1;

! number of cells you want to add
dd = 10.

! adding more cells edge
nx = nx+2*dd
ny = ny+2*dd
nz = nz+2*dd

! set xLo and xHi
xLo =(/minX-dd*dx,minY-dd*dx,minZ-dd*dx/)
xHi =(/maxX+dd*dx,maxY+dd*dx,maxZ+dd*dx/)

! allocate phi
ALLOCATE(phi(0:nx,0:ny,0:nz))
phi = 1.

! allocate a grid of x,y,z points 
ALLOCATE(gridX(0:nx,0:ny,0:nz,3))
DO i = 0,nx
   DO j = 0,ny
      DO k = 0,nz
         gridX(i,j,k,1) = xLo(1) + i*dx;
         gridX(i,j,k,2) = xLo(2) + j*dx;
         gridX(i,j,k,3) = xLo(3) + k*dx;
      END DO
   END DO
END DO

!*************************************************************************************!
! Determine Inside and Outside of Surface
!*************************************************************************************!

! cut down on excess and use minimum nodes
im = floor((minX-xLo(1))/dx)-3
jm = floor((minY-xLo(2))/dx)-3
km = floor((minZ-xLo(3))/dx)-3
! cut down on excess and use maximum nodes
ip = floor((maxX-xLo(1))/dx)+3
jp = floor((maxY-xLo(2))/dx)+3
kp = floor((maxZ-xLo(3))/dx)+3

! print out grid spacing
PRINT*, " Setting Grid Size "
PRINT*, " Grid Size: nx =",nx,", ny =",ny,",nz =",nz 
PRINT*, " Grid Spacing: dx =", dx
PRINT*,
PRINT*, " Determining Inside and Outside of Geometry "
PRINT*, 


! allocate centroid
ALLOCATE(centroid(nSurfElem,4))

DO n = 1,nSurfElem
   n1 = surfElem(n,1)
   n2 = surfElem(n,2)
   n3 = surfElem(n,3)
   pX1 = surfX(n1,1)
   pY1 = surfX(n1,2)
   pZ1 = surfX(n1,3)
   pX2 = surfX(n2,1)
   pY2 = surfX(n2,2)
   pZ2 = surfX(n2,3)
   pX3 = surfX(n3,1)
   pY3 = surfX(n3,2)
   pZ3 = surfX(n3,3)
   centroid(n,1) = (pX1+pX2+pX3)/3.
   centroid(n,2) = (pY1+pY2+pY3)/3.
   centroid(n,3) = (pZ1+pZ2+pZ3)/3.
END DO

! find which nodes are inside and outside
DO i = im,ip 
   DO j = jm,jp  
      DO k = km,kp  

         ! search through all surface elements to find closest element
         minD = 100000.;
         DO n = 1,nSurfElem
            pX = centroid(n,1)
            pY = centroid(n,2)
            pZ = centroid(n,3)
            gX = gridX(i,j,k,1)
            gY = gridX(i,j,k,2)
            gZ = gridX(i,j,k,3)
            dis = sqrt((pX-gX)*(pX-gX) + (pY-gY)*(pY-gY) + (pZ-gZ)*(pZ-gZ))
            IF (dis < minD) THEN
               minD = dis;
               fN   = n;
            END IF
         END DO

      ! create three vectors from our point to the surfElem points
      n1 = surfElem(fN,1)
      n2 = surfElem(fN,2)
      n3 = surfElem(fN,3)
      A1 = surfX(n1,1) - gridX(i,j,k,1)
      A2 = surfX(n1,2) - gridX(i,j,k,2)
      A3 = surfX(n1,3) - gridX(i,j,k,3)
      B1 = surfX(n2,1) - gridX(i,j,k,1)
      B2 = surfX(n2,2) - gridX(i,j,k,2)
      B3 = surfX(n2,3) - gridX(i,j,k,3)
      C1 = surfX(n3,1) - gridX(i,j,k,1)
      C2 = surfX(n3,2) - gridX(i,j,k,2)
      C3 = surfX(n3,3) - gridX(i,j,k,3)
    
      ! cross product two of the vectors
      pSx = A2*B3-A3*B2;
      pSy = -(A1*B3-B1*A3);
      pSz = A1*B2-B1*A2;

      ! dot product last vector with your cross product result
      pS =-(pSx*C1+pSy*C2+pSz*C3);
      
      gM = 1.
      ! return sign of the dot product
      CALL phiSign(pS,sgn,dx,gM)
              
      phi(i,j,k) = sgn

      END DO
   END DO
END DO


! deallocate surface mesh data
DEALLOCATE(surfX)
!DEALLOCATE(gridX)
DEALLOCATE(surfElem)
DEALLOCATE(centroid)

CALL cpu_time(t2)
PRINT*, " Search Run Time: ",t2-t1," Seconds"
PRINT*,

!*************************************************************************************!
! Fast Marching Method
!*************************************************************************************!

! allocate arrays used for FMM
ALLOCATE(phiS(0:nx,0:ny,0:nz))
ALLOCATE(phiO(0:nx,0:ny,0:nz))
ALLOCATE(phiN(0:nx,0:ny,0:nz))
ALLOCATE(gradPhiMag(0:nx,0:ny,0:nz))
ALLOCATE(phiE(0:nx,0:ny,0:nz))
ALLOCATE(phi1(0:nx,0:ny,0:nz))
ALLOCATE(phi2(0:nx,0:ny,0:nz))
ALLOCATE(phi3(0:nx,0:ny,0:nz))
ALLOCATE(S(0:nx,0:ny,0:nz))

! source term is zero
S = 0.


PRINT*, " Level Set Time Integration "
PRINT*, 


! set the phi sign array
phiS = phi

! number of iterations
iter = 10000 !1500 ! 10000

! normalized dx
dxx = dx/sqrt(ddx*ddx+ddy*ddy+ddz*ddz)

! time step
CFL = .1
h = CFL*dxx

! iterate
DO n=0,iter

   !********************* Explicit Forward Euler Scheme **************************!

   DO i = 1,nx-1
      DO j = 1,ny-1
         DO k = 1,nz-1
            CALL weno(gM,i,j,k,nx,ny,nz,dx,phi,dPlus,dMinus)
            CALL phiSign(phiS(i,j,k),sgn,dx,gM)
            k1  = sgn*(1.-gM) 
            phiN(i,j,k) = phi(i,j,k)+h*k1
         END DO
      END DO
   END DO


   !****************** Explicit Third Order TVD RK Scheme ************************!

!   DO i = 1,nx-1
!      DO j = 1,ny-1
!         DO k = 1,nz-1
!            CALL weno(gM,i,j,k,nx,ny,nz,dx,phi,dPlus,dMinus)
!            CALL phiSign(phiS(i,j,k),sgn,dxx,gM)
!            k1  = sgn*(1.-gM) 
!            !h = CFL*dxx*abs(1.-gM)
!            phi1(i,j,k) = phi(i,j,k)+h*k1
!         END DO
!      END DO
!   END DO
!
!   DO i = 1,nx-1
!      DO j = 1,ny-1
!         DO k = 1,nz-1
!            CALL weno(gM,i,j,k,nx,ny,nz,dx,phi,dPlus,dMinus)
!            CALL phiSign(phiS(i,j,k),sgn,dxx,gM)
!            k2  = sgn*(1.-gM) 
!            !h = CFL*dxx*abs(1.-gM)
!            phi2(i,j,k) = 3./4.*phi(i,j,k)+1./4.*phi1(i,j,k)+1./4.*h*k2
!         END DO
!      END DO
!   END DO
!
!
!   DO i = 1,nx-1
!      DO j = 1,ny-1
!         DO k = 1,nz-1
!            CALL weno(gM,i,j,k,nx,ny,nz,dx,phi,dPlus,dMinus)
!            CALL phiSign(phiS(i,j,k),sgn,dxx,gM)
!            k3  = sgn*(1.-gM)
!           !h = CFL*dxx*abs(1.-gM)
!            phiN(i,j,k) = 1./3.*phi(i,j,k)+2./3.*phi2(i,j,k)+2./3.*h*k3
!         END DO
!      END DO
!   END DO
 

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

! set original phi
phiO = phi

! print out run time
CALL cpu_time(t3)
PRINT*, " Initialization Run Time: ",t3-t1," Seconds"
PRINT*,


!*************************************************************************************!
! Paraview Output
!*************************************************************************************!


PRINT*, " Writing Out Cartesian Grid to Paraview Format "
PRINT*, 


! output to Paraview
WRITE(extent,'(3(A3,I6))')' 0 ',nx,' 0 ',ny,' 0 ',nz
WRITE(origin,'(3(F20.8,A1))')xLo(1),' ',xLo(2),' ',xLo(3),' '
WRITE(spacing,'(3(F20.8,A1))')dx,' ',dx,' ',dx,' '
nbytePhi =(nx+1)**3*24
offset = 0
WRITE(coffset,'(I16)')offset


sUnit = 11
OPEN(UNIT=sUnit,FILE='blockBefore.vti',FORM='unformatted',ACCESS='stream',STATUS='replace')
WRITE(sUnit)'<?xml version="1.0"?>'//lf
WRITE(sUnit)'<VTKFile type="ImageData" version="0.1" byte_order="LittleEndian">'//lf
WRITE(sUnit)'<ImageData WholeExtent="',TRIM(extent),'" Origin="',TRIM(origin),'" Spacing="',TRIM(spacing),'">'//lf
WRITE(sUnit)'<Piece Extent="',TRIM(extent),'">'//lf
WRITE(sUnit)'<PointData Scalars="phi">'//lf
WRITE(sUnit)'<DataArray type="Float64" Name="phi" format="appended" offset="',TRIM(coffset),'"/>'//lf
WRITE(sUnit)'</PointData>'//lf
WRITE(sUnit)'</Piece>'//lf
WRITE(sUnit)'</ImageData>'//lf
WRITE(sUnit)'<AppendedData encoding="raw">'//lf
WRITE(sUnit)'_'
WRITE(sUnit)nbytePhi,(((phi(i,j,k),i=0,nx),j=0,ny),k=0,nz)
WRITE(sUnit)lf//'</AppendedData>'//lf
WRITE(sUnit)'</VTKFile>'//lf
CLOSE(sUnit)

!*************************************************************************************!
! Determine Narrow Band
!*************************************************************************************!

ALLOCATE(phiSB(0:nx,0:ny,0:nz))
ALLOCATE(phiNB(0:nx,0:ny,0:nz))

CALL narrowBand(nx,ny,nz,dx,phi,phiNB,phiSB)

!*************************************************************************************!
! Initialize Gradients
!*************************************************************************************!


ALLOCATE(gradPhi(0:nx,0:ny,0:nz,3))
ALLOCATE(grad2Phi(0:nx,0:ny,0:nz,3))
ALLOCATE(gradMixPhi(0:nx,0:ny,0:nz,3))

gradPhi = 0.
grad2Phi = 0.
gradMixPhi = 0.
gradPhiMag = 0.

phiN = phi

orderUp = 1
order1 = 2
order2 = 2


!*************************************************************************************!
! Min/Max Flow
!*************************************************************************************!

iter = 2000 !20000
h1 = .1
h1 = h1*dx*dx

DO n = 1,iter


   !****************** Explicit Third Order TVD RK Stage 1 ***********************!

   ! Calculate first derivative if it falls within stencil band
   DO i = 0,nx
      DO j = 0,ny
         DO k = 0,nz
            IF (phiSB(i,j,k) == 1) THEN
               CALL firstDeriv(i,j,k,nx,ny,nz,dx,phi,phiX,phiY,phiZ,order1,gMM)   
               gradPhi(i,j,k,1) = phiX
               gradPhi(i,j,k,2) = phiY
               gradPhi(i,j,k,3) = phiZ   
            END IF
         END DO
      END DO
   END DO

   ! Calculate second derivative flow if it is in the narrow band
   DO i = 0,nx
      DO j = 0,ny
         DO k = 0,nz
            IF (phiNB(i,j,k) == 1) THEN
               CALL secondDeriv(i,j,k,nx,ny,nz,dx,phi,gradPhi,phiXX,phiYY,phiZZ,phiXY,phiXZ,phiYZ,order2)
               grad2Phi(i,j,k,1) = phiXX
               grad2Phi(i,j,k,2) = phiYY
               grad2Phi(i,j,k,3) = phiZZ
               gradMixPhi(i,j,k,1) = phiXY
               gradMixPhi(i,j,k,2) = phiXZ
               gradMixPhi(i,j,k,3) = phiYZ          
            END IF

         END DO
      END DO
   END DO

   ! Caluclate the min/max flow
   DO i = 0,nx
      DO j = 0,ny
         DO k = 0,nz
            IF (phiNB(i,j,k) == 1) THEN
               CALL firstDeriv(i,j,k,nx,ny,nz,dx,phi,phiX,phiY,phiZ,order1,gMM)   
               CALL minMax(i,j,k,nx,ny,nz,dx,phi,gradPhi,grad2Phi,gradMixPhi,F,gridX)
               k1 = F !F*gMM !max(F,0.)*dPlus + min(F,0.)*dMinus ! F*gM  
               phi1(i,j,k) = phi(i,j,k) + h1*k1
            END IF

         END DO
      END DO
   END DO


   !****************** Explicit Third Order TVD RK Stage 2 ***********************!

   ! Calculate first derivative if it falls within stencil band
   DO i = 0,nx
      DO j = 0,ny
         DO k = 0,nz
            IF (phiSB(i,j,k) == 1) THEN
               CALL firstDeriv(i,j,k,nx,ny,nz,dx,phi1,phiX,phiY,phiZ,order1,gMM)   
               gradPhi(i,j,k,1) = phiX
               gradPhi(i,j,k,2) = phiY
               gradPhi(i,j,k,3) = phiZ   
            END IF
         END DO
      END DO
   END DO

   ! Calculate second derivative flow if it is in the narrow band
   DO i = 0,nx
      DO j = 0,ny
         DO k = 0,nz
            IF (phiNB(i,j,k) == 1) THEN
               CALL secondDeriv(i,j,k,nx,ny,nz,dx,phi1,gradPhi,phiXX,phiYY,phiZZ,phiXY,phiXZ,phiYZ,order2)
               grad2Phi(i,j,k,1) = phiXX
               grad2Phi(i,j,k,2) = phiYY
               grad2Phi(i,j,k,3) = phiZZ
               gradMixPhi(i,j,k,1) = phiXY
               gradMixPhi(i,j,k,2) = phiXZ
               gradMixPhi(i,j,k,3) = phiYZ          
            END IF

         END DO
      END DO
   END DO

   ! Caluclate the min/max flow
   DO i = 0,nx
      DO j = 0,ny
         DO k = 0,nz
            IF (phiNB(i,j,k) == 1) THEN
               CALL firstDeriv(i,j,k,nx,ny,nz,dx,phi1,phiX,phiY,phiZ,order1,gMM)   
               CALL minMax(i,j,k,nx,ny,nz,dx,phi1,gradPhi,grad2Phi,gradMixPhi,F,gridX)
               k2 = F !F*gMM !max(F,0.)*dPlus + min(F,0.)*dMinus ! F*gM  
               phi2(i,j,k) = 3./4.*phi(i,j,k) + 1./4.*phi1(i,j,k) + 1./4.*h1*k2
            END IF

         END DO
      END DO
   END DO

   !****************** Explicit Third Order TVD RK Stage 3 ***********************!

   ! Calculate first derivative if it falls within stencil band
   DO i = 0,nx
      DO j = 0,ny
         DO k = 0,nz
            IF (phiSB(i,j,k) == 1) THEN
               CALL firstDeriv(i,j,k,nx,ny,nz,dx,phi2,phiX,phiY,phiZ,order1,gMM)   
               gradPhi(i,j,k,1) = phiX
               gradPhi(i,j,k,2) = phiY
               gradPhi(i,j,k,3) = phiZ   
            END IF
         END DO
      END DO
   END DO

   ! Calculate second derivative flow if it is in the narrow band
   DO i = 0,nx
      DO j = 0,ny
         DO k = 0,nz
            IF (phiNB(i,j,k) == 1) THEN
               CALL secondDeriv(i,j,k,nx,ny,nz,dx,phi2,gradPhi,phiXX,phiYY,phiZZ,phiXY,phiXZ,phiYZ,order2)
               grad2Phi(i,j,k,1) = phiXX
               grad2Phi(i,j,k,2) = phiYY
               grad2Phi(i,j,k,3) = phiZZ
               gradMixPhi(i,j,k,1) = phiXY
               gradMixPhi(i,j,k,2) = phiXZ
               gradMixPhi(i,j,k,3) = phiYZ          
            END IF

         END DO
      END DO
   END DO

   ! Caluclate the min/max flow
   DO i = 0,nx
      DO j = 0,ny
         DO k = 0,nz
            IF (phiNB(i,j,k) == 1) THEN
               CALL firstDeriv(i,j,k,nx,ny,nz,dx,phi2,phiX,phiY,phiZ,order1,gMM)   
               CALL minMax(i,j,k,nx,ny,nz,dx,phi2,gradPhi,grad2Phi,gradMixPhi,F,gridX)
               k3 = F !F*gMM !max(F,0.)*dPlus + min(F,0.)*dMinus ! F*gM  
               phiN(i,j,k) = 1./3.*phi(i,j,k) + 2./3.*phi2(i,j,k) + 2./3.*h1*k3 
            END IF

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
   IF (phiErr < 1.E-7) THEN
      PRINT*, " Distance function time integration has reached steady state "
      EXIT
   END IF
   
   ! set phi to new value
   phi = phiN
   
   PRINT*, " Iteration: ",n," ", " RMS Error: ",phiErr
  
   ! check for a NAN
   IF (isnan(phiErr)) STOP

   
   !**************************** Reinitialization *******************************!
 
   ! Currently not used  

   !****************************** Narrow Band **********************************!

   CALL narrowBand(nx,ny,nz,dx,phi,phiNB,phiSB)


END DO
print*,

! print out run time
CALL cpu_time(t4)
PRINT*, " Total Run Time: ",t4-t1," Seconds"
PRINT*,
   

!*************************************************************************************!
! Asymptotic Error
!*************************************************************************************!

phiErr = 0.

! calculate RMS
DO i = 0,nx
   DO j = 0,ny
      DO k = 0,nz
         phiErr = phiErr + (phi(i,j,k)-phiO(i,j,k))*(phi(i,j,k)-phiO(i,j,k))
      END DO
   END DO
END DO

phiErr = sqrt(phiErr/(nx*ny*nz))

PRINT*, " Asymptotic Error: ",phiErr


!*************************************************************************************!
! Output Grad Phi Mag
!*************************************************************************************!


! grad phi
order1 = 2
DO i = 0,nx
   DO j = 0,ny
      DO k = 0,nz
         CALL firstDeriv(i,j,k,nx,ny,nz,dx,phi,phiX,phiY,phiZ,order1,gMM)
         gradPhiMag(i,j,k) = gMM
      END DO
   END DO
END DO



!*************************************************************************************!
! Paraview Output
!*************************************************************************************!


PRINT*, " Writing Out Cartesian Grid to Paraview Format "
PRINT*, 


! output to Paraview
WRITE(extent,'(3(A3,I6))')' 0 ',nx,' 0 ',ny,' 0 ',nz
WRITE(origin,'(3(F20.8,A1))')xLo(1),' ',xLo(2),' ',xLo(3),' '
WRITE(spacing,'(3(F20.8,A1))')dx,' ',dx,' ',dx,' '
nbytePhi =(nx+1)**3*24
offset = 0
WRITE(coffset,'(I16)')offset

sUnit = 11
OPEN(UNIT=sUnit,FILE='blockAfter.vti',FORM='unformatted',ACCESS='stream',STATUS='replace')
WRITE(sUnit)'<?xml version="1.0"?>'//lf
WRITE(sUnit)'<VTKFile type="ImageData" version="0.1" byte_order="LittleEndian">'//lf
WRITE(sUnit)'<ImageData WholeExtent="',TRIM(extent),'" Origin="',TRIM(origin),'" Spacing="',TRIM(spacing),'">'//lf
WRITE(sUnit)'<Piece Extent="',TRIM(extent),'">'//lf
WRITE(sUnit)'<PointData Scalars="phi">'//lf
WRITE(sUnit)'<DataArray type="Float64" Name="phi" format="appended" offset="',TRIM(coffset),'"/>'//lf
WRITE(sUnit)'</PointData>'//lf
WRITE(sUnit)'</Piece>'//lf
WRITE(sUnit)'</ImageData>'//lf
WRITE(sUnit)'<AppendedData encoding="raw">'//lf
WRITE(sUnit)'_'
WRITE(sUnit)nbytePhi,(((phi(i,j,k),i=0,nx),j=0,ny),k=0,nz)
WRITE(sUnit)lf//'</AppendedData>'//lf
WRITE(sUnit)'</VTKFile>'//lf
CLOSE(sUnit)

sUnit = 12
OPEN(UNIT=sUnit,FILE='blockAfterGrad.vti',FORM='unformatted',ACCESS='stream',STATUS='replace')
WRITE(sUnit)'<?xml version="1.0"?>'//lf
WRITE(sUnit)'<VTKFile type="ImageData" version="0.1" byte_order="LittleEndian">'//lf
WRITE(sUnit)'<ImageData WholeExtent="',TRIM(extent),'" Origin="',TRIM(origin),'" Spacing="',TRIM(spacing),'">'//lf
WRITE(sUnit)'<Piece Extent="',TRIM(extent),'">'//lf
WRITE(sUnit)'<PointData Scalars="phi">'//lf
WRITE(sUnit)'<DataArray type="Float64" Name="gradPhiMag" format="appended" offset="',TRIM(coffset),'"/>'//lf
WRITE(sUnit)'</PointData>'//lf
WRITE(sUnit)'</Piece>'//lf
WRITE(sUnit)'</ImageData>'//lf
WRITE(sUnit)'<AppendedData encoding="raw">'//lf
WRITE(sUnit)'_'
WRITE(sUnit)nbytePhi,(((gradPhiMag(i,j,k),i=0,nx),j=0,ny),k=0,nz)
WRITE(sUnit)lf//'</AppendedData>'//lf
WRITE(sUnit)'</VTKFile>'//lf
CLOSE(sUnit)


!*************************************************************************************!
!
! Program End
!
!*************************************************************************************!

END PROGRAM set3d

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
         IF (abs(phi(i,j,k)) < 5.1*dx) THEN
            phiNB(i,j,k) = 1
         END IF

         ! stencil band
         IF (abs(phi(i,j,k)) < 10.1*dx) THEN
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


INTEGER,INTENT(IN) :: i,j,k,order
REAL,INTENT(IN) :: dx
REAL,DIMENSION(0:nx,0:ny,0:nz),INTENT(IN) :: phi
REAL,INTENT(OUT) :: phiX,phiY,phiZ,gMM
REAL :: aa1,aa2,aa3,aa4,aa5,aa6,aa7
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

SUBROUTINE secondDeriv(i,j,k,nx,ny,nz,dx,phi,gradPhi,phiXX,phiYY,phiZZ,phiXY,phiXZ,phiYZ,order)



REAL :: dxx,dx12
INTEGER,INTENT(IN) :: i,j,k,order
REAL,INTENT(IN) :: dx
REAL,DIMENSION(0:nx,0:ny,0:nz),INTENT(IN) :: phi
REAL,DIMENSION(0:nx,0:ny,0:nz,3),INTENT(IN) :: gradPhi
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

SUBROUTINE minMax(i,j,k,nx,ny,nz,dx,phi,gradPhi,grad2Phi,gradMixPhi,F,gridX)

REAL :: dxx,curv,b1,b2,b3,b4,b5,b6,denom,gM,HH
INTEGER :: h
INTEGER,INTENT(IN) :: i,j,k
REAL,INTENT(IN) :: dx
REAL,DIMENSION(0:nx,0:ny,0:nz),INTENT(IN) :: phi
REAL,DIMENSION(0:nx,0:ny,0:nz,3),INTENT(IN) :: gradPhi,grad2Phi,gradMixPhi,gridX
REAL,INTENT(OUT) :: F
       
phiX = gradPhi(i,j,k,1)
phiY = gradPhi(i,j,k,2)
phiZ = gradPhi(i,j,k,3)

phiXX = grad2Phi(i,j,k,1)
phiYY = grad2Phi(i,j,k,2)
phiZZ = grad2Phi(i,j,k,3)

phiXY = gradMixPhi(i,j,k,1)
phiXZ = gradMixPhi(i,j,k,2)
phiYZ = gradMixPhi(i,j,k,3)

! gradient magnitude
gradMag2 = phiX*phiX+phiY*phiY+phiZ*phiZ             
IF (gradMag2 < 1.E-13) THEN
   gM = 0.
ELSE
   gM = sqrt(gradMag2)
END IF  

! calculate curvature
IF (gM < 1.E-13) THEN
   curv = 0.0;
ELSE 
   curv = phiXX+phiYY+phiZZ
END IF


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

SUBROUTINE weno(gM,i,j,k,nx,ny,nz,dx,phi,dPlus,dMinus)

IMPLICIT NONE

INTEGER,INTENT(IN) :: nx,ny,nz,i,j,k
REAL,INTENT(IN) :: dx
REAL,DIMENSION(0:nx,0:ny,0:nz),INTENT(IN) :: phi
REAL,INTENT(OUT):: gM,dPlus,dMinus
REAL :: gradX,gradY,gradZ,pa,pb,pc,pd,pe,pf,na,nb,nc,nd,ne,nf,a,b,c,d,e,f,uzm1,uzm2,uzm3
REAL :: u,ap,am,bp,bm,cp,cm,dp,dm,IS0p,IS0m,IS1p,IS1m,IS2p,IS2m,epsp,epsm
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

dPlus = max(a,0.)**2+min(b,0.)**2+max(c,0.)**2+min(d,0.)**2+max(e,0.)**2+min(f,0.)**2
dPlus = sqrt(dPlus)

dMinus = min(a,0.)**2+max(b,0.)**2+min(c,0.)**2+max(d,0.)**2+min(e,0.)**2+max(f,0.)**2
dMinus = sqrt(dMinus)


END SUBROUTINE weno 



!*************************************************************************************!
! End of Code 
!*************************************************************************************!

