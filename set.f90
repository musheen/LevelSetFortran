PROGRAM set

IMPLICIT NONE

INTEGER :: n,i,j,nx,sUnit,nbytePhi,offset,nxs
REAL :: c,a0,a1,a2,a3,a4,t,b,xLo(3),xHi(3),dx,x(2),xs,ys,yp,ypp,g,gp,x0,x1,dxs
REAL,ALLOCATABLE,DIMENSION(:,:) :: phi
CHARACTER(len=1) :: lf=char(10)
CHARACTER(len=1024) :: extent,origin,spacing,coffset


! define 0012 airfoil coefficients
c  = 1.
t  = 0.12
b  = 5.*t*c
a0 = b*0.2969/SQRT(c)
a1 =-b*0.1260/c
a2 =-b*0.3516/(c*c)
a3 = b*0.2843/(c*c*c)
a4 =-b*0.1036/(c*c*c*c)



! define Cartesian grid and allocate phi
nx  = 256
xLo =(/-1.,-1.5,0./)
xHi =(/2.,1.5,0./)
dx  =(xHi(1)-xLo(1))/REAL(nx)
ALLOCATE(phi(0:nx,0:nx))
phi = 1.D14
 

! find phi at each Cartesian grid point
x0  = 0.
x1  = 1.
nxs = 10000
dxs =(x1-x0)/REAL(nxs)
DO j=0,nx
   x(2) = ABS(xLo(2)+REAL(j)*dx)
DO i=0,nx
   x(1) = xLo(1)+REAL(i)*dx
search:DO n=0,nxs
   xs       = x0+REAL(n)*dxs
   ys       = a0*SQRT(xs)+a1*xs+a2*xs*xs+a3*xs*xs*xs+a4*xs*xs*xs*xs
   phi(i,j) = MIN(phi(i,j),SQRT((x(1)-xs)*(x(1)-xs)+(x(2)-ys)*(x(2)-ys)))
END DO search
END DO
END DO

DO j=0,nx
   ys = ABS(xLo(2)+REAL(j)*dx)
DO i=0,nx
   xs = xLo(1)+REAL(i)*dx
   IF (xs >= 0. .AND. xs <= 1.) THEN
      IF (ys < a0*SQRT(xs)+a1*xs+a2*xs*xs+a3*xs*xs*xs+a4*xs*xs*xs*xs) &
           phi(i,j) =-phi(i,j)
   END IF
END DO
END DO


! output to Paraview
WRITE(extent,'(3(A3,I6))')' 0 ',nx,' 0 ',nx,' 0 ',0
WRITE(origin,'(3(F20.8,A1))')xLo(1),' ',xLo(2),' ',xLo(3),' '
WRITE(spacing,'(3(F20.8,A1))')dx,' ',dx,' ',0.,' '
nbytePhi =(nx+1)**2*SIZEOF(c)
offset = 0
WRITE(coffset,'(I16)')offset

sUnit = 11
OPEN(UNIT=sUnit,FILE='block.vti',FORM='unformatted',ACCESS='stream',STATUS='replace')
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
WRITE(sUnit)nbytePhi,((phi(i,j),i=0,nx),j=0,nx)
WRITE(sUnit)lf//'</AppendedData>'//lf
WRITE(sUnit)'</VTKFile>'//lf
CLOSE(sUnit)


END PROGRAM set
