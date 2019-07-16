# levelSetFortran

Three dimensional level set methods code. Reads a .stl file and creates a signed distance function field using WENO5. Also runs min/max flow on the level set field to smooth the geometry using 2nd order central differencing. 

I have included two sample .stl files to get you started. One is a cube, and the other is two cubes with a unit spacing of 10 between them. 

The code will output 2 .vti files: signedDistanceFunction.vti, which is the initialized level set field, smoothedDistanceFucntion.vti, which is the level set field after min/max has been run. I suggest ParaView to open these files. 

The code will also output a .s3d file. This is a mesh file type for the CFD code Strand3dFC. 

Compile using the make file. Should be simple to follow. Working on adding a namelist for inputs.

Couple of notes:
- The code is pretty sensitve to the time step. 
- Check the size of the geometry inside the .stl file. Ensure your dx is small or large enough to have at least 10 cells inside the zero level set. 
- Currently has no capability to do moving geometry.
- Currently in serial. Parallel version is in the works.
- Currently uses uniform grid spacing. 
- If you want to only run the signed distance function part, set the min/max iterations to zero. 
- The .stl files need to be clean. No pierced faces or dirty CAD.
- Currently uses .stl only right now.

