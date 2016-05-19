# levelSetFortran
Level set methods code. Reads an .stl file and creates a level set field using WENO5. Also runs min/max flow on the level set field to smooth the geometry using 2nd order central differencing. 

I have included two sample .stl files to get you started. One is a cube, and the other is two cubes with a unit spacing of 10 between them. 

The code will output 3 .vti files: blockBefore, which is the initialized level set field, blockAfter, which is the level set field after min/max has been run, and blockAfterGrad, which is just for debugging purposes. I suggest ParaView to open these files.

Compile with:

gfortran -fdefault-real-8 -O3 set3d.f90 -o set3d.exec

Couple of notes:
- Check the size of the geometry inside the .stl file. Ensure your dx is small or large enough to have at least 10 cells inside the zero level set. 
- The WENO5 is pretty sensitive to time step changes. 
- Currently has no capability to do moving geometry, that is in the works.
- Currently in serial. Parallel version is in the works.
- Currently uses uniform grid spacing. 
- The .stl files need to be clean. No pierced faces or dirty CAD.
- Current code is pretty messy. I will try and clean it up later when I have time.
