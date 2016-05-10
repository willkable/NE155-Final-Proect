# NE155-Final-Proect
2D Diffusion Equation Solver


For this program I have created a python tkinter GUI that can take inputs and solve the fixed source 2D diffusion equation.
To run the program simply import all the files in this repository into a folder. Then you can simply run the python file
and it should work.

Instructions:
1. Simply Run the 2D Diffusion Solver.py file and a window should appear.
2. Input the mesh size.
3. Input your upper a lower boundary positions for the problem. (A default set should already be entered)
4. Input your source information. If the source is uniform throughout the system simply check the box and input the
   source strength. 
5. If the source is not uniform, input the source strength and specify the source location in mesh coordinates. 
   (i.e. Don’t input Cartesian coordinates of the System) The mesh is of size (nxn) and the source is located at some mesh position (i, j))
   Where i,j are integers in the range(n). You may simply enter these mesh coordinates in the form 'x1,y1' in the box labeled 
   'Source Mesh Location'. An example should already be in the box.
6. If there are multiple sources you may enter them as well. In the box labeled 'Source(s)' enter your source strengths in
   the form 'S1,S2,....Sn' for as many sources as you have. Again an example should already be in the box. To places these
   Sources on the mesh, do the same as you did in part 5. However now you can input in the form 'x1, y1,x2,y2,....,xn,yn'
7. For the macroscopic cross sections simply input some non-zero/non-negative values for the absorption and transport cross section
8. You may also input a tolerance to be used in the Gauss-Seidel iteration method. The default is set at 1e-7.
9. Once all information is entered you may create a Surface/Heat/Side Plot of the Flux Solution. You may also export the flux solution
   to a file you specify under the 'Export Data' button. A default 'Out.txt' is already specified.

The 3 test output files included were create I created with the system conditions in the 'Test_Inputs.txt' file.


