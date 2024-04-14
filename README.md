# Finite-Element-Method-Project
The program capable of calculating the non-stationary solution of heat transfer problem using Finite Element Method.

## NOTICES
- The project has been translated to english, so the polish PDF including code screenshots might differ.
- The project already includes 4 grid files: Test1_4_4.txt, Test2_4_4_MixGrid.txt, Test3_31_31_square.txt and Test4_31_31_trapezium.txt. Any other grid files should work as well but shall be placed in proper directory.

## About the program

### FEM introduction

  The Finite Element Method (FEM) is a technique for solving boundary problems encountered in physics and engineering. It involves transforming any continuous quantity into its discrete model. The discrete model is characterized by a limited number of nodes, which in turn define a limited number of finite elements. In practice, this means dividing a given area into smaller elements (these sub-areas are the aforementioned finite elements), making it easier to conduct calculations and achieve the most accurate result possible.

  The implemented program performs calculations for thermal problems. The FEM algorithm for such problems takes into account a limited number of points. These points are the aforementioned nodes, which form the mesh of the model (the mesh of finite elements). The values to be determined at each node are temperature values. This stems from the nature of the chosen boundary problem. In the case of the presented program, these are thermal problems, but for other tasks, they could be values of a different function. The area should also be divided into a limited number of sub-areas. Finite elements with common nodes approximate the shape of the computational area.

  The algorithm also involves interpolating the temperature within each finite element using a polynomial, which can be determined using the temperature values at the nodes. However, it's important to maintain the continuity condition of temperature at the boundaries of elements and to select the polynomial for each element accordingly. It's also crucial to select the nodal temperatures in such a way as to ensure the best approximation of the temperature field compared to reality. For this purpose, Fourier's equation is utilized, which is a differential equation of heat conduction. The number of Fourier equations needed to obtain a result depends on the number of unknown nodal temperature values and must be equal to it.

### The idea of FEM

  The solution to our problem of unsteady heat exchange in a two-dimensional system involves solving a system of equations given by the formula:
([H]+[C]/Δτ){t1} - ([C]/Δτ){t0} + {P} = 0
  By solving the system of equations, we aim to determine {t1}, which represents the nodal temperatures after a time interval ∆𝜏. To calculate these temperatures, we will also need:
- the matrix [H] - it describes the flow of heat through individual nodes of the grid,
- the matrix [C] - it describes the nodes' capacity for storing thermal energy,
- the vector {t0} - representing the nodal temperatures at time 𝜏 = 0,
- the vector {P} - describing the influence of ambient temperature on individual nodes.

### Purpose of the program and what it consists of

  The program is meant to calculate the results of thermal problems using FEM algorithm. The obtained result will allow to picture how heat spreads across a given surface over time. Given result can be used for Paraview program.
  My project is a multi-file program consisting of "main.cpp" along with two other source files and their corresponding header files. "Main.cpp" contains basic operations like reading the grid file and storing its values into appropriate structures, as well as function calls. 
  
  Files "Dispaying.cpp" and "Displaying.h" are responsible for displaying the results. There are 3 available options - displaying the calculations (the results of calculated H matrixes, results of aggregations etc.), displaying the results of heat diffusion (calculated for each step time) and displaying the results of calculations for each step time (so it can be used for Paraview visualisation).
  
  Files "Calculations.cpp" and "Calculations.h" are the core of the program. They contain all the functions and structures responsible of calculating the results. The structures I used to store information about the grids loaded from the file are GlobalData and Element.
  
  GlobalData stores general information about the grid, including the simulation time, time step value, thermal conductivity coefficient value, heat exchange coefficient value, ambient temperature, starting temperature, material density value, specific heat, number of nodes, and the number of elements.
  
  The Element structure, on the other hand, is responsible for storing information about the node ID numbers of a specific element. They are stored in a double vector named ID. Vectors x and y are responsible for storing the coordinates x and y of the loaded nodes respectively, and the vector BC_xy stores information about which nodes have boundary conditions. There are also many functions that make the program working, all of them are described in the program files. 
