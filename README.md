# Finite Element Method
The program calculates temperatures on a 2D grid generated in Abaqus for a steady-state heat flow process and simulates temperature changes over time during unsteady-state heat flow. The behavior of the plate is obtained by solving a system of FEM equations derived from Fourier's equation.

([H] + [C]/Δτ) * {t₁} - ([C]/Δτ){t₀} + [P] = 0

H - Stiffness matrix (thermal conductivity matrix, with boundary condition)

C - Heat capacity matrix

Δτ - Current time interval

t1 - Temperature vector at the next time step

t0 - Temperature vector at the previous time step

P - Load vector (representing heat sources, boundary condition)
<p align="center">
  <img src="assets/animation.gif" width="600"><br>
  <em>Heat simulation over time for 31x31 square grid</em>
</p>

After simulation is finished, ready gif file with simulation will save in folder assets, and plot with figure after end of simulation will render.

## Install dependencies
You should create and activate new virtual environment, then in project folder use command:
```bash
pip install -r req.txt
```

## Run
Use the grids located in the test_grid_data folder. You can also upload your own grid file there. For better precision use 4 integration points schema.
```bash
python -u main.py <grid_file_name.txt> <number_of_integration_points(2 or 3 or 4)>
```
