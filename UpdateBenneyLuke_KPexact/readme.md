
Update Benney-Luke demonstration for Firedrake website (Onno: works 15-02-2021 on MacBook):
- bennylukefb.py reproduces Figure 3 in Bokhove and Kalogirou 2016 [To do: add that Figure 3 as well for comparison]
- limited output to times shown in Figure
- output of energy in energy.txt; made Python program energyplotBL.py to plot these data and reproduce lower panel of Figure 3. [To do: add that new output figure on energy as well for comparison]
- Paraview: Open, apply choose eta and eta exact; then under Filter -> Data Analysis -> Plot over line; choose cross-section along mid channel (tick use normal and XY plane before overlaying a line); to do: lines for different times in Plot over line? To do: put in correct "exact solution" for reflection?

Update exact KP-solution 9A_analytical.py (does not work anymore 22-02-2021; fix in progress):
- less output; smaller domain; added Petsc include to find maximum value; printed in file.
- Paraview of 3D exact solution: Open, Filters: alphabetical -> warp by scalar; View: colormap/bar editor; then top-right icons placing of bar etc; range of colorbar: to right of colorbar range: 2nd box on right from top; go in help and type axes (Alphabetical -> Axes) then find adjustments.

Under firedrake and firedrake-new updated bennylukefb.py and 9A_analytical.py work fine (27-04-2020 output generated).





