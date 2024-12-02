# Zero Pressure Gradient Boundary Layer

# Overview
## Project Goals and Benefits

This project aims to deepen the understanding of boundary layer formation by analyzing key parameters along a flat plate. By employing numerical methods, it seeks to solve the Navier-Stokes equations—five fully coupled, time-dependent 3D partial differential equations—leveraging Prandtl's theoretical framework to simplify the problem. This approach provides insights into the steady-state evolution of boundary layer dynamics, which is crucial for optimizing designs in engineering fields like aerodynamics.

The investigation focuses on the evolution of various boundary layer parameters in a planar, steady, incompressible boundary layer over a flat plate with a zero-pressure gradient, as described by boundary layer theory. Specifically, the analysis will address:

* Velocity Profile
* Boundary Layer Thickness
* Displacement Thickness
* Wall Friction


# Scientific Details
## Introduction

The project builds upon Ludwig Prandtl's foundational boundary layer theory from 1904, simplifying the Navier-Stokes equations by dividing the flow domain into two distinct regions: the outer flow region, where inviscid assumptions can be applied, and the thin boundary layer near the surface, where viscous effects are significant. 

Within the boundary layer, the fluid motion is governed by the Prandtl boundary layer equations, which describe the velocity profile and other relevant parameters. These equations are:

Continuity Equation:
    
$$
    \frac{\partial u}{\partial x} + \frac{\partial v}{\partial y} = 0  
$$

Momentum Equation:
    
$$
    \frac{\partial uu}{\partial x} + \frac{\partial uv}{\partial y} = -\frac{\partial p}{\partial x} + \frac{\partial \tau}{\partial y} 
$$

Shear Stress Equation:
    
$$
    \tau = \frac{1}{Re} \frac{\partial u}{\partial y}
$$

Where:

* u and v are the velocity components in the x (streamwise) and y (wall-normal) directions, respectively.
* p is the pressure field.
* τ is the shear stress.
* Re is the Reynolds number.

### Boundary Conditions
To solve the boundary layer equations, we must specify appropriate boundary conditions for the velocity components at different locations within the flow domain. These boundary conditions are:

#### Initial conditions (x = 0)

The velocity component in the x-direction, u, is given by an initial velocity profile

$$ u(0, y) = u_0(y) $$

The velocity component in the y-direction, v, is also specified by an initial profile

$$ v(0, y) = v_0(y) $$

#### At the surface of the plate (y = 0)

The velocities in the boundary layer at the wall are zero due to the no-slip condition.

$$ u(x, 0) = v(x, 0) = 0 $$

#### Far away from the wall (as y → ∞)

The velocity in the x-direction u approaches the free-stream velocity.

$$ u(x, y\to\infty) = u_e(x) $$

These boundary conditions ensure the system is well-posed for numerical solution.

### Physical parameters

The physical parameters used for the simulation are:

    Reynolds number (Re) = 1000,
    Plate length (L) = 10,
    Boundary layer height (H) = 2,

For the initial velocity profiles, we use the following functions:
 
$$ u_0(y) = tanh⁡(5y) $$

$$ v_0(y) = 0 $$

Additionally, the free-stream velocity, is constant due to the zero-pressure gradient:

$$ u_e(x) = 1 $$

All variables are dimensionless, and the Reynolds number is based on the free-stream velocity at x = 0.

## Numerical Scheme

The numerical method employed to solve the governing equations utilizes a uniform grid in both the x-direction (along the length of the flat plate) and the y-direction (normal to the surface). The flow domain is discretized as follows:

* The x-coordinate spans the range [0, L] and is divided into N=100 intervals.
* The y-coordinate spans the range [0, H] and is divided into M=600 intervals.

The resulting grid provides a systematic representation of the flow domain, where each grid point corresponds to specific values of the flow variables.

    y
    ^
    |
    |     H +---+---+-/ /-+---+---+
    |       |   |   |     |   |   | 
    |       +---+---+-/ /-+---+---+
    |       |   |   |     |   |   |
    |
    |       |   |   |     |   |   |
    |       +---+---+-/ /-+---+---+
    |       |   |   |     |   |   |
    |     0 +---+---+-/ /-+---+---+
    |       0                     L
    +------------------------------------>
                    x

To solve the boundary layer equations, the Keller Box Algorithm, a finite-difference method, is used. This algorithm discretizes the flow domain into a grid of cells and applies a central differencing scheme. The values of key variables (such as the velocity components u and v) at the center of each cell are computed through linear interpolation from the surrounding grid points. The method uses finite difference approximations to convert the partial differential equations into algebraic equations, which are solved iteratively. The algorithm proceeds by iterating through the grid, progressively refining the solution until convergence is achieved.

## Simulation Data

The generated simulation data is stored in a structured text format in the file simulation-data.csv. The file is divided into two main sections:

### First Section

This section contains general information about the flow along the streamwise coordinate (in the x-direction). For each x-value (ranging from 0 to N), the following parameters are recorded:

    'Um':     The mean streamwise velocity.
    'Vm':     The mean normal velocity.
    'd':      The displacement thickness.
    'd99':    The boundary layer thickness at which the velocity reaches 99% of the free-stream velocity.
    'tw':     The wall shear stress.

This section is formatted as follows:

    x Um Vm d d99 tw
    0 0.929966 0.000000 0.138634 0.529339 0.005000
    1 0.928717 0.045704 0.141137 0.532872 0.004728
    ...

### Second Section

This section contains detailed data for each point in the grid, including both the x- and y-coordinates. For each x-coordinate, the following is recorded:

    'y':          The vertical grid point.
    'u':          The streamwise velocity at the grid point (x,y).
    'v':          The normal velocity at the grid point (x,y).
    'tau_xy':     The shear stress at the grid point (x,y).

This section is formatted as follows:

    x 0
    y u v tau_xy
    0 0.000000 0.000000 0.005000
    1 0.016665 0.000000 0.004999
    ...
    #
    x 1
    y u v tau_xy
    0 0.000000 0.000000 0.004728
    1 0.015761 0.000030 0.004729
    ...
    #
    x 2
    y u v tau_xy
    0 0.000000 0.000000 0.004632
    1 0.015440 -0.000019 0.004631

### Extracting the Data

To extract the simulation data, the following pseudo code can be used:

```plaintext
    (1) Open the file raw_simulation_data.csv for reading.

    (2) Initialize variables to store the extracted data:
        'x_values':     A list to store the grid index x-coordinates (from the first section).
        'Um_values':    A list to store the streamwise velocity values (Um) from the first section.
        'Vm_values':    A list to store the normal velocity values (Vm) from the first section.
        'd_values':     A list to store the displacement thickness values from the first section.
        'd99_values':   A list to store the boundary layer thickness at 99% velocity values from the first section.
        'tw_values':    A list to store the wall shear stress values from the first section.
        'flow_data':    A dictionary to store detailed flow data, indexed by grid x and y: 
                                flow_data[x] = {y: [u, v, tau_xy]}

    (3) While reading the file: 
        a. Read each line. 
        b. If the line starts with "x Um Vm d d99 tw":
            This marks the start of the first section. Skip this line and continue to the next. 
        c. If the line contains a data row like "x Um Vm d d99 tw":
            Parse the values and append them to the corresponding lists: 'x_values', 'Um_values', 'Vm_values', 'd_values', 'd99_values', and 'tw_values'. 
        d. If the line contains "x <value>":
            Parse the grid index x-coordinate from the line and initialize a new entry in flow_data for this grid index. 
        e. If the line contains "y u v tau_xy":
            This marks the start of the detailed grid data. Skip this line and continue to the next. 
        f. If the line contains data like "i u v tau_xy":
            Parse the grid index y value and the corresponding values for u, v, and tau_xy​.
            Store these values in 'flow_data[x][y]'.

    (4) Important Transformation:
        When processing the raw simulation data, the grid points (x and y) represent discrete indices. To convert these into physical coordinates, use the following formulas:

            xtrue = x / N * L
            ytrue = y / M * H

        Where N and M are the grid sizes, and L and H are the physical dimensions of the domain. The converted values correspond to the actual coordinates in the flow domain.

    (5) Close the file.

    The extracted data is now organized as follows:
        * 'x_values', 'Um_values', 'Vm_values', 'd_values', 'd99_values', and 'tw_values' contain the general information about the flow along the streamwise coordiante.

        * 'flow_data' is a dictionary where each grid index x is a key, and the corresponding values for u, v, and τxy​ at each grid index y are stored.
```

This structure provides easy access to both the overall flow characteristics and the detailed distribution of velocity and shear stress, with the transformation to physical coordinates prepared for further analysis or visualization.

# Technical Details
## Technical Requirements

    Programming Language: Python [Version>=3.12.2]

    Libraries: NumPy [Version>=1.26.4], Matplotlib [Version>=3.8.0] (for numerical methods and plotting)

    Tools: Any text editor or IDE (e.g., VS Code, PyCharm)


## Usage

To run the simulation, execute the main.py script. The script discretizes the grid, formulates the governing equations into an algebraic system, and initializes the boundary layer solver. It then solves the system iteratively and generates the raw_simulation_data.csv file, which contains the raw simulation data for further analysis.

#### Set up the virtual environment:

```bash
    python -m venv env
```
#### Activate the environment

(a) On Windows:

```bash
    .\env\Scripts\activate
```

(b) On maxOS/Linux:

```bash
    source env/bin/activate
```

#### Install dependencies
```bash
    pip install -r requirements.txt
```

#### Run the script
```bash
    python src/main.py
```

#### Deactivate the environment (when done):
```bash
    deactivate
```

## Folder Structure

The folder structure is designed to keep the project organized, with the source code in the src directory, data files in the data folder, and essential project information such as the license, README, and dependencies in the root directory.

```bash
    /ZPGBL
    ├── src
    │   ├── main.py
    ├── data
    │   ├── simulation-data.csv
    ├── LICENSE
    ├── README.md
    ├── requirements.txt
```

# General Information
## Authors

Kofler Joshua
- [GitHub](https://www.github.com/joshuakofler)
- [ORCID](https://orcid.org/0009-0008-7710-5392)

## Contributing

If you wish to contribute to this project, feel free to fork the repository and submit pull requests. For bug reports or feature requests, please open an issue on the GitHub repository.

## FAQ or Troubleshooting

For any issues or questions, please refer to the GitHub repository and check the issues section or open a new issue.

## References

* Spurk, Aksel. Strömungslehre. Springer, 2016. This book provides a short detailed explanation of the basics of the boundary layer theory.

## License

This project is licensed under the [MIT License](https://choosealicense.com/licenses/mit/) - see the LICENSE file for details.