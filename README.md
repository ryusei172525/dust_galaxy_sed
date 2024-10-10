# Galaxy SED Library

A Python library for calculating galaxy SEDs using C++ for internal computations. Users can specify parameter ranges like maximum galaxy age and mass, and the library uses MCMC to optimize parameters.

## Project Structure

```
galaxy_sed/
├── README.md
├── setup.py
├── galaxy_sed/
│   ├── __init__.py
│   ├── sed_model.py       # Python wrapper for SED calculations
│   ├── mcmc.py            # MCMC implementation
│   ├── utils.py           # Utility functions
│   ├── cpp/
│   │   ├── sed_calculator.cpp  # SED calculations in C++
│   │   ├── sed_calculator.hpp  # Header file
│   │   ├── CMakeLists.txt  # Build configuration for C++ code
│   └── bindings/
│       ├── pybind_wrapper.cpp  # Wrapper using pybind11
│       └── CMakeLists.txt  # Build configuration for pybind11
```

## Required Libraries

- **pybind11**: Used for interfacing between Python and C++.
- **NumPy**: Supports numerical calculations.
- **scipy**: Used for MCMC sampling.
- **cmake**: Needed for building the C++ components.

## Getting Started

1. Build the C++ components:

   ```bash
   mkdir galaxy_sed/bindings/build
   cd galaxy_sed/bindings/build
   cmake ..
   make
   ```

2. If the `make` command completes successfully, the library is ready to use, and you can run the `example.ipynb` file to see the functionality.

## Usage

After building, you can use the library to perform SED calculations and MCMC optimization directly from Python. More detailed examples can be found in `example.ipynb`.
