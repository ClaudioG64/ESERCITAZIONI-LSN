# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Repository Overview

This is a collection of C++ and Jupyter notebook exercises for the Numerical Simulation Laboratory (NSL) course from the Physics Department at Università degli Studi di Milano. The repository contains exercises focused on Monte Carlo methods, statistical analysis, and numerical simulation techniques.

## Code Architecture

### Random Number Generation Framework
- **Core Component**: `random.h` and `random.cpp` - RANNYU algorithm-based random number generator used across all exercises
- **Initialization Pattern**: Uses `Primes` file for parameters and `seed.in` for random seed configuration
- **Usage**: All C++ exercises depend on this shared random number generation system

### Exercise Structure
Each exercise directory (`Esercitazione X`) contains:
- **C++ Implementation**: Main simulation code (`es1_1.cpp`, `main.cpp`, etc.)
- **Jupyter Notebook**: Analysis and visualization (`LSN_Exercises_XX.ipynb`)  
- **Configuration Files**: `seed.in`, `Primes` for random number setup
- **Output Data**: `.dat` files with simulation results

### Advanced Simulator (Exercise 4+)
The NSL_SIMULATOR is a comprehensive Monte Carlo simulation framework:
- **System Class**: Core simulation engine (`system.h`, `system.cpp`) supporting multiple simulation types
- **Particle Class**: Particle representation with position/velocity management
- **Architecture**: Modular design with separate INPUT, OUTPUT, SOURCE, and CONFIG directories
- **Dependencies**: Uses Armadillo C++ library for linear algebra operations

## Common Development Commands

### Building C++ Code

**Basic Exercises (1-3)**:
```bash
# In exercise directory
make
# Or for specific exercise
make es1_1
```

**NSL Simulator (Exercise 4+)**:
```bash
cd "Esercitazione X/X_Y/NSL_SIMULATOR/SOURCE"
make
# Creates simulator.exe
```

**Cleaning Build Files**:
```bash
make clean
# For NSL simulator, also removes output files:
make remove
```

### Running Simulations

**Basic Exercises**:
```bash
./es1_1  # or ./main.exe
```

**NSL Simulator**:
```bash
cd "NSL_SIMULATOR/SOURCE"
./simulator.exe
# Results written to ../OUTPUT/
```

## File Naming Conventions

- Exercise source files: `es1_1.cpp`, `es1_2.cpp`, etc.
- Main programs: `main.cpp`
- Output data: `output1_1.dat`, `option_prices.dat`, etc.
- Configuration: `input.dat`, `properties.dat`
- Random seed files: `seed.in` (input), `seed.out` (output)

## Development Notes

### Compilation Settings
- Standard: C++11 minimum, C++14 for advanced exercises
- Optimization: `-O3` for performance-critical simulations
- Armadillo library integration for matrix operations in NSL_SIMULATOR

### Data Analysis Workflow
1. Run C++ simulation to generate `.dat` files
2. Use Jupyter notebooks for data visualization and analysis
3. Statistical uncertainty estimation via blocking method is standard practice

### Random Number Generation
All simulations use the custom RANNYU algorithm. Proper initialization requires:
- `Primes` file with two prime numbers
- `seed.in` with 4-element seed array
- Call to `rnd.SaveSeed()` to preserve state