# Clausing Factor Monte Carlo Simulation in Rust

This project is a Rust implementation of a Monte Carlo simulation for calculating the **Clausing factor** for Charge Exchange (CEX). It was converted from the original VBA (Visual Basic for Applications) spreadsheet code.

## Overview

The Clausing factor is used in vacuum physics and ion propulsion systems to characterize the conductance of cylindrical apertures and grid systems. This Monte Carlo simulation tracks particles through a dual-cylinder geometry to calculate:

- **Clausing Factor**: The transmission probability through the grid system
- **Downstream Correction Factor**: A velocity-based correction factor
- **Max Count**: Maximum number of collisions for any particle
- **Particles Lost**: Number of particles that exceeded the iteration limit

## Code Structure

The implementation consists of:

- **`ClausingParams`**: Input parameters structure
  - `thick_screen`: Thickness of the screen grid
  - `thick_accel`: Thickness of the acceleration grid
  - `r_screen`: Radius of the screen grid
  - `r_accel`: Radius of the acceleration grid
  - `grid_space`: Space between grids
  - `npart`: Number of particles to simulate

- **`ClausingResults`**: Output results structure
  - `clausing_factor`: The calculated Clausing factor
  - `max_count`: Maximum collision count
  - `nlost`: Number of lost particles
  - `den_cor`: Downstream correction factor

- **`clausing()`**: Main Monte Carlo simulation function

## Building and Running

### Prerequisites

- Rust (1.93.0 or later)
- Cargo (comes with Rust)

### Build

```bash
cargo build --release
```

### Run

```bash
cargo run --release
```

Or run the compiled binary directly:

```bash
./target/release/clausing
```

### Test

```bash
cargo test
```

## Example Output

```
Running Clausing factor calculation...
Parameters: ClausingParams { thick_screen: 1.0, thick_accel: 0.5, r_screen: 2.0, r_accel: 1.0, grid_space: 0.3, npart: 10000 }

Results:
  Clausing Factor: 0.775200
  Max Count: 12
  Particles Lost: 0
  Downstream Correction Factor: 0.937916
```

## Algorithm Description

The simulation uses a Monte Carlo approach:

1. **Particle Launch**: Particles are launched from the bottom of the screen grid with random positions and velocities following a cosine distribution
2. **Ray Tracing**: Each particle is traced through the geometry, calculating intersections with cylinder walls
3. **Diffuse Reflection**: When a particle hits a wall, it is re-emitted with a new random velocity direction
4. **Exit Conditions**: Particles either escape through the top, return through the bottom, or exceed the iteration limit
5. **Statistical Analysis**: The Clausing factor is calculated from the fraction of particles that escape

## Key Differences from VBA

- **Type Safety**: Rust's strong type system prevents many common errors
- **Memory Safety**: No risk of buffer overflows or memory leaks
- **Performance**: Compiled Rust code is significantly faster than interpreted VBA
- **Modern Syntax**: Cleaner, more maintainable code structure
- **Testing**: Built-in unit testing framework

## Customization

To use different parameters, modify the `main()` function:

```rust
let params = ClausingParams {
    thick_screen: 1.0,
    thick_accel: 0.5,
    r_screen: 2.0,
    r_accel: 1.0,
    grid_space: 0.3,
    npart: 10000,  // Increase for better accuracy
};
```

## Dependencies

- `rand = "0.8"`: Random number generation for Monte Carlo simulation

## License

This code is provided as-is for educational and research purposes.
