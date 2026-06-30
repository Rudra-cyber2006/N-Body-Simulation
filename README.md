# High-Performance N-Body Galaxy Simulation (Barnes–Hut)

This repository presents a high-performance two-dimensional N-body gravitational simulation developed in C++. The system models a galaxy consisting of over 5,000 stars orbiting a central supermassive black hole, with a strong emphasis on computational efficiency and physical accuracy.

To address the inherent scalability challenges of gravitational simulations, the project implements the Barnes–Hut algorithm, reducing computational complexity from (O(N^2)) to (O(N \log N)). This enables efficient large-scale simulations while maintaining reasonable accuracy. The engine is further supported by custom memory management strategies, robust numerical integration techniques, and Python-based visualization tools for analysis and rendering.

*(Note: Insert the `star_1_high_res.gif` file here to showcase the simulation output.)*

---

## Key Features and Performance Enhancements

Designing a large-scale physics simulation requires careful consideration of both numerical stability and hardware efficiency. This project incorporates several optimizations to achieve high performance:

**Barnes–Hut Tree ((O(N \log N)))**
The simulation space is recursively subdivided into quadrants using a hierarchical tree structure. Distant clusters of stars that satisfy a threshold criterion ((\theta = 0.5)) are approximated as a single mass, significantly reducing the number of force computations per time step.

**Custom Memory Allocation**
To eliminate the overhead associated with frequent dynamic memory allocation (`new`/`delete`), the engine employs a pre-allocated contiguous memory pool of 500,000 nodes. This approach improves cache locality and minimizes interaction with the operating system’s memory manager.

**Velocity Verlet Integration**
The simulation uses the symplectic Velocity Verlet integration scheme instead of the Forward Euler method. This ensures improved long-term stability and better conservation of energy, which is critical for accurately modeling orbital dynamics.

**Optimized Mathematical Operations**
Distance comparisons are performed using squared values, thereby avoiding unnecessary square root computations. The expensive `std::sqrt` operation is invoked only when strictly required, reducing computational overhead.

**Recursion Control Mechanisms**
To maintain robustness, the implementation includes recursion depth limits and skips zero-mass nodes. These safeguards prevent potential computational bottlenecks and ensure stability in edge cases such as particle overlap.

**Proximity Constraints**
To avoid numerical singularities and extreme velocity spikes, particles are prevented from initializing too close to the central black hole. This maintains stability in force calculations and preserves the integrity of the simulation domain.

---

## Project Structure

* **`Barnes_Hut.cpp`**
  Core implementation of the optimized C++ physics engine.

* **`phase2.cpp`**
  An earlier version of the simulation, retained for comparison between brute-force and tree-based approaches.

* **`orbits.csv`**
  Output dataset containing positional data of particles over time.

* **`visuals.ipynb`**
  Jupyter Notebook for plotting orbital trajectories and analyzing galactic structure using Python, Pandas, and Matplotlib.

* **`anim.ipynb`**
  Notebook for generating high-resolution animations (GIFs) from simulation data.

---

## Build and Execution Instructions

### Compile the Simulation Engine

To achieve optimal performance, compile the code with high-level compiler optimizations enabled:

```bash
g++ -O3 Barnes_Hut.cpp -o sim
```
