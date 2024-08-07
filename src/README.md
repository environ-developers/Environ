# Environ Source Code

    Primary codebase for Environ.

### Changes in v3.0

---

- Added API for use by host software
- Made atom species input more flexible (now by label, number, or weight)
- Redesigned boundaries in OOP hierarchal structure
- Improved performance when using soft spheres (ionic boundary)
- Added ms-gcs functionality

### Changes in v2.0

---

- Added double-cell and Mott-Schottky functionalities
- Forces now computed with a Gaussian-spread description of ions
- Redesigned modules as classes with type-bound procedures
- Redesigned solvers, cores, and functions using OOP inheritance
- Packaged calculation routines in a calculator class
- Packaged cleaning routines in a destructor class
- Packaged simulation parameters and numerical base in setup class
- Isolated environment objects (including potentials) in environment class
