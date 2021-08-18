# patches

    Shell scripts for patching QE with Environ.

### Changes in v2.0

---

- Redesigned patching to work with (decoupled) Environ installation scheme
  - No longer calling QE/install/addsonpatch script
  - Adds references to Environ/src and Environ/libs in target/Makefile
  - Modifies QE/install/makedeps.sh to add Environ/src dependency for target package
- Added QE-version check in PW patch to account for nspin dependency in QE < 6.4