.. Environ documentation api methods file.
   Created by Edan Bainglass on Sun Jul 3 2022.

Methods
=======

``init_interface``
##################
    - used to initialize the interface

``init_io``
###########
    - used to initialize Environ's I/O
    - parameters
        - ``ionode`` (``LOGICAL``)
            - if the I/O node processor
        - ``ionode_id`` (``INTEGER``)
            - the I/O node processor number
        - ``comm`` (``INTEGER``)
            - the MPI communicator
        - ``program_unit`` (``INTEGER``)
            - the host's stdout unit
        - ``lstdout`` (``LOGICAL``)
            - if Environ is allowed to write to stdout

``read_input``
##############
    - used to read Environ's input file
    - parameters
        - ``filename`` (``CHAR``) optional
            - the name of Environ's input file
            - defaults to ``environ.in``
        - ``nsx`` (``*INTEGER``) optional
            - maximum number of species
            - defaults to 10
            - used for species array allocations

``destroy``
###########
    - used to destroy Environ's objects
    - parameters
        - ``level`` (``INTEGER``) optional
            - destruction scheme (calculation-dependent)
            - defaults to 1

``update_cell``
###############
    - used to set Environ's simulation cells at the start of a calculation
    - note that Environ does not support cell updating during the calculation
    - parameters
        - ``at`` (``REAL(DP)`` | ``SHAPE(3, 3)``)
            - the cell lattice

``update_electrons``
####################
    - used to update Environ's electrons at each SCF step
    - parameters
        - ``rho_in`` (``REAL(DP)`` | ``SHAPE(:)``)
            - the electronic density array
        - ``nelec`` (``INTEGER``) optional
            - the expected number of electrons
            - if present, used as a sanity check against the integrated density
        - ``lscatter`` (``LOGICAL``) optional
            - if the density should be scattered over multiple processors
            - defaults to ``.FALSE.``

``get_verbosity``
#################
    - ``FUNCTION``
    - getter for Environ's verbosity level

``calc_potential``
##################
    - used to compute Environ's potential contribution at each SCF step
    - parameters
        - ``update`` (``LOGICAL``)
            - if Environ should start computing its potential
        - ``potential`` (``READ(DP)`` | ``SHAPE(:)``)
            - the returned calculated potential
        - ``local_verbose`` (``INTEGER``) optional
            - if present, added to Environ's verbosity level
        - ``lgather`` (``LOGICAL``) optional
            - if the potential should be gathered from multiple processors
            - defaults to ``.FALSE.``
