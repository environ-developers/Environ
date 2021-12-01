<font size="7">STYLE GUIDE</font>

<br>

# Keywords

- uppercase all Fortran built-in keywords

# Indentation

- 4 spaces
- blocks

```
IF (condition) THEN
    code
ELSE IF
    code
END IF
```

- wrapped lines

```
IF (condition) &
    code
```

# ! spacing (before and after a line)

- return/exit/stop
- implicit none
- single if-statements, except when only line in a block
- do/select/if/etc. blocks (including nested blocks)

```
!
code...
!
IF (condition) do something
!
code...
!
some block
  IF (condition) do something
end block
!
```

- for select blocks, group case (with code) in separate blocks

```
SELECT CASE (x)
    !
CASE (1)
    code
    !
CASE (2)
    code
    !
END SELECT
```

- wrapped lines (when using & continuation)
- for clarity (attach comment)
- between IN, OUT/INOUT, and local-variable blocks in subroutines

```
!
........INTENT(IN)
........INTENT(IN)
........INTENT(IN), OPTIONAL
!
........INTENT(INOUT)
........INTENT(OUT)
........INTENT(OUT), OPTIONAL
!
local vars
!
```

- CALL statements, except when only line in a block

```
!
CALL something
!
CALL something
!
some block
    CALL something
end block
!
```

<br>

# COMMENTS

# Attached to single-line code (lowercase)

- < col 90 (including comment):
  - on the same line, separated by a single blank space from code

```
code ! comment...
```

- \>= col 90:
  - directly below code (\*above code if line is continued across multiple lines)

```
code
! comment...

! comment
code...... &
    ......
```

# Section comments (sentence-case)

- see `!---------- sectioning`

# !---------- sectioning

- up to, but not including col 90
- use to section code blocks (multiple commands of common purpose)

```
!
!------------------------------------------------------------------------------------
! Header comment describing the following few lines of code
!
code...
code...
code...
!
!------------------------------------------------------------------------------------
!
```

- use dashed line to separate USE statements, variable declarations, and body
- section-closing dashed line not required if final block
  - for IF, DO, SELECT blocks, final means before END statement

# Docstrings

- docstrings should be added to the following blocks:
  - module
  - interface
  - type definition
  - subroutine
  - function
- the following template should be used (supports doxygen automated documentation):

```
    !------------------------------------------------------------------------------------
    !>
    !! docstring...
    !! ...
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE name(param)
        !--------------------------------------------------------------------------------
        !
        code
```

- leave !> line and last !! line blank

<br>

# OPERATORS

# Conditional operators

- use `> >= == \= <= <` instead of `.LT. .LE. .EQ. .NE. .GE. .GT.`
- use `.EQV. .NEQV.` when comparing logical values, e.g. `(1 < 2) .EQV. (3 \= 4)`

# Operator spacing

- +-_/ single space on both ends ----> ```1 _ 2 + 2 / 4 - 3```
  - no space around `-` if used as a negative sign ----> `-7`
- no space around powers ----> `3**2`
- single space on both ends of all conditionals ----> `a == b .AND. c < d`

# Short-circuit

- DO NOT assume that FORTRAN will apply short-circuit logic

<br>

# SELECT TYPE/CASE

- Use a default case to catch unintended behavior
  - SELECT TYPE -> CLASS DEFAULT
  - SELECT CASE -> CASE DEFAULT
- Exception: omit default case if the code should simply proceed (ignore other cases)

```
SELECT CASE (deriv_method)
!
CASE IS ('fft')
    ...
    !
CASE IS ('chain')
    ...
    !
END SELECT
```

- In the above case, any other option is meant to be ignored

<br>

# VARIABLE DECLARATION

- define same-type variables in a single (or continued) line using a single TYPE descriptor
  - exceptions:
    - if initialized
    - for clarity - must be accompanied by a comment (see `TYPE2`)

```
TYPE1 :: name1, name2, name3 &
         name4, name5, name6...
!
TYPE2 :: name7 ! comment
TYPE2 :: name8 ! comment
TYPE3 :: name9, name10
```

- use `TYPE, ATTRIBUTE, INTENT :: var` order for argument declarations
  - `ATTRIBUTE` is alphabetized, e.g. (O)PTIONAL, (T)ARGET, etc.

# Dimensions

- `TYPE :: name(dim)`
  - when defining a single variable
- `TYPE :: name(dim1), name(dim2)`
  - when defining multiple variables of the same TYPE, but different dimensions
- `TYPE, DIMENSION(dim) :: name1, name2`
  - when defining multiple variables of the same TYPE and dimensions, or

<br>

# Continuation-line &

- ! spacing on both ends of wrapped block
- use 4-space tabs to align

```
IF (condition) &
    code...
```

- try to keep products/quotients together

```
some_variable = var1 + var2 - &
                var3 * var4 / var5
```

- use prefixed & for strings (required) and to improve expression readability

```
some_string = "blah blah blah &
              &blah blah blah"
```

```
sum_ = 1 + var1 + &
       & + var2 * var3
```

<br>

# subroutines/functions

- always use IMPLICIT-NONE
- separate declarations from body (see `!---------- sectioning`)
- group INTENT(IN) and INTENT(OUT) separately, with OPTIONALs at the bottom of each block
- match order of declaration with order of parameters as best/close as possible

<br>

# Imports

- unless importing a class, use explicit imports (DO NOT rely on nested imports!!!)
- grouped at the top of the module, !-separated into subgroups (in order):
  - utils imports
  - FFTs imports
  - constants
  - representation classes (cell, density, etc.)
  - core classes (containers, then cores; no !-separator)
  - solver classes (setup, then solvers; no !-separator)
  - physical-quantity classes
  - I/O
  - debugging
- each block is sorted alphabetically

<br>

# Private/public

- use empty private (to prevent import inheritance)
- explicitly define all that is publicly accessible

<br>

# Error/info/warning messages

- use `io%error` to raise an error
  - requires sub_name and error code
  - sentence-cased message
- use `io%write` for indented messages
- use `io%warning` for non-indented, non-terminating warning messages
  - outputs `'ENVIRON WARNING (ierr): message'`
  - lowercased message

<br>

# ASSOCIATE BLOCKS

- use ASSOCIATE blocks instead of many POINTERs

```
ASSOCIATE (item1 => this%item1, &
           item2 => this%item2)
```

- pointers can point to expressions

```
ASSOCIATE (item1 => SUM(some_items))
```

- pointers cannot reference one another

  incorrect

  ```
  ASSOCIATE (item1 => this%item1, &
             attr1 => item1%some_attr)
  ```

  correct

  ```
  ASSOCIATE (item1 => this%item1, &
             attr1 => this%item1%some_attr)
  ```

- ASSOCIATE block pointers do not require declaration
- ASSOCIATE block targets do not require the TARGET attribute

<br>

# GENERAL SYNTAX

- avoid GOTO statements

- DO NOT exceed col 90

  - exceptions for long strings, but not past col 130

- use `array = ...` when operating on the entire array (no need for `array(:)`)

  - `array(1, :) = ...` is necessary when operating on all members along specific dimensions

- use `i`, `j`, `k`, `l`, `m`, `n` for loop indices whenever possible

- do not use `TRIM(ADJUSTL())` for input parameter checks
  - input parameters are pre-tested against set values in `checkin` routines
