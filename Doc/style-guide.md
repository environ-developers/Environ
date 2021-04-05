# Style Guide

## keywords

- uppercase all Fortran built-in keywords

## indentation

- 4 spaces
- if/do blocks

```
IF (condition) THEN
    code
ELSE IF
    code
END IF

DO x = 1, 10
    code
END DO
```

- select blocks

```
SELECT CASE (x)
CASE (1)
    code
CASE (2)
    code
END SELECT
```

- wrapped lines

```
IF (condition) &
    code
```

## ! spacing (before and after a line)

- return/exit/stop
- implicit none
- do/select/if/etc. blocks (including nested blocks)
- single if-statements, except when only line in a block

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

- wrapped lines (when using & continuation)
- for clarity (attach comment)
- between IN, OUT/INOUT, and local-variable blocks in subroutines

```
!
........(IN)
........(IN)
!
........(OUT)
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

## comments

### attached to single-line code

- <= col 90 (including comment):
  - on the same line, separated by a single blank space from code

```
code ! comment...
```

- col 91+:
  - directly below code (\*above code if line is continued across multiple lines)

```
code
! comment...

! comment
code...... &
    ......
```

### !---------- sectioning

- up to, but not including col 90
- use dashed line to separate variable declaration from subroutine/function body
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

- use dashed line to separate variable declaration from subroutine/function body
- section-closing dashed line not required if final block
  - for IF, DO, SELECT blocks, final means before END statement

## docstrings

- docstrings should be added to the following blocks:
  - module
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

- leave !> line and last line blank

## operator spacing

- +-_/ single space on both ends ----> ```1 _ 2 + 2 / 4 - 3```
  - no space around `-` if used as a negative sign ----> `-7`
- no space around powers ----> `3**2`
- single space on both ends of all conditionals ----> `a .EQ. b .AND. c .LT. d`

## arguments/parameters

- keep intent(OUT) and intent(INOUT) variables at the end of the list

## variable definition

- define same-type variables in a single (or continued) line using a single TYPE descriptor
  - exceptions:
    - if initialized
    - for clarity - must be accompanied by a comment (see `TYPE2`)

```
TYPE1 :: name1, name2, name3 &
         name4, name5, name6...

TYPE2 :: name7 ! comment
TYPE2 :: name8 ! comment
TYPE3 :: name9, name10
```

### dimensions

- `TYPE :: name(dim)`
  - when defining a single variable
- `TYPE :: name(dim1), name(dim2)`
  - when defining multiple variables of the same TYPE, but different dimensions
- `TYPE, DIMENSION(dim) :: name1, name2`
  - when defining multiple variables of the same TYPE and dimensions, or

## continuation-line &

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
