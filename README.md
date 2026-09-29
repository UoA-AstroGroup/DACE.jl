# DACE.jl

DACE.jl is a Julia interface to the [DACE](https://github.com/dacelib/dace) library.

Documentation can be found [here](https://UoA-AstroGroup.github.io/DACE.jl/).

## Example usage

DACE.jl can be installed using the Julia package manager. From the Julia REPL,
type `]` to enter the Pkg REPL mode and run

```
pkg> add DACE
```

Run a simple example:

```
$ julia examples/sine.jl
x
     I  COEFFICIENT              ORDER EXPONENTS
     1    1.0000000000000000e+00   1   1
------------------------------------------------
y = sin(x)
     I  COEFFICIENT              ORDER EXPONENTS
     1    1.0000000000000000e+00   1   1
     2   -1.6666666666666666e-01   3   3
     3    8.3333333333333332e-03   5   5
     4   -1.9841269841269841e-04   7   7
     5    2.7557319223985893e-06   9   9
     6   -2.5052108385441720e-08  11  11
     7    1.6059043836821616e-10  13  13
     8   -7.6471637318198174e-13  15  15
     9    2.8114572543455210e-15  17  17
    10   -8.2206352466243310e-18  19  19
------------------------------------------------
y(1.0)   = 0.8414709848078965
sin(1.0) = 0.8414709848078965
```

## Running the tests

More examples of how to use *DACE.jl* can be found in the tests, for example [validation_2.jl](test/validation_2.jl).

You can run the validation tests using the Julia package manager.
Clone this git repository and change into the cloned directory.
Then enter the Julia REPL using `julia --project` and
type `]` to enter the Pkg REPL mode and instantiate the package:

```
(DACE) pkg> instantiate
```

Then run the tests:

```
(DACE) pkg> test
```
