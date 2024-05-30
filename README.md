# Weighted-blow-up-of-Monomial-ideals
There is both a notebook, and a python file with code for doing Weighted blow-up of Monomial ideals as described in https://arxiv.org/pdf/1906.07106 by DAN ABRAMOVICH, MICHAEL TEMKIN, AND JAROSLAW WLODARCZYK.

The algorithm used is however constructed concurrently with my Master's thesis (Jonas Pedersen), and it is using methods I could not find anyone else describing, but probably similar to "Newton polyhedra without coordinates" by Boris Youssin.
There is no proof as of yet that the invariant produced in this algorithm corresponds to one due to ABRAMOVICH et al. however it seems very likely that these coincide.

Examples of use of the python file:
Change 
```
if __name__ == "__main__":
    out = multiblowup_parallel(make_monomial(minors2d(4,MM)))
    with open("Output.txt", "w") as text_file:
        text_file.write(out)
    print(out)
```
To for example 
```
if __name__ == "__main__":
    out = multiblowup_parallel("b22*b33, b33*t, b22**2 * t, b23*t, b12*b22*t, b22* t**2, t**3")
    with open("Output.txt", "w") as text_file:
        text_file.write(out)
    print(out)
```
It is also possible to run the code with cmd arguments. For example
```
.../speciale-parallel.py outfile --stop_depth 10 ---ideal t**2*x**2,t**3*y
```

In general the code can only blow up at the origin, but it can handle monomial ideals written as a string using "*" between two variables that should be multiplied. Variables can be on any form of letter plus numbers, i.e. "b32" or "x1". No underscores.
There is also a notebook, which is easier to play with, but it does not support multicore processing as the python file does.
