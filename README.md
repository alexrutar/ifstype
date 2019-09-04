# IFSType

This library is a simulation library iterated function systems of contractive similitudes on the real line.
See (TODO: add link) for a description of the theory behind this project.
The primary purpose of the library is to generate the transition graph of such an IFS, and to provide tools for working with and analyzing the transition graph.

Main features:
- Exact precision rational and algebraic numbers.
- Symbolic computation of transition matrices.
- Visualization tools: transition graph drawing, and IFS visualization through intervals.

## Basic Usage
The first thing do to is to define an IFS.
Here is an example:

```python
from ifstype import CtrFunc
from ifstype.exact import Rational

a=Rational(1,5)
b=Rational(1,8))

ifs = IFS([CtrFunc(a,0),
           CtrFunc(b,a-a*b),
           CtrFunc(a,1-a-b+a*b),
           CtrFunc(b,1-b)])

```
Here we have an example of a family if IFS, parametrized by values `a` and `b`.

A CtrFunc instance represents a contraction function `CtrFunc(r,d)` is equivalent to the affine contraction function `f(x)=r*x+d`.

The decorator `ifs_family` is a convenience used to streamline the creation of an IFS instance.
We then this function to generate a specific IFS.
This library also provides two general functions: `verify_wft` and `run_ifs`.
Called without any keyword arguments, `run_ifs` takes two parameters: a foldername, and an IFS instance.
```python
from ifstype import run_ifs

tr_graph = run_ifs(ifs, "output")
```
This will generate a transition graph instance, and also print the graph visualization to the file "output/graph.pdf" and IFS information to "output/info.txt".
Note that this information is less informative (interesting) if the IFS is not weak finite type.
Thus the other general function `verify_wft` is a streamlined version of `run_ifs` which generates no output, but will determine whether or not the IFS is in fact weak finite type.
To use this, for example,
```python
from ifstype import verify_wft

verify_wft(ifs)
```
