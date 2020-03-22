Introduction
============

.. currentmodule:: ifstype

.. contents::
   :local:

Installation
------------

TODO: currently not hosted on PyPI ... when this is done, installation will be easy.


Basic Usage
-----------

The first thing do to is to define an IFS.
Here is an example::

   from ifstype import AffineFunc
   from ifstype.exact import Fraction
   
   a = Fraction(1,5)
   b = Fraction(1,8))
   
   ifs = IFS([AffineFunc(a,0),
              AffineFunc(b,a-a*b),
              AffineFunc(a,1-a-b+a*b),
              AffineFunc(b,1-b)])

An :class:`ifstype.AffineFunc` instance represents a real affine function.
For example, ``AffineFunc(r,d)`` is equivalent to the affine contraction function ``f(x)=r*x+d``.

For basic usage, the module provides two general functions: :func:`ifstype.verify_wft` and :func:`ifstype.run_ifs`.
Called without any keyword arguments, :func:`ifstype.run_ifs` takes two parameters: an IFS instance, and a folder name::

   from ifstype import run_ifs

   tr_graph = run_ifs(ifs, "output")

This will generate a transition graph instance, and also print the graph visualization to the file ``output/graph.pdf`` and IFS information to ``output/info.txt``.
It is often also useful to visualize the interval structure of the IFS: as a result, :func:`ifstype.run_ifs` also accepts the boolean keyword argument ``with_gens=True``, which will print the generations visualization to the file ``output/gens.pdf``.
This can be really slow if the IFS has a lot of net intervals!

Note that this information is less informative (interesting) if the IFS does not satisfy the finite neighbour condition.
Thus another general function :func:`verify_fnc` is a streamlined version of :func:`run_ifs` which generates no output, but will determine whether or not the IFS does in fact satisfy the finite neighbour condition.
To use this, for example::

   from ifstype import verify_wft

   verify_wft(ifs)

Finally, it is often useful to generate diagrams highlighting the complete behaviour of the generations of the IFS within the first k generations.
We thus have the function :func:`run_ifs_gens`, which takes three parameters: an IFS instance, a folder name, and a depth value.
For example::
   
   from ifstype import run_ifs_gens

   run_ifs_gens(ifs,"output",3)

The depth value indicates how many children to compute beneath each net interval, starting from the root net interval [0,1].

