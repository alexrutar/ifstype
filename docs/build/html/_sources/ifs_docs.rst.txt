.. automodule:: ifstype.ifs

.. autoclass:: ifstype.AffineFunc

   .. autoattribute:: r
   .. autoattribute:: d
   .. automethod:: fixed_point
   .. automethod:: interval
   .. automethod:: __call__
   .. automethod:: compose
   .. automethod:: inverse

.. autoclass:: ifstype.IFS

   .. autoattribute:: funcs
   .. autoattribute:: probs
   .. automethod:: __init__
   .. automethod:: __str__
   .. automethod:: extend
   .. automethod:: invariant_convex_hull

.. autodecorator:: ifstype.ifs_family

.. autoclass:: ifstype.ifs.Neighbour
   :show-inheritance:

   .. autoattribute:: a
   .. autoattribute:: L
   .. automethod:: from_aff

.. autoclass:: ifstype.ifs.NeighbourSet

   .. autoattribute:: neighbours
   .. autoattribute:: lmax
   .. automethod:: __iter__
   .. automethod:: __str__
   .. automethod:: __contains__

.. autoclass:: ifstype.ifs.NetInterval
   :show-inheritance:

   .. automethod:: from_funcs
   .. automethod:: transition_gen
   .. automethod:: normalization_func
   .. automethod:: containing_funcs

.. autoclass:: ifstype.ifs.TransitionMatrix
   :show-inheritance:

   .. automethod:: pos_row
   .. automethod:: spectral_radius
