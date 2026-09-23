.. _public_api_policy:

Public API and compatibility policy
===================================

PyFiberModes has one stable public import namespace: the names listed in
``PyFiberModes.__all__`` and imported directly from ``PyFiberModes``. For
example::

   from PyFiberModes import Fiber, LP01, MaterialRegistry, SolverResult

These names are covered by the compatibility policy below. Their documented
signatures, return semantics, and import locations are treated as public API.
The public surface is stored in ``tests/api_snapshot.json`` and checked in CI.
An intentional API change must update that snapshot, the changelog, and the
:doc:`migration_guide` in the same pull request.

Internal modules
----------------

Imports from implementation modules are not stable unless the object is also
exported at the package root. In particular, ``PyFiberModes.solver.*``,
``PyFiberModes.services``, and private names beginning with an underscore are
internal machinery. They may change between minor versions. Advanced users
may inspect or extend them, but should isolate those dependencies in their own
adapter layer.

Deprecation lifecycle
---------------------

When a supported name or behavior must change, PyFiberModes will:

#. introduce the replacement and document it in the migration guide;
#. keep the former API working while emitting ``DeprecationWarning``;
#. retain it for at least two minor releases and six months, whichever is
   longer; and
#. remove it in the next eligible major release, except when an urgent safety
   or correctness issue makes continued support impractical.

Deprecation messages identify the replacement and earliest removal version.
Patch releases do not intentionally break supported imports or behavior.
Minor releases may add public API and deprecate existing API. Major releases
may remove deprecated API and make documented incompatible changes.

Experimental features
---------------------

An API is experimental only when its documentation labels it as such.
Experimental APIs are excluded from the compatibility guarantee and use an
``experimental`` module or an explicit warning in their documentation. An
unlabelled root export is always considered stable.

Downstream recommendations
--------------------------

Applications should import supported objects from ``PyFiberModes`` instead of
their defining modules, test with deprecation warnings enabled, and constrain
major versions for reproducible deployments. Solver plugins should implement
the published protocols while avoiding inheritance from internal backends.
