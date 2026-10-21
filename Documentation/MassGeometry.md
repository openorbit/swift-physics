# Rigid mass components

`MassComponent` is a value tree, independent of scene nodes and ephemeris frames.
Each component defines its own mass, centre and intrinsic inertia, a rigid pose
relative to its parent, and child components. Structural nodes may have no mass.
Identifiers are unique throughout an assembly. Intrinsic inertia excludes children.
Cylinder helpers use the local X axis as the cylinder axis.

`massGeometry()` rotates each intrinsic tensor and applies the parallel-axis
theorem about the assembled centre. The result is expressed in the root parent's
axes and origin, including the root's own placement.

`detaching(id:)` returns retained and detached value trees atomically. The detached
root carries the composed parent transforms, so both results remain in the same
structural coordinate system. Recompute both geometries, then call
`offsetKinematics(from:angularVelocity:)` for each result. Rotate these offsets
into the simulation reference frame and add the original assembly position and
velocity. Both bodies initially inherit the original orientation and angular
velocity; their subsequent dynamics use their respective inertia tensors.

Separation springs and motors require explicit impulses with matching recoil.
This API does not integrate flexible structures, slosh, articulated joints,
propellant momentum transport, or time-dependent internal motion. Fuel depletion
can rebuild the value tree at a step boundary; instantaneous reference-point
kinematics alone is not a complete variable-mass dynamics model.
