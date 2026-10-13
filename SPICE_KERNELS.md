# Native Swift SPICE kernels

## Supported first increment

`SpiceKernelSet` is a caller-owned, synchronized text kernel pool. It loads
numeric/string assignments in identified `KPL/PCK`, `KPL/FK`, `KPL/LSK` and
`KPL/MK` files. These identifiers classify files; they do not yet enable physical
orientation, frame or leap-second interpretation. Binary DAF/SPK and other
architectures produce explicit unsupported-format errors.

```swift
let kernels = try SpiceKernelSet()
let load = try kernels.load(URL(fileURLWithPath: "/assets/solar-system.tm"))
let view = kernels.snapshot()
if case .numbers(let radii)? = view["BODY399_RADII"]?.value {
    // Raw PCK radii are conventionally kilometres, not library SI states.
    let radiiMetres = radii.map { $0 * 1_000 }
}
try kernels.unload(load)
```

Supported syntax includes comment/data blocks, case-sensitive names of at most
32 ASCII characters, direct and incremental assignments, parenthesized vectors,
commas or whitespace, finite decimal numbers with E/D exponents, and single-quoted
strings with doubled quotes. String elements are limited to 80 characters and
cannot cross lines. Numeric and string elements cannot be mixed within a vector.
An incremental assignment appends to a preceding value of the same type; changing
the type replaces that value. Raw string arrays retain continuation markers.
Empty vectors, non-ASCII text, non-finite numbers and `@` date literals are rejected.
In particular, ordinary LSK files containing date literals are not yet loadable.
The parser is a documented subset, not a claim of complete toolkit compatibility.

Meta-kernels support `KERNELS_TO_LOAD`, `PATH_SYMBOLS`, `PATH_VALUES` and `+` string
continuation. Symbols in this subset contain ASCII letters, digits or underscores.
Undefined or duplicate symbols, mismatched arrays, dangling continuations, missing
files and dependency cycles produce errors. Substitutions are single-pass, without
shell/environment expansion. Nested manifests are supported with bounded depth;
all manifest variables remain local to loading. Canonical paths detect symlink
cycles. Dependency paths are relative to the containing manifest rather than the
process working directory. Spaces in local paths are allowed. Those last policies
are intentional library conventions rather than full CSPICE path compatibility.

## Ownership, precedence and invalidation

Every successful load creates an opaque ownership group, including its manifest
and dependencies in traversal order. Later groups have higher priority. Within a
file, later direct assignments replace earlier values and incremental assignments
append. Inspection returns contributing source paths and immediate parent manifests.
Repeated loading creates another group, taking a fresh snapshot of file contents.
Unloading by load ID removes exactly that group's contributions, rebuilding the
pool from retained earlier groups without rereading files. Shared dependencies
loaded separately remain present. Physical files are never deleted.

Loads commit atomically. A parse, dependency or capacity failure leaves the previous
set and cache scope intact. This differs deliberately from toolkit partial loading.
Each successful load/unload/clear creates a new `EphemerisCacheScope`. External
caches must namespace by that scope; there are currently no internal state caches.
Immutable snapshots retain their contents and scope after mutation. Do not combine
values from different snapshots when constructing a configured provider.

An internal lock serializes mutation and snapshot creation, including file reads.
Evaluation of immutable, `Sendable` snapshots requires no lock. There is no global
pool and no callback-based reader. Concurrent load order follows lock acquisition;
applications needing reproducible precedence should issue loads sequentially.
Caller-selected limits bound file bytes, total loaded file occurrences, dependency
depth, assignments and values. Unknown unload IDs are errors.

## Specification basis and validation

This implementation was independently written from the public
[NAIF Kernel Required Reading](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/req/kernel.html),
revision **2021 DEC 29**, consulted 2026-10-13: Kernel Type Identification;
Text Kernel Specifications (names, assignments, values, additional syntax);
Additional Meta-kernel Specifications; and Kernel Management (priority, path
symbols, relative paths, reload, unload). No NASA/NAIF implementation source was
consulted or translated, and there is no CSPICE dependency.

Swift Testing fixtures are synthesized from that specification. They cover
replacement/append precedence, quoted strings, block comments, manifest traversal,
continued paths, provenance, independent sets, reload/unload, retained snapshots,
transactional rollback, missing dependencies, recursive manifests, unsupported
formats and bounded parsing. They require neither large assets nor CSPICE.

## Remaining work

Next: validated DAF inspection and bounded random-access reads, followed by SPK
2/3 interpolation, centre chains and geometric SI ephemeris adapters. Date literals
and LSK, text PCK/FK interpretation remain separate increments. Binary PCK, CK/SCLK,
file-backed fallback integration and independent black-box reference fixtures follow
as tracked in `ROADMAP.md`. Kernel data provenance, licensing and redistribution
must be assessed per supplied asset; loading does not grant redistribution rights.
