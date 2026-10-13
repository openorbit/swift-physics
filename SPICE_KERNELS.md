# Native Swift SPICE kernels

## Supported kernel loading

`SpiceKernelSet` is a caller-owned, synchronized text kernel pool. It loads
numeric/string assignments in identified `KPL/PCK`, `KPL/FK`, `KPL/LSK`, `KPL/SCLK` and
`KPL/MK` files. Explicit snapshot adapters supply leap-second tables, sourced
constants, PCK orientations, fixed FK frame providers and type-1 clock correlation. Modern IEEE DAF/SPK and DAF/PCK
files can also be loaded for inspection and bounded array reads.
Primitive geometric state evaluation supports SPK types 2 and 3; other architectures
are rejected by the kernel set.

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
Empty vectors, non-ASCII text and non-finite numbers are rejected.
Date literals support `@YYYY-MON-D` or `@D-MON-YYYY`, with optional clock
components introduced by `-`, `/` or `//`. Hours may stand alone; minutes and
seconds are optional, and seconds may have a fraction. Four-digit years
0001–9999, abbreviated case-insensitive months and Gregorian validation are used. They become unscaled nominal seconds from J2000 noon in the numeric
pool, without UTC/TDB conversion. Other date spellings are diagnosed as unsupported.
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
and dependencies in traversal order. SPK readers appear in `snapshot.binaryKernels`. Later groups have higher priority. Within a
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
Caller-selected limits bound text file bytes, total loaded file occurrences, dependency
depth, assignments and values. `SpiceDAFLimits` separately bounds binary metadata
and per-query array reads; binary assets are not copied into memory in full. Unknown unload IDs are errors.

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

Next: CK loading and attitude providers using the explicit SCLK clock. Additional binary PCK types,
file-backed fallback integration and independent black-box reference fixtures follow
as tracked in `ROADMAP.md`. Kernel data provenance, licensing and redistribution
must be assessed per supplied asset; loading does not grant redistribution rights.

## DAF binary inspection

`SpiceDAFFile` opens modern `DAF/SPK`, `DAF/CK` and `DAF/PCK` files marked
`LTL-IEEE` or `BIG-IEEE`, with an intact published FTP integrity string. Legacy
`NAIF/DAF`, absent format markers, VAX representations, DAS and transfer formats
are unsupported. `SpiceKernelSet` registers SPKs and binary PCKs; low-level CK
inspection does not imply attitude/orientation evaluation support.

Metadata includes internal name, format, size, summary dimensions, ordered arrays
and raw double/integer summary components. SPK descriptors additionally expose NAIF
target and centre IDs, native frame ID, segment type and inclusive coverage in TDB
seconds past J2000. Body centres and barycentres retain distinct integer identifiers.
Every positive segment type remains visible, including those not evaluable.
Inspection does not select or evaluate a segment, validate its type-specific payload,
or convert units/time/frames. Primitive type 2/3 evaluation is a separate layer below.

Array data is read on demand through `readWords(array:offset:count:)`. Array indices
and offsets are zero based; metadata addresses remain one based, inclusive DAF
word addresses. Reads may cross physical record boundaries and cannot escape the
selected array or configured read limit. Raw payload doubles are preserved, including
non-finite values; the future type-specific reader must validate coefficient content.
File seek/read operations are synchronized, allowing concurrent queries on one
reader. Readers close when their final reference is released, so old kernel snapshots
retain usable file handles after unload. Binary assets must remain unchanged while
readers exist; reload refreshes metadata but does not snapshot payload bytes. In-place
asset replacement therefore requires discarding existing views/readers first. Truncated
reads fail explicitly, but equal-size external edits cannot be reliably detected.

The initial profile requires complete 1024-byte physical records and an initial summary
at or before subsequent linked summaries. It checks header dimensions, pointers,
FTP integrity, exact integral summary controls, backward links, tail agreement,
summary/name collisions, finite coverage, array bounds and array/structural overlaps.
Empty summary records are allowed. Limits bound summary traversal, array count,
reserved bytes and each data request. No coefficient caching is introduced yet.

Records between the header and initial summary are exposed verbatim as
`metadata.reservedRecords`, including padding. The public DAF and SPC required
readings describe the comment area and its interfaces but do not specify the
byte-level line/end delimiters. This documentation gap is explicit: the loader
preserves those bytes and does not claim to decode comment lines. Neither the
absence of decoded comments nor unsupported segment types is silently treated as
missing ephemeris data.

Implementation sources consulted 2026-10-13:

- [DAF Required Reading](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/req/daf.html),
  revision 2017 MAR 22: Format Concept, Array Addresses, Organization (File Record,
  Comment Area, Summary Records, Name Records, Element Records), Format Extension.
  FTP bytes come from the published file-record octal example.
- [SPK Required Reading](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/req/spk.html),
  revision 2021 AUG 31: SPK File Structure, Segments and Segment Order/Priority.
- [SPC Required Reading](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/req/spc.html),
  revision 2004 DEC 26: description of the comment area and extraction interfaces.

Synthetic Swift Testing fixtures independently encode headers, multiple linked
summaries, names, reserved bytes and element records in both byte orders. Tests
cover descriptors and unsupported types, cross-record reads, concurrent access,
malformed pointers/coverage/ranges, integrity damage, resource limits, truncation,
kernel ownership and transactional manifest failure. No toolkit implementation
source, toolkit writer or runtime was used.

## SPK type 2/3 primitive states

```swift
let view = kernels.snapshot()
let state = try view.spkState(targetID: 399, atSecondsTDB: 0)
// Relative to state.descriptor.centreID, in state.descriptor.frameID axes.
let positionKilometres = state.positionKilometres
let velocityKilometresPerTDBSecond = state.velocityKilometresPerTDBSecond
```

`SpiceSPKSegment(file:array:)` validates a caller-selected type 2/3 segment's
four-word trailing directory: initial TDB seconds past J2000 (`INIT`), positive
interval length (`INTLEN`), integral record size (`RSIZE`) and positive integral
record count (`N`). The payload must consist of exactly `N * RSIZE` words plus
the directory. A type 2 record contains midpoint, radius and three equal-length
coefficient sets; type 3 contains six sets ordered X, Y, Z, VX, VY, VZ. Degree zero
is supported. Record size is constrained by the DAF reader's configured word limit.
Descriptor coverage may be a subinterval of the coefficient grid.

Evaluation reads only the selected record, validates finite coefficients and a
positive radius consistent with the directory, maps time to [-1,1], and evaluates
all supplied coefficients with Clenshaw summation. The constant coefficient is
used in full. Type 2 analytically differentiates the recurrence and divides by
record radius; type 3 independently evaluates the supplied velocity polynomials.
It does not replace those velocities with derivatives of the position fit.
Overflow or unrepresentable numerical results fail explicitly.

The library's deterministic boundary rule selects the later record at a shared
internal boundary and the last record at the final endpoint. Quotient rounding
near a boundary is checked against the actual directory grid. Descriptor bounds
are inclusive and enforced without broadening the query coverage. Up to eight
ULPs of the larger time/interval scale are allowed for stored record/directory
consistency; normalization clamps only this storage-rounding allowance at the
coefficient endpoints. This is a numerical validation tolerance, not a physical
accuracy bound. Pathologically unresolvable time grids are rejected. No attempt
is made to impose continuity between separately fitted records.

Snapshot target queries search files in reverse load order and arrays in reverse
DAF order, selecting by target and descriptor coverage. They return that primitive
state in its native centre/frame. Unsupported selected types, malformed selected
directories/records and I/O failures propagate, even if earlier data could answer.
Segments that do not cover the epoch do not mask earlier data. The returned state
reports source path, segment name, array index, record index and descriptor. A
missing covering target produces `SpiceSPKError.noSegment`; invalid epochs and
direct out-of-coverage queries have distinct errors. The primitive API does not resolve centres or convert axes/time/units. The
relative-state overload below resolves centres; the SI provider further below
converts units, coordinate-time rates and explicit frame bindings.

Retain a `SpiceSPKSegment` for repeated queries to reuse its validated directory.
Convenience file/snapshot queries validate the selected directory on each call.
No record or result cache is introduced. Initialization reads four words; each
query reads one bounded record and performs work linear in coefficient count.
Unqueried records remain unvalidated, so loading/inspection does not certify the
whole payload. The existing immutable-asset and concurrent-reader contracts apply.

Format basis: [SPK Required Reading, Type 2 and Type 3](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/req/spk.html#Type%202:%20Chebyshev%20(position%20only)),
revision 2021 AUG 31, consulted 2026-10-13, together with its Segment Order and
Priority sections. The recurrence was implemented independently from the
mathematical Chebyshev series definition, without consulting toolkit source.
Synthetic Swift Testing fixtures cover both types and byte orders, independent
ordinary-power position/velocity references, deliberately independent type 3
velocities, degree-zero and degree-20 series, endpoint/near-boundary selection,
subinterval coverage, reads across physical records, malformed payloads, limits,
file/array precedence and corrupt/unsupported high-priority rejection. These tests
require no assets or CSPICE. Black-box toolkit reference validation remains open;
no cross-toolkit precision or boundary compatibility claim is made yet.

## Relative geometric states and centre chains

```swift
let relative = try kernels.snapshot().spkState(
    targetID: 301, centreID: 399, frameID: 1, atSecondsTDB: 0)
// Moon minus Earth, in kernel frame 1 axes, kilometres and TDB-second rates.
for term in relative.contributions {
    print(term.side, term.kernelIndex, term.state.segmentName)
}
```

The relative-state overload resolves both bodies' precedence-selected centre links
at one finite TDB epoch. Exact NAIF integer identifiers distinguish, for example,
Earth (399) from the Earth–Moon barycentre (3). Negative spacecraft IDs are retained.
ID 0 is the implicit solar-system barycentre and terminates traversal without a
segment. A target queried relative to itself returns zero, without kernel data.

Traversal alternates between target and requested centre using descriptor metadata
until it finds a shared centre. Only the necessary target prefix is added and the
necessary centre prefix subtracted. Thus Moon-relative-to-Earth works with their
barycentre-relative segments alone. No barycentre-to-SSB payload is required or
subtracted. Each contributing position and velocity is accumulated with Neumaier
compensation per axis; non-finite intermediate or final sums fail explicitly.

The algorithm applies later-file/later-array precedence independently at each body.
A higher-priority link may change a body's centre and hence the route through the
graph. Missing data cannot cause a lower-priority covering segment to be selected.
A contributing corrupt or unsupported type propagates its error. Shared ancestors
are not evaluated, so their unsupported/corrupt payloads do not invalidate a local
relative query. Kernels above the common centre are not globally certified by a
successful query. Cycles irrelevant to the required relative path may similarly
remain unused; unresolved cyclic chains fail with the repeated-body path.

`SpiceSPKResolutionError` distinguishes missing body links (no descriptor for the
body), coverage gaps (descriptors exist but none covers the query), cycles, work
limits and frame mismatches. If disconnected branches reveal multiple failures,
cycles take diagnostic priority, followed by coverage gaps and then missing links;
within each category the target branch is reported first. No partial state is
returned. The selected route, rather than whole-file minimum/maximum epochs,
determines whether a query has coverage.

Every contributing segment must have the requested `frameID` (default 1); this
layer does not sum vectors from mismatched axes or infer a rotation. Any nonzero
native frame ID can be explicitly requested, but this does not establish its
physical definition. The SI adapter below resolves transformations using an explicitly supplied library
FrameTree. Text PCK orientations and fixed FK definitions now connect to FrameTree through explicit axes/centre bindings. Automatic graph construction remains deferred. Light-time and stellar aberration corrections remain
outside this geometric API. Quantities remain kilometres and per-TDB-second rates.

The result retains requested IDs/time/axes, common-centre ID, contributing signed
primitive states with source path/array/record/descriptor, the binary-kernel snapshot
index (distinguishing repeated loads), cache scope and metadata work counters.
Neither the result nor evaluation mutates the snapshot. Unloading changes a new
snapshot's selected route; existing snapshots remain independently usable.

`SpiceSPKResolutionOptions` bounds total selected metadata edges explored across
both branches (64 by default) and inspected descriptors (1,000,000 by default).
Exploration can include unused ancestors while finding a common centre; they count
against these budgets, although their coefficient data is never read. Each evaluated
record remains subject to its reader's word limit. Linear metadata searches are
currently uncached; indexing/caching can be introduced separately without changing
precedence or the ownership contract.

Specification basis remains the public [SPK Required Reading](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/req/spk.html),
revision 2021 AUG 31, particularly Introduction (centre graph), State definitions,
Computing States, Data Precedence, Segments and Primitive States, consulted
2026-10-13. The bounded traversal and summation are independent implementations
of relative-vector composition, not translations of toolkit implementation code.

Additional synthetic Swift Testing cases cover multi-link spacecraft/Moon/Earth
states, forward/reverse queries, body-versus-barycentre distinction, identity/SSB
queries, unavailable or unsupported ancestors, precedence changes of centre,
retained snapshots after unload, missing links, coverage gaps, cycles, native-frame
contracts, both work limits, compensated cancellation and overflow. Existing
single-file array-precedence and type 2/3 tests also exercise the shared selector.
Black-box toolkit reference validation remains pending.

## SI ephemeris provider and frame bindings

```swift
let assets = try PhysicalDataProvenance(
    source: "My verified ICRF planetary kernel set", version: "asset-pack-1")
let axes = try SpiceFrameBindings([
    SpiceFrameBinding.icrfJ2000(provenance: assets)
])
let provider = try SpiceSPKStateProvider(
    snapshot: kernels.snapshot(), body: .earth,
    frameBindings: axes, coverage: usableAtomicEpochRange,
    tdbModel: timeModel, dataProvenance: assets, accuracy: .balanced)
let state = try provider.checkedState(at: simulationEpoch)
let detail = try provider.evaluation(at: simulationEpoch)
// detail.chain retains the original selected segments, records, units and frame IDs.
```

`SpiceSPKStateProvider` implements `AccuracyAwareEphemerisProvider`, so it can be
queried directly or placed in an `EphemerisCollection`. Target/origin identifiers
must be `.naif`; negative spacecraft identifiers work without special catalogue
entries. Body centres and barycentres are never aliased. A non-SSB origin requires
an explicit frame centred on that origin; SSB defaults to `CelestialFrames.icrf`.
An explicit output frame may use any connected orientation, including rotating axes.
The caller must ensure that its physical origin agrees with `metadata.origin`;
frame names and NAIF axes bindings alone cannot prove an origin definition.

Frame bindings are mandatory and associate each required NAIF axes ID with an
exact library `Frame` and supplied provenance. Unknown IDs fail rather than acquire
an assumed rotation. Duplicate IDs and ID zero are rejected. SPICE frame 1 is
labelled J2000, but JPL DE4xx and compatible modern products use that label for
ICRF axes. Consequently no general frame-1-to-ICRF assumption is made here.
`icrfJ2000(provenance:)` is an explicit caller assertion that the chosen dataset
uses ICRF axes. Historical mean-equator/equinox data may instead bind ID 1 to
`CelestialFrames.meanEquatorJ2000` and be converted through the frame system.
Bindings apply to the whole retained snapshot; datasets using the same ID with
incompatible physical axes must not be combined under a single binding.

An injected `TDBTimeModel` supplies the TT/atomic instant's TDB offset and rate.
The argument is formed directly from split atomic seconds plus that correction,
without subtracting large Julian dates, before conversion to the scalar Double ET
required by SPK. Metres are `1000 * kilometres`; velocity is
`1000 * velocityKilometresPerTDBSecond * dTDB/dTT`. TT advances at the same elapsed
rate as atomic time in this epoch representation. Large ET values still have the
scalar precision limit inherent in the file format and query argument.

`spkChain` exposes only the necessary primitive links with independent native frame
IDs. The provider converts each link to SI/atomic derivatives and then to output
axes before signed compensated composition. It uses the registered frame tree's
rotation and angular derivatives but discards origin translation and its derivatives,
since each link is already a relative displacement. In particular the transformed
velocity includes `omega × (R * displacement)`. This avoids subtracting an origin
twice and allows chains with mixed native axes. Native bindings and output frames
must be registered in the supplied tree whenever their identities differ.
If every binding uses the output frame itself, no tree is needed.

The frame tree is weakly referenced to avoid cycles when the tree retains an
SPK-backed centre provider. Callers must retain it while transformations are needed.
An unavailable tree produces an explicit error. Frame-transform dependencies can
require origin acceleration even when this adapter uses only orientation terms;
select any existing origin acceleration policy explicitly when building such edges.
No acceleration is fabricated by the SPK provider: output acceleration is nil,
metadata marks it unavailable, and velocity is described as supplied to cover both
type 2 differentiated fits and type 3 independently supplied velocity fits.

Provider metadata records caller-selected atomic coverage, TDB convention, asset
version/provenance and paths, binding provenance, time-model provenance and requested
accuracy. All supplied coefficients are evaluated, so effective accuracy is always
`full`, satisfying lower hints. Construction checks time-model coverage and evaluates
both interval endpoints. This does not certify continuous coverage between endpoints:
interior gaps, unsupported/corrupt higher-priority segments, missing frame mappings
and frame/time-provider failures propagate on the actual query. Select narrower
intervals when a continuous usable route is required. No fallback is automatic.

`evaluation(at:)` returns the SI state, original primitive chain and TDB correction.
The snapshot is retained; unload/reload on the kernel set does not change an existing
provider. Binary assets must remain immutable. No result cache is introduced. Snapshot
cache scope identifies kernel contents only: cache keys for configured providers must
also distinguish target/origin, output frame, frame/time-provider settings and any
mutations of those dependencies (existing cache wrappers supply private scopes).
File reads are synchronized, but concurrent provider evaluation additionally depends
on the supplied time model and frame providers; no blanket thread-safety claim is made.

The axes interpretation follows [NAIF Frames Required Reading, ICRF vs J2000](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/req/frames.html#ICRF%20vs%20J2000),
revision 2021 DEC 31, consulted 2026-10-13. Unit conversion and coordinate-time/rotation derivatives are
independent mathematical implementations using existing library APIs. No NASA/NAIF
implementation source or runtime is used. Synthetic tests verify ordinary-power
position references and numerical atomic-time derivatives under a non-unit TDB clock,
type 3 supplied velocities, large-epoch fractional precision, mixed/rotating native
axes, real geocentric frame construction without duplicate origin subtraction,
EphemerisCollection integration, accuracy/provenance contracts, snapshot retention,
weak tree lifetime, unknown bindings, stale time samples and interior coverage gaps.


## Validated text-kernel definitions

Snapshots now expose `leapSecondData()`, `bodyRadiiMetres(bodyID:)`,
`bodyOrientationParameters(bodyID:)` and `fixedFrameDefinition(frameID:)`.
Decoding is explicit and uses effective pool precedence; semantic errors do not
mutate the pool. Retained snapshots and decoded values survive unload.

LSK decoding validates the complete numeric DELTET parameter set and integral,
strictly ordered post-1971 offset history through the existing `LeapSecondTable`.
`DELTET/DELTA_T_A` must equal the library's fixed 32.184-second TT–TAI relation.
The periodic TDB coefficients are retained rather than automatically replacing an
injected TDB model. Construct the table with a caller-supplied version and exclusive
`validUntilModifiedJulianDay`; an LSK does not certify future leap-second coverage.
Existing UTC conversion handles positive and negative leap seconds and rejects
queries outside that explicit coverage.

PCK radii use conventional kilometres and return positive finite metres.
Orientation decoding retains three polynomial coefficients each for pole RA,
declination and prime meridian, in original degree/day/Julian-century units, plus
available body-specific periodic amplitudes, long-axis offset and epoch/frame
parameters. Planetary-system phase-angle arrays remain available through raw pool
inspection. This accessor decodes parameters; the orientation provider described below performs evaluation with system-level reference metadata and phase angles. Provenance is retained on orientation
and LSK definitions; radii provenance is also included by the sourced physical-value adapter below.
GM conversion requires the explicit conventional-unit/convention assertion described below.

FK decoding supports explicit nonzero IDs, class 4, matching class ID, numeric
centre, named relative axes, and `MATRIX`, `ANGLES` or `QUATERNION` specification addressed by numeric ID.
The nine values are column-major and map local vectors to relative axes by left
multiplication. Matrices must be orthonormal and have determinant +1 (tolerance
1e-10). Name-addressed TK parameters, built-in lookup,
frame aliases, dependency graph construction and other frame classes remain
deferred; unsupported classes/representations produce explicit errors.
The relative name remains unresolved until frame-provider construction.

Independent implementation references, consulted 2026-10-13:

- [Kernel Required Reading](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/req/kernel.html),
  revision 2021 DEC 29, Variable Value Rules (unscaled date literals).
- [Time Required Reading](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/req/time.html),
  revision 2021 DEC 23, Relationship between TDB and UTC (LSK parameter meanings).
- [PCK Required Reading](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/req/pck.html),
  revision 2021 DEC 24, Text PCK Contents and Text PCK Kernel Variable Names.
- [Frames Required Reading](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/req/frames.html),
  revision 2021 DEC 31, frame specification, TK Frames and matrix definitions.

Synthetic Swift Testing fixtures verify independent date offsets, leap-second
labels and expiry, malformed histories, radius units, retained coefficients,
column-major rotation direction, handedness and unsupported definitions. No
NASA/NAIF implementation source or runtime dependency is used.


## Text-kernel providers

`bodyReferenceRadii(bodyID:axesDescription:provenance:)` returns the existing
`BodyReferenceRadii` type. `bodyGravitationalParameter(bodyID:convention:provenance:)`
returns `GravitationalParameter`, converting conventional kernel km³/s² to m³/s².
The caller confirms those input units and chooses the coordinate-time convention;
no mass is inferred. Both retain contributing paths plus caller asset version and
reference. Radii axes descriptions remain explicit, including any ellipsoid offset.

`SpicePCKOrientationProvider` conforms to `OrientationProvider` and works with
`FrameTree.makeBodyFixedFrame`. Supply the body ID, its system ID (planetary
barycentre for a planet/satellite, body ID otherwise), explicit `SpiceFrameBindings`,
atomic epoch coverage, TDB model and asset provenance. The system ID is required,
so unusual NAIF relationships are never guessed.

The provider evaluates the full quadratic pole/prime-meridian polynomials and
optional sine/cosine amplitudes. System phase polynomials support degrees 0–3,
with a default degree of 1; malformed arrays or amplitudes exceeding available
phases fail explicitly. System-level `CONSTANTS_REF_FRAME`/`CONSTS_REF_FRAME` and
`CONSTANTS_JED_EPOCH`/`CONSTS_JED_EPOCH` are honored, defaulting to frame 1 and
JD 2451545.0. Conflicting long/short aliases fail instead of silently choosing one.
The inertial axes ID must have an explicit binding. Body-specific overrides for
planets/satellites do not replace system metadata, as specified in PCK Appendix B.

Orientations map reference vectors into body axes. Rates and accelerations are
analytic in elapsed atomic seconds, including both TDB derivative factors and the
composition of changing pole and spin axes. Coverage is caller-selected, bounded
by the TDB model, and evaluated at construction endpoints. It is not a certified
physical accuracy range. Retained coefficients remain valid after kernel unload;
`cacheScope` identifies their snapshot. No implicit fallback or global registration
occurs, and injected time-model thread safety remains the caller's responsibility.
`LONG_AXIS` describes shape axes and does not alter prime-meridian body axes.

FK Euler definitions support all seven documented angular units. Their matrices
use the documented coordinate-rotation signs and angle order. Quaternion definitions
require a unit, scalar-first Hamilton quaternion. All representations normalize to
validated column-major local-to-relative matrices; the orientation adapter inverts
this mapping for the library's reference-to-child convention.

`SpiceFixedOrientationProvider` binds a fixed definition to exact relative axes.
`definition.makeFrame(in:relativeAxes:centre:provenance:accelerationPolicy:)` creates
a child with those axes and the declared object's physical centre. The caller
asserts that the supplied axes represent `relativeFrameName` and that the centre
provider represents `centreID` in those exact axes. Origin derivatives use the
existing explicit acceleration policy; missing acceleration is not silently zeroed.
Callers construct dependencies in order and can add resulting identities to SPK
frame bindings. No owning-tree references are captured by these providers.

Additional public specification sections used: PCK Models for the Sun, Planets,
and Satellites, phase polynomials and Appendix B; Frames Euler/quaternion TK
definitions; and [Rotations Required Reading](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/req/rotation.html),
revision 2017 MAR 09, consulted 2026-10-13, Quaternion arithmetic and
quaternion-to-matrix definitions. These are public
mathematical/behavior descriptions, not toolkit implementation source.

Synthetic tests cover SI constants/provenance, coupled angular derivatives,
non-unit TDB clock rates, periodic phase degree 3, independent basis rotations,
alternate system epochs/axes, malformed dependencies, all FK angle units, Euler
order, quaternion layout, matrix inversion including half turns, moving-centre
translation and snapshot lifetime after unload. Additional binary PCK types and CSPICE
oracle references remain subsequent work.


## Binary PCK type 2

The kernel set accepts modern IEEE `DAF/PCK` assets alongside SPKs, directly or
through manifests. DAF PCK files require ND=2/NI=5; typed descriptors expose the
frame **class** ID, inertial base-axes ID, data type and inclusive TDB coverage.
Class IDs are selected explicitly; they are not inferred from a body or frame ID.
Unknown positive segment types remain inspectable, but selected types other than
2 fail with `unsupportedType`.

`SpiceBinaryPCKSegment` validates the four-word INIT/INTLEN/RSIZE/N directory,
record shape, descriptor containment and finite/resolvable time grid. Queries read
only the selected coefficient record and validate MID/RADIUS and finite values.
Evaluation uses the full constant coefficient and analytically differentiated
Chebyshev recurrences through second order. Record reads use DAF resource limits;
no whole-asset copy or result cache is introduced. Internal boundaries select the
later record; the inclusive final endpoint uses the last record. Storage-rounding
tolerance applies to record metadata, not extrapolation outside descriptors.

`file.pckAngles(array:atSecondsTDB:)` returns native angles in radians, their first
and second derivatives per TDB second, and path/array/name/record provenance.
`snapshot.pckAngles(frameClassID:atSecondsTDB:maximumDescriptorChecks:)` selects
covering PCK files and arrays in reverse load/list order. Missing class IDs and
coverage gaps have separate errors. A corrupt or unsupported selected segment
never silently substitutes older binary or text data. Mixed SPK/PCK pools preserve
independent segment selection; PCK arrays do not consume SPK descriptor budgets.

### Binary angle convention and documentation limits

The PCK Required Reading type-2 diagram labels the series RA/DEC/W, but those
labels alone do not specify the exact Euler-coordinate mapping or angle units.
The public type-2 [PCKEUL Detailed Output contract](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/FORTRAN/spicelib/pckeul.html#Detailed_Output)
(version 2.1.0, 2021 AUG 20, consulted 2026-10-13) defines the mapping as
`[angle3]3 [angle2]1 [angle1]3`. The public
[Rotations Required Reading](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/req/rotation.html)
(revision 2017 MAR 09) specifies radians and coordinate-rotation signs for this
notation. The reader treats the stored values as direct Euler angles; it does not
add the text model's RA/node or declination/co-latitude offsets. The type-2
implementation uses these public behavior contracts and independently implemented
mathematics, without consulting toolkit implementation source.

PCK Required Reading revision 2021 DEC 24, consulted 2026-10-13, supplies the
binary descriptor/directory/record layout and latest-file/latest-segment rules.
Type 3's variable-interval layout is not fully specified there; it remains deferred
rather than borrowing an implementation. Type 20 has a different rate/integration
format and remains a separate extension. No decoded comment delimiter format or
unprovided derivative data is inferred.

### Orientation provider

`SpiceBinaryPCKOrientationProvider` implements `OrientationProvider` over a retained
snapshot. Supply a frame class ID, exact output reference frame, explicit dataset
axes bindings, caller-selected atomic coverage, TDB model and asset provenance.
Construction validates both coverage endpoints; interior gaps still fail on query.
`evaluation(at:)` retains native selected-record provenance and the TDB correction.
The provider computes the 3-1-3 rotation and composed angular velocity/acceleration,
including TDB rate and acceleration. It uses every coefficient supplied by the asset.

Selected base axes may differ across segments. An externally retained FrameTree
converts the output reference axes into those selected axes before the Euler
rotation, including frame rotation rates. Origin translations are omitted because
this provider supplies orientation only. Output/base axes must be explicitly
pseudo-inertial; their physical interpretation follows caller bindings. The tree
reference is weak to avoid owning-tree cycles. Standard `makeBodyFixedFrame` then
adds the caller's centre state using the existing explicit acceleration policy.

Asset paths and caller versions, axes interpretation and TDB model provenance are
retained. Use the snapshot cache scope for external invalidation. Old providers
remain valid after unload while assets remain immutable. File reads are synchronized;
frame/time provider concurrency follows the injected components. Binary-to-text
orientation fallback is not selected automatically by this adapter.

Synthetic Swift Testing fixtures cover both byte orders and odd NI padding,
independent polynomial values and first/second derivatives, degree 0/60,
cross-storage-record reads, boundaries, invalid descriptors/directories/records,
resource bounds, missing/gapped class IDs, file/array precedence, mixed SPK/PCK
loads, unload/snapshot lifetime, Euler matrix signs, coupled angular derivatives,
non-unit TDB rates and explicit base-frame conversion without origin translation.
Large local assets and CSPICE remain optional; black-box oracle fixtures are not
yet included.


## Type-1 spacecraft clocks

`SpiceSCLKClock` decodes effective pool variables for an explicit negative NAIF
spacecraft ID. Kernel variable suffixes use the negated ID (for example, spacecraft
-77 selects `_77`). `KPL/SCLK` loads directly or through a manifest, with the same
atomic-load, precedence, unload and snapshot semantics as other text kernels.
Common SCLK metadata dates such as `@04-SEP-1990//4:23:00` are supported; these
pool dates remain unscaled and do not determine clock coverage or asset version.

The supported format uses 1–16 fields, positive integral moduli, nonnegative
integral offsets and all five documented output delimiters. Readings are explicit
`SpiceSCLKReading(partition:fields:)` values; partitions are one-based. Formatting
and parsing require an explicit partition and all canonical integral fields.
Omitted fields/partitions, empty fields and overflowing field values are rejected;
this is a deliberately stricter subset than CSPICE's permissive string parser.
Typed and string readings both validate against the selected partition and
caller coverage. No repeated clock value is assigned to a partition by guessing.

Partition starts/ends are raw clock ticks. Encoded ticks accumulate `end-start`
from preceding partitions, without adding an extra tick: the end of one partition
and start of the next represent the same encoded instant. Both representations
are available when a caller explicitly selects their partition. Civil decoding
requires integral ticks; continuous correlation APIs preserve fractional ticks.
Counts and cumulative partition lengths are bounded to the exact integer domain
of Double (up to 2^53); this does not imply arbitrary fractional precision at large
counts. Zero-duration partitions are outside the supported subset.

Clock coefficients are triplets of encoded ticks, parallel seconds and rate.
The rate is seconds per most-significant-field count and is divided by the product
of lower-field moduli to obtain seconds per tick. Tick anchors must be strictly
increasing and rates strictly positive. Later anchors apply at exact breakpoints.
The final rate applies through the caller's explicit `encodedTickCoverage`, which
must lie within declared partition history and begin at or after the first anchor.
No coverage before the first anchor or outside that selected interval is inferred.
Parallel time magnitudes are bounded below 1e15 seconds as a representation policy.

`sample(atEncodedTicks:)` exposes the parallel label, local rate and coefficient
index. Piecewise correlations may be discontinuous. `sample(atParallelSeconds:)`
checks all covered branches and requires a unique solution: gaps and duplicate
parallel times are explicit errors. This differs from the permissive NAIF inverse
mapping for non-invertible correlations. Forward selection uses binary search;
inverse work is linear in the bounded pool's coefficient count.

Parallel system 1 is TDB (also the documented default when absent); system 2 is
TDT/TT. `epoch(atEncodedTicks:using:)` and `sample(at:using:)` bridge these labels
to atomic `Epoch`. TDB conversion requires an explicit `TDBTimeModel` and validates
its coverage/sample contracts. TT conversion uses the existing atomic TT relation.
No UTC conversion, LSK lookup or guessed time model is performed. Other clock types
and parallel systems fail explicitly. Angular data and CK interpolation are not
part of this clock adapter.

```swift
let clock = try SpiceSCLKClock(
    snapshot: kernels.snapshot(), spacecraftID: -77,
    encodedTickCoverage: 0...200, provenance: clockAssetProvenance)
let reading = try clock.reading(from: "1/15:1")
let ticks = try clock.encodedTicks(for: reading)
let instant = try clock.epoch(atEncodedTicks: ticks, using: tdbModel)
```

The decoded model is immutable and Sendable, retaining effective source paths,
caller asset version/reference and snapshot cache scope. It remains valid after
unload without rereading files; time-model concurrency follows the injected model.

Public implementation references, consulted 2026-10-13:

- [SCLK Required Reading](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/req/sclk.html),
  revision 2021 SEP 04: Encoded SCLK, Partitions, Type 1 SCLK conversion algorithms,
  parallel time system, format constants and coefficient assignments.
- [SCENCD public Detailed Output contract](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/scencd_c.html#Detailed_Output):
  adjacent partition endpoints have identical encoded ticks. This public behavior
  description supplies the boundary rule; no toolkit implementation source is used.

Synthetic tests independently verify mixed-radix fields and offsets, reset
partitions/shared endpoints, explicit string formats, rate normalization,
fractional ticks and breakpoint selection, discontinuity gaps and ambiguous
inverse mappings, TDB/TT Epoch conversion, stale time-model contracts, malformed
constants, representational bounds, manifest loading, precedence/unload and retained
clock provenance. CK formats and attitude/frame integration remain next.


## Discrete spacecraft attitude: CK type 1

CK files now participate in kernel-set loading, manifests, snapshots and unloads.
Their typed descriptors contain encoded SCLK bounds, instrument ID, reference
frame ID, segment type and angular-velocity flag. Modern little- and big-endian
IEEE DAF files are supported. CK types other than 1, 2 and 3 are inspectable but evaluation
throws an explicit unsupported-type error.

`SpiceCKSegment` implements the public type-1 layout: scalar-first Hamilton
quaternions, optional three-component angular velocity, strictly increasing
encoded times, midpoint directory entries and the final record count. The reader
validates layout, times and directory, then reads only the selected pointing
record. Record count is bounded (100,000 by default), with DAF read limits still
applied. Non-unit/non-finite selected quaternions and invalid rates fail explicitly.

`SpiceKernelSnapshot.ckPointing` searches latest loaded files and latest segments
first. A discrete gap permits an older segment; malformed or unsupported selected
segments fail. Exact matching is the default. Optional nearest-sample tolerance is
in encoded ticks, applies within descriptor bounds, and equal-distance ties select
the later record (an explicit library policy). The result preserves both requested
and selected times. There is no interpolation or extrapolation, and queries outside
a descriptor are not expanded by tolerance. This bounded subset is not a complete
replacement for all toolkit tolerance-search behavior.

The returned quaternion maps base-reference coordinates into instrument coordinates.
Optional physical angular velocity is in base-reference axes, radians per second.
`SpiceCKAttitudeProvider` requires an explicit `SpiceSCLKClock`, instrument/body
association, native axes bindings and rate time-scale convention (`atomic` or `tdb`).
The caller is responsible for associating that clock with the instrument; instrument
IDs do not imply spacecraft IDs. TDB clock/rate conversion requires a TDB model.
The provider accepts exact discrete samples only, converts rates into body axes and
atomic seconds, and produces the inverse body-to-reference attitude quaternion.
Alternate inertial base axes use the frame tree; retain that tree externally.

Missing angular velocity stays nil in primitive pointing and throws from the
`AttitudeProvider` adapter. Type-1 samples cannot establish angular acceleration;
for type 1 the adapter returns nil acceleration, so the existing attitude-to-orientation/frame
adapter rejects them instead of inventing a zero derivative. Nearest-sample rendering
should use the primitive result and retain its selected timestamp. Type-2 continuous
rotation is supported below; type-3 interpolation intervals are supported below.

Swift Testing fixtures independently construct both byte orders and verify direction,
nearest-sample timestamps, directories, malformed records, unsupported types,
clock conversion, snapshot retention and derivative availability. These are synthetic
public-format checks, not a claim of independent real-kernel/CSPICE validation.

Format and behavior reference: [NAIF CK Required Reading](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/req/ck.html),
particularly C-matrix/velocity conventions and type-1 segment structure. No NAIF
implementation source or foreign runtime dependency was used.


## Continuous spacecraft attitude: CK type 2

`SpiceCKType2Segment` reads eight-word pointing records, interval start/stop
arrays and the midpoint directory. There is no count footer: record count follows
from the published segment-length formula. Reads and record counts remain bounded.
Intervals must have positive duration and cannot overlap except at endpoints;
the interval starting at a shared endpoint wins. Within a covered interval,
reference-axis physical spin advances the quaternion analytically. Zero spin is
supported. Gaps remain uncovered, even with a nonzero primitive query tolerance;
this subset does not snap type-2 gaps to endpoints. Snapshot priority and selected
malformed-data errors follow the type-1 policy.

Primitive pointing preserves the stored radians-per-model-second angular velocity
and exposes `secondsPerEncodedTick`. The attitude adapter differentiates the actual
quaternion trajectory: angular velocity is multiplied by the record seconds/tick,
divided by the selected SCLK correlation's parallel seconds/tick, and multiplied by
dTDB/dTT when the clock's parallel system is TDB. Angular acceleration includes the
TDB second derivative. The explicit `angularRateTimeScale` argument continues to
apply to discrete type-1 rates; type 2 uses its record factor and clock correlation.
This matters when nominal CK factors differ from measured clock correlations.

These derivatives are local to the selected CK interval and SCLK correlation branch.
At shared boundaries they describe the selected branch, not a claim that the attitude
or its derivatives are continuous across a rate change, clock jump or interval boundary.
Type-2 acceleration can therefore feed `AttitudeOrientationProvider` and the existing
body-fixed frame machinery; type-1 acceleration remains unavailable. Frame-axis
rotation derivatives are composed using the existing transform system.

Independent synthetic Swift Testing fixtures cover both byte orders, spin direction,
shared endpoints, gaps, directory/resource limits, malformed intervals and rate factors,
and finite-difference verification of orientation and angular acceleration under TT and
nonlinear TDB clocks. Real-kernel oracle checks remain roadmap work.
The implementation uses only the [public NAIF CK type-2 specification](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/req/ck.html#Data%20Type%202),
without NAIF implementation source or external C/Fortran libraries. Type 3 is supported below.


## Interpolated spacecraft attitude: CK type 3

`SpiceCKType3Segment` supports the public record layout, both time directories,
interval starts and the two count footers. Times are strictly increasing and every
interval start must match a pointing sample; the first interval starts at the first
sample. Each interval ends at the sample before the next interval starts, including
single-sample intervals. Layout, membership, directory and resource limits are checked
before selected quaternion/rate records are read. Other CK types remain unsupported.

Pointing uses shortest-arc fixed-axis rotation between adjacent samples in an interval.
Quaternion sign changes do not change the path; exact half-turn ambiguity uses a
canonical quaternion sign. Exact samples preserve their stored orientation. No
interpolation crosses a gap. With nonzero tolerance, the closest gap endpoint can be
returned, preserving its selected encoded time; ties use the later endpoint. As with
type 1, tolerance search stays within descriptor bounds, rather than implementing the
entire toolkit's expanded descriptor search behavior.

Optional stored angular velocities are linearly interpolated in reference axes,
independently of pointing, as required by the public specification. They may disagree
with the derivative of the interpolated quaternion. The primitive keeps them optional;
`SpiceCKAttitudeProvider` requires available rates and the explicit atomic/TDB rate
convention. It returns nil angular acceleration, so the existing orientation/frame
adapter rejects that incomplete derivative set. This avoids presenting independent
stored rates as a derivative-consistent rotating frame. A future opt-in provider could
derive rates from the pointing interpolation, with a separate contract and explicit
handling of interpolation knots and singleton intervals.

Synthetic Swift Testing fixtures cover both byte orders, independent rate interpolation,
shortest arcs, quaternion sign changes, tiny rotations, noncommuting initial attitudes,
interval gaps and singleton intervals, directory validation and malformed membership,
and the attitude adapter's derivative-availability contract. These do not replace the
pending real-kernel/oracle validation.

Reference: [NAIF CK type-3 format and interpolation algorithm](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/req/ck.html#Data%20Type%203).
Only public format and behavior documentation was used; no toolkit implementation
source or external native library was introduced.


## Ephemeris fallback integration

`FallbackEphemeris` now accepts an optional SPK ephemeris alongside optional DE441,
model and Kepler providers. The caller selects `.de441ThenSPK` or `.spkThenDE441`
within the file-backed tier. Missing assets are represented by omitted providers;
no files are loaded automatically. Exact body, origin and frame identities are required.

The default SPK policy in this overload permits only missing chain links and segment
coverage gaps to fall through. Strict `.declaredCoverageOnly` propagates these too.
Corrupt/unsupported data, cycles, limits, frame and time-model failures, and state
contract errors always propagate. DE441/model/Kepler query failures retain their existing
strict behavior. Results report the selected source and skipped sources/reasons; if no
backup succeeds, the last SPK availability error is preserved.

See [EPHEMERIDES.md](EPHEMERIDES.md) for the API and provenance/coverage contracts.
Swift Testing covers explicit file priority, absent sources, gap diagnostics, strict
error handling and real synthetic SPK interior gaps, including retained snapshots and
a selected unsupported segment. No optional large asset or CSPICE installation is needed.

## Pinned black-box validation

Small checked-in synthetic kernels and outputs from CSPICE N0067 now independently
check SPK types 2/3 and chained centres, CK types 1/2/3 (including missing rates), binary
PCK type-2 orientation, leap labels and fractional SCLK correlations. The normal Swift
Testing suite needs neither CSPICE nor network access. Existing independent public-format
fixtures retain coverage of byte orders, malformed inputs, precedence, gaps and unload/reload.

See [fixture provenance and tolerances](Tests/CelMekTests/Fixtures/spice-oracle/README.md)
and [optional generator](Scripts/generate_spice_oracle.py). Toolkit/wrapper versions,
asset hashes and a reference-output hash are pinned. The optional generator uses only
public writer/query APIs with a binary oracle; no toolkit implementation source was used.
These tests validate supported synthetic formats and conventions, not deferred types or
all real mission assets. No CSPICE runtime dependency was added to the library.
