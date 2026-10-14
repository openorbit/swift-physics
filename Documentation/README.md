# Documentation

Public documentation is authored in five DocC catalogs, alongside the module sources:

| Module | Catalog |
| --- | --- |
| CelMek | [CelMek.docc](../Sources/CelMek/CelMek.docc/CelMek.md) |
| Math | [Math.docc](../Sources/Math/Math.docc/Math.md) |
| Aerodynamics | [Aerodynamics.docc](../Sources/Aerodynamics/Aerodynamics.docc/Aerodynamics.md) |
| Atmospheres | [Atmospheres.docc](../Sources/Atmospheres/Atmospheres.docc/Atmospheres.md) |
| MagneticFields | [MagneticFields.docc](../Sources/MagneticFields/MagneticFields.docc/MagneticFields.md) |

The module landing pages curate conceptual guides and API symbols. Internal article
links use DocC references; repository code links point to the GitHub source tree.
Roadmaps, working audits, fixture notes and developer architecture policy remain
repository documents. The former Antora source pages/configuration have been replaced;
the original bibliography and ERFA license notice remain in this directory.

## Build one module

Use the existing Swift-DocC package plugin:

```sh
mkdir -p .build/documentation
swift package generate-documentation --target CelMek --disable-indexing \
    --output-path .build/documentation/CelMek.doccarchive
```

In Xcode, select the package scheme and use Product → Build Documentation.

## Build the package archive

With a Swift toolchain that includes `docc merge` (Xcode 16 or newer):

```sh
python3 Tools/Documentation/build.py
```

This generates one archive for each library module and merges them into
`.build/documentation/Physics.doccarchive`, with a package landing page. Build logs
are retained next to the archives. DocC warnings are treated as build failures.
The script selects the navigator-index option supported by the installed DocC toolchain. Open the combined archive in Xcode to browse it.
Output stays under `.build`; no generated archives are checked into the repository.

## Browse in a web browser

Serve the archive over HTTP:

```sh
python3 Tools/Documentation/preview.py
```

Open [the local documentation](http://127.0.0.1:8080/documentation). The preview binds
only to localhost; stop it with Ctrl-C. Use `--port 8081` if the default port is busy,
or `--archive /path/to/Physics.doccarchive` for another combined archive.

Opening the archive's `index.html` directly with `file://` produces a blank page:
the DocC web viewer needs HTTP to load its JavaScript and JSON assets from the
archive root. Open the archive itself in Xcode, or use this browser preview.

When running inside an existing restricted sandbox, pass `--disable-sandbox` to
avoid nesting SwiftPM's plugin sandbox. Use `--output /path/to/directory` to select
another output directory. The script grants the plugin write access only to that
selected output directory.

## Editing guides

Add articles to the owning module's `.docc` catalog and curate them under a `## Topics`
section on its module landing page. Use `<doc:ArticleName>` for conceptual guides and
double-backtick symbol links for API references. Build the affected module and check
its DocC diagnostics before submitting changes.
