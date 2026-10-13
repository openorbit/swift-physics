// swift-tools-version:5.9
// The swift-tools-version declares the minimum version of Swift required to build this package.

import PackageDescription

let package = Package(
    name: "swift-physics",
    products: [
        // Products define the executables and libraries produced by a package, and make them visible to other packages.
        .library(
            name: "Physics",
            targets: ["CelMek", "Aerodynamics", "Math"]),
        .library(name: "CelMek", targets: ["CelMek"]),
        .library(name: "Math", targets: ["Math"]),
        .library(name: "Aerodynamics", targets: ["Aerodynamics"]),
    ],
    dependencies: [
        // Dependencies declare other packages that this package depends on.
        // .package(url: /* package url */, from: "1.0.0"),
      .package(url: "https://github.com/jpsim/Yams.git", from: "5.0.6"),
      .package(url: "https://github.com/apple/swift-docc-plugin", from: "1.1.0"),
    ],
    targets: [
        // Targets are the basic building blocks of a package. A target can define a module or a test suite.
        // Targets can depend on other targets in this package, and on products in packages which this package depends on.
        .target(
            name: "CelMek",
            dependencies: ["Math"]),
        .target(
            name: "Math",
            dependencies: []),
        .target(
            name: "Aerodynamics",
            dependencies: ["Yams"]),
        .executableTarget(
            name: "TwoBodyGameDemo",
            dependencies: ["CelMek"],
            path: "Examples/TwoBodyGameDemo"),
        .executableTarget(
            name: "AttitudeDynamicsDemo",
            dependencies: ["CelMek", "Math"],
            path: "Examples/AttitudeDynamicsDemo"),
        .executableTarget(
            name: "Voyager2GameDemo",
            dependencies: ["CelMek", "Math"],
            path: "Examples/Voyager2GameDemo",
            exclude: ["README.md"]),
        .executableTarget(
            name: "VSOP87Characterization",
            dependencies: ["CelMek", "Math"],
            path: "Tools/VSOP87Characterization"),
        .testTarget(
            name: "CelMekTests",
            dependencies: ["CelMek"],
            resources: [.copy("Fixtures")]),
        .testTarget(
            name: "AerodynamicsTests",
            dependencies: ["Aerodynamics"]),
        .testTarget(
            name: "MathTests",
            dependencies: ["Math"]),

    ]
)
