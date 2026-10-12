//
// SPDX-License-Identifier: Apache-2.0
//
// Copyright 2026 Mattias Holm
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
// http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
//
extension BodyCatalogue {
  /// Selected values from NAIF gm_de440.tpc, 2022-12-14 snapshot. The DE440 values
  /// apply to codes 1...10, 199, 299, 301, 399. Other planet-centre values come from
  /// the satellite solutions identified in that file, not from DE440 system GMs.
  /// All are converted from km³/s² to m³/s². No masses or radii are inferred.
  public static let jplGM2022: BodyCatalogue = {
    let reference = "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/pck/gm_de440.tpc"
    let identity = try! PhysicalDataProvenance(source: "NAIF body codes",version: "named Solar System subset",
      reference: "https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/req/naif_ids.html")
    let de440 = try! PhysicalDataProvenance(source: "JPL DE440 ASTRO-VALUES via NAIF gm_de440.tpc",
      version: "2022-12-14",reference: reference)
    let satellites = try! PhysicalDataProvenance(source: "JPL satellite solutions via NAIF gm_de440.tpc",
      version: "2022-12-14; source 2 in kernel",reference: reference)
    let values: [(BodyIdentifier,String,CelestialBodyKind,Double)] = [
      (.mercuryBarycentre,"Mercury system barycentre",.barycentre,2.2031868551400003e4),
      (.venusBarycentre,"Venus system barycentre",.barycentre,3.2485859200000000e5),
      (.earthMoonBarycentre,"Earth–Moon barycentre",.barycentre,4.0350323562548019e5),
      (.marsBarycentre,"Mars system barycentre",.barycentre,4.2828375815756102e4),
      (.jupiterBarycentre,"Jupiter system barycentre",.barycentre,1.2671276409999998e8),
      (.saturnBarycentre,"Saturn system barycentre",.barycentre,3.7940584841799997e7),
      (.uranusBarycentre,"Uranus system barycentre",.barycentre,5.7945563999999985e6),
      (.neptuneBarycentre,"Neptune system barycentre",.barycentre,6.8365271005803989e6),
      (.plutoBarycentre,"Pluto system barycentre",.barycentre,9.7550000000000000e2),
      (.sun,"Sun",.star,1.3271244004127942e11),
      (.mercury,"Mercury",.planet,2.2031868551400003e4),
      (.venus,"Venus",.planet,3.2485859200000000e5),
      (.moon,"Moon",.naturalSatellite,4.9028001184575496e3),
      (.earth,"Earth",.planet,3.9860043550702266e5),
      (.mars,"Mars",.planet,4.282837362069909e4),
      (.jupiter,"Jupiter",.planet,1.266865319003704e8),
      (.saturn,"Saturn",.planet,3.793120623436167e7),
      (.uranus,"Uranus",.planet,5.793951256527211e6),
      (.neptune,"Neptune",.planet,6.835103145462294e6),
      (.pluto,"Pluto",.dwarfPlanet,8.696138177608748e2)
    ]
    var bodies = [try! CelestialBody(identifier: .solarSystemBarycentre,name: "Solar System barycentre",
      kind: .barycentre,provenance: identity)]
    for (id,name,kind,gm) in values {
      let planetaryCentre = [BodyIdentifier.mars,.jupiter,.saturn,.uranus,.neptune,.pluto].contains(id)
      let parameter = try! GravitationalParameter(cubicMetresPerSecondSquared: gm*1e9,
        convention: .tdbCompatible,provenance: planetaryCentre ? satellites : de440)
      bodies.append(try! CelestialBody(identifier: id,name: name,kind: kind,
        properties: BodyPhysicalProperties(gravitationalParameter: parameter),provenance: identity))
    }
    return try! BodyCatalogue(bodies: bodies)
  }()
}
