// SPDX-License-Identifier: Apache-2.0
import Foundation
import CelMek
import Math

// swift run --disable-sandbox AttitudeDynamicsDemo
// Newtonian Earth-orbit example with a continuously firing, off-centre engine.
let reference = Frame(name: "Earth inertial",pseudoInertial: true)
let body = Frame(parent: reference,name: "Spacecraft axes",pseudoInertial: false)
let tree = try FrameTree(root: reference)
let gm = try GravitationalParameter(cubicMetresPerSecondSquared: 3.986004418e14,
    convention: .newtonian,provenance: PhysicalDataProvenance(source: "Earth game model",version: "1"))
let context = try RotationalEvaluationContext(inertia: RigidBodyInertia(xx: 200,yy: 300,zz: 400),bodyFrame: body,
    referenceContext: ForceEvaluationContext(frame: reference,origin: .earth))
let spacecraft = try SpacecraftState(cartesian: CartesianState(epoch: .j2000,frame: reference,
    position: SIMD3(7e6,0,0),velocity: SIMD3(0,sqrt(gm.cubicMetresPerSecondSquared/7e6),0)),
    mass: SpacecraftMassProperties(dryMassKilograms: 1000,propellantMassKilograms: 100))
let attitude = try AttitudeState(epoch: .j2000,referenceFrame: reference,bodyFrame: body,
    bodyToReference: Quaternion(angle: .pi/2,axis: SIMD3(0,0,1)))
let engine = try BodyMountedThruster(bodyFrame: body,thrustNewtons: 1,directionBody: SIMD3(1,0,0),
    positionBodyMetres: SIMD3(0,1,0),effectiveExhaustVelocityMetresPerSecond: 3000)
let dynamics = try CoupledSpacecraftDynamics(context: context,
    forceModels: [CentralGravity(body: .earth,gravitationalParameter: gm)],wrenchModels: [engine])
let trajectory = try CoupledSpacecraftPropagator(initialState: CoupledSpacecraftState(spacecraft: spacecraft,attitude: attitude),
    dynamics: dynamics,attitudeOptions: .preset(for: .full),numericalOptions: .preset(for: .full))
try tree.registerSpacecraftFrame(spacecraft: trajectory.spacecraftProvider,attitude: trajectory.attitudeProvider)
struct Sample: Encodable {
    let seconds: Double
    let positionMetres: [Double]
    let forceReferenceNewtons: [Double]
    let angularVelocityBodyRadiansPerSecond: [Double]
    let quaternionBodyToReferenceXYZW: [Double]
    let totalMassKilograms: Double
}
func array(_ vector: SIMD3<Double>) -> [Double] { [vector.x,vector.y,vector.z] }
let samples = try [0.0,5,10,20].map { seconds in
    let epoch = try Epoch.j2000.advanced(by: CelMek.Duration(seconds: seconds))
    let state = try trajectory.state(at: epoch), force = try engine.contribution(for: state,in: context)
    let q = state.attitude.bodyToReference
    return Sample(seconds: seconds,positionMetres: array(state.spacecraft.cartesian.position),
        forceReferenceNewtons: array(force.forceReferenceNewtons),
        angularVelocityBodyRadiansPerSecond: array(state.attitude.angularVelocityBodyRadiansPerSecond),
        quaternionBodyToReferenceXYZW: [q.x,q.y,q.z,q.w],totalMassKilograms: state.spacecraft.mass.totalMassKilograms)
}
let encoder = JSONEncoder()
encoder.outputFormatting = [.prettyPrinted,.sortedKeys]
print(String(decoding: try encoder.encode(samples),as: UTF8.self))
