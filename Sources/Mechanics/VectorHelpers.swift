// SPDX-License-Identifier: Apache-2.0
// Copyright 2026 Mattias Holm
import Foundation
func mechanicsFinite(_ p:SIMD3<Double>)->Bool {[p.x,p.y,p.z].allSatisfy(\.isFinite)}
func mechanicsCross(_ a:SIMD3<Double>,_ b:SIMD3<Double>)->SIMD3<Double> {SIMD3(a.y*b.z-a.z*b.y,a.z*b.x-a.x*b.z,a.x*b.y-a.y*b.x)}
