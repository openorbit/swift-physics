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
public struct EphemerisRequest {
  public let body: BodyIdentifier
  public let epoch: Epoch
  public init(body: BodyIdentifier, epoch: Epoch) { self.body = body; self.epoch = epoch }
}

extension EphemerisProvider {
  /// Ordered exact queries. Duplicates are preserved; no sorting or time quantization.
  /// Fail-fast: the first error propagates and no partial array is returned. Earlier
  /// successful evaluations may populate an explicitly selected cache.
  public func states(at epochs: [Epoch]) throws -> [CartesianState] {
    try epochs.map { try checkedState(at: $0) }
  }
}

extension EphemerisCollection {
  public func states(of body: BodyIdentifier, at epochs: [Epoch]) throws -> [CartesianState] {
    try provider(for: body).states(at: epochs)
  }
  /// Results match request order even across different origins, frames or coverage.
  /// Inspect each body's metadata before combining coordinates.
  public func states(for requests: [EphemerisRequest]) throws -> [CartesianState] {
    try requests.map { try state(of: $0.body,at: $0.epoch) }
  }
}
