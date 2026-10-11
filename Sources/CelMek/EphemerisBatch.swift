// SPDX-License-Identifier: Apache-2.0

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
