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

/// A tolerance hint: models may use a higher accuracy than requested.
/// Shared caller intent, not a universal physical error bound or runtime guarantee.
/// Each algorithm must document its settings, supported intervals and measured error.
public enum AccuracyLevel: String, CaseIterable, Hashable, Sendable {
  /// Prioritize computational cost, allowing the model's largest supported approximation.
  case low
  /// A model-specific compromise between computational cost and approximation error.
  case balanced
  /// Prefer smaller approximation error while permitting documented simplifications.
  case high
  /// Complete evaluation of the implemented model, including any intrinsic truncation.
  /// Does not imply exact physics, full original published theory, or the newest data.
  case full
}

/// Opt-in contract reporting requested hints and effective evaluation modes. Algorithm-specific
/// controls (such as derivative intervals and TDB models) remain explicit and independent.
public protocol AccuracyAwareEphemerisProvider: EphemerisProvider {
  /// Requested tolerance hint.
  var accuracy: AccuracyLevel { get }
  /// Actual model mode; may be more accurate than the hint.
  var effectiveAccuracy: AccuracyLevel { get }
  /// Implemented evaluation modes, not a restriction on accepted hints.
  var availableAccuracyLevels: Set<AccuracyLevel> { get }
}
