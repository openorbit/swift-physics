// SPDX-License-Identifier: Apache-2.0
import Foundation

/// Kernel-set provenance uses resolved paths, including expanded Windows 8.3
/// aliases and symlinks. Compare the complete resolved path, not just a basename.
func resolvedSpiceTestPath(_ url: URL) -> String {
    url.standardizedFileURL.resolvingSymlinksInPath().path
}
