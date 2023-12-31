//
// SPDX-License-Identifier: Apache-2.0
//
// Copyright 2020 Mattias Holm
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

import XCTest
@testable import CelMek;

final class goffin2000_tests : XCTestCase {
  func testGoffin2000() {
    
    let res = goffin2000(2448908.5);
    
    XCTAssertEqual(res.x, 232.74009.deg, accuracy:  0.0000005);
    XCTAssertEqual(res.y, 14.58769.deg, accuracy:  0.0000005);
    XCTAssertEqual(res.z, 29.711383, accuracy: 0.000001);
  }
}


