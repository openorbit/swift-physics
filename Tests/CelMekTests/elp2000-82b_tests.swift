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

import Testing
@testable import CelMek;

@Suite
struct elp2000_82b_tests {
  @Test
  func testELP2000() {
    
    let res = elp2000_82b(2448724.5);
    
    #expect(abs((res.x) - (133.162659.deg)) <= 0.00000002)
    #expect(abs((res.y) - (-3.229127.deg)) <= 0.000000005)
    #expect(abs((res.z) - (368409.7e3)) <= 100.0)
  }
}




