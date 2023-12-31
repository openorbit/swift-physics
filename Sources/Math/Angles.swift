//
// SPDX-License-Identifier: Apache-2.0
//
// Copyright 2020-2022 Mattias Holm
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
import Foundation

let MINUTES_OF_ARC_PER_DEGREE : Double = 60.0
let SECONDS_OF_ARC_PER_DEGREE : Double = 3600.0
let DEGREES_PER_HOUR : Double = 15.0
let SECONDS_PER_HOUR : Double = 3600.0
let MINUTES_PER_HOUR : Double = 60.0
let SECONDS_PER_MINUTE : Double = 60.0

public extension Double {
  var mdeg: Double {
    return self / 1000.0 * (.pi / 180.0)
  }
  var udeg: Double {
    return self.µdeg
  }
  var µdeg: Double {
    return self / 1000_000.0 * (.pi / 180.0)
  }
  var deg: Double {
    return self * (.pi / 180.0)
  }
  var rad: Double {
    return self
  }
  var arcsec: Double {
    return (self / SECONDS_OF_ARC_PER_DEGREE).deg
  }

  var asDeg: Double {
    return self * (180.0 / .pi)
  }
  var degAsHours: Double {
    return self / 15.0
  }
  var hoursAsDeg: Double {
    return self * 15.0
  }
  var asArcSec: Double {
    return self.asDeg * SECONDS_OF_ARC_PER_DEGREE
  }

}

public func normalize(secondsInDay: Double) -> Double {
  let tmp = secondsInDay.remainder(dividingBy: 3600*24)
  if tmp < 0 {
    return tmp + 3600*24
  }
  return tmp
}

public func normalize(radians: Double) -> Double {
  let tmp = radians.remainder(dividingBy: 2.0 * .pi)
  if tmp < 0.0 {
    return tmp + 2.0 * .pi
  }
  return tmp
}

public func normalize(degrees: Double) -> Double {
  let tmp = degrees.remainder(dividingBy: 360.0)
  if tmp < 0.0 {
    return tmp + 360.0
  }
  return tmp
}

public func normalize(seconds: Double) -> Double {
  let tmp = seconds.remainder(dividingBy: 60.0)
  if tmp < 0.0 {
    return tmp + 60.0
  }
  return tmp
}

public func normalize(minutes: Int) -> Int {
  let tmp = minutes % 60
  if tmp < 0 {
    return tmp + 60
  }
  return tmp
}

public func normalize(hours: Int) -> Int {
  let tmp = hours % 24
  if tmp < 0 {
    return tmp + 24
  }
  return tmp
}

@frozen
public struct HourAngle {
  public var hours: Int
  public var minutes: Int
  public var seconds: Double
  public init(hours: Int, minutes: Int, seconds: Double) {
    self.hours = hours
    self.minutes = minutes
    self.seconds = seconds
  }
}


public extension HourAngle {
  init(seconds: Double) {
    let tmpSeconds = normalize(secondsInDay: seconds)
    self.seconds = tmpSeconds
    self.minutes = Int(tmpSeconds / 60)
    self.hours = self.minutes / 60
    self.minutes = normalize(minutes: self.minutes)
    self.seconds = normalize(seconds: self.seconds)
  }

  init(degrees: Double) {
    let h = degrees / DEGREES_PER_HOUR
    let m = modf(h).1 * 60.0
    let s = modf(m).1 * 60.0

    hours = Int(h)
    minutes = Int(m)
    seconds = s
  }

  static func +(lhs: HourAngle, rhs: HourAngle) -> HourAngle {
    var ha = HourAngle(hours: lhs.hours + rhs.hours, minutes: lhs.minutes + rhs.minutes, seconds: lhs.seconds + rhs.seconds)

    ha.minutes += Int(ha.seconds / 60)
    ha.seconds = normalize(seconds: ha.seconds)
    ha.hours += Int(ha.minutes / 60)
    ha.minutes = normalize(minutes: ha.minutes)
    ha.hours = normalize(hours: ha.hours)
    return ha
  }

  func toDeg() -> Double {
    let totalHours = Double(abs(self.hours)) + Double(self.minutes) / MINUTES_PER_HOUR + self.seconds / SECONDS_PER_HOUR

    if self.hours >= 0 {
      return totalHours.hoursAsDeg
    } else {
      return -totalHours.hoursAsDeg
    }
  }

  func toRad() -> Double {
    return toDeg().deg
  }
}

@frozen
public struct AngleOfArc {
  public var degrees: Int
  public var minutes: Int
  public var seconds: Double

  public init(degrees: Int, minutes: Int, seconds: Double) {
    self.degrees = degrees
    self.minutes = minutes
    self.seconds = seconds
  }
}

public extension AngleOfArc {
  func toDeg() -> Double {
    let totalDegrees = Double(abs(self.degrees)) + Double(self.minutes) / MINUTES_OF_ARC_PER_DEGREE + self.seconds / SECONDS_OF_ARC_PER_DEGREE
    if self.degrees >= 0 {
      return totalDegrees
    } else {
      return -totalDegrees
    }
  }

  func toRad() -> Double {
    return toDeg().deg
  }
}
