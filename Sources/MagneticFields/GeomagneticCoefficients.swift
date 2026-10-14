// SPDX-License-Identifier: Apache-2.0
import Foundation

public enum GeomagneticCoefficientError: Error, Equatable {
    case malformedLine(Int), invalidHeader, duplicateCoefficient, incompleteCoefficients, unsupportedDegree
}

/// Immutable Schmidt semi-normalized coefficients in nT and nT/year.
/// Parsers accept standard WMM .COF and IGRF g/h epoch-table text, independently of file paths.
struct GeomagneticCoefficients: Sendable {
    let epochs: [Double]
    let values: [[Double]] // Alternating g,h at triangular index.
    let secularVariation: [Double]
    let endYear: Double
    let degree: Int
    let name: String
    static func index(_ n: Int, _ m: Int) -> Int { n*(n+1)/2+m }
    static func validDegree(_ n: Int, _ m: Int) throws {
        guard n >= 1, n <= 64 else { throw GeomagneticCoefficientError.unsupportedDegree }
        guard m >= 0, m <= n else { throw GeomagneticCoefficientError.incompleteCoefficients }
    }
    func at(_ year: Double) -> [Double] {
        if year >= epochs.last! {
            return zip(values.last!,secularVariation).map { $0 + (year-epochs.last!)*$1 }
        }
        let lower = epochs.lastIndex(where: { $0 <= year })!
        let weight = (year-epochs[lower])/(epochs[lower+1]-epochs[lower])
        return zip(values[lower],values[lower+1]).map { $0 + weight*($1-$0) }
    }
    static func lines(_ text: String) -> [(Int,[Substring])] {
        text.split(whereSeparator: \.isNewline).enumerated().compactMap { i,line in
            let tokens = line.split(whereSeparator: \.isWhitespace)
            guard let first = tokens.first, !first.hasPrefix("#") else { return nil }
            return (i+1,tokens)
        }
    }
    static func wmm(_ text: String, validityYears: Double) throws -> Self {
        let rows = lines(text)
        guard let header = rows.first?.1, header.count >= 2, let epoch = Double(header[0]), epoch.isFinite,
              validityYears.isFinite, validityYears > 0, (epoch+validityYears).isFinite else {
            throw GeomagneticCoefficientError.invalidHeader
        }
        var records: [Int:(Double,Double,Double,Double)] = [:], degree = 0
        var ended = false
        for (line,t) in rows.dropFirst() {
            if t[0].hasPrefix("9999") { ended = true; continue }
            guard !ended, t.count == 6, let n = Int(t[0]), let m = Int(t[1]) else {
                throw GeomagneticCoefficientError.malformedLine(line)
            }
            try validDegree(n,m)
            let numbers = t.dropFirst(2).compactMap { Double($0) }
            guard numbers.count == 4, numbers.allSatisfy(\.isFinite), m != 0 || (numbers[1] == 0 && numbers[3] == 0) else {
                throw GeomagneticCoefficientError.malformedLine(line)
            }
            let key = index(n,m)
            guard records[key] == nil else { throw GeomagneticCoefficientError.duplicateCoefficient }
            records[key] = (numbers[0],numbers[1],numbers[2],numbers[3]); degree = max(degree,n)
        }
        guard ended, degree > 0, records.count == index(degree,degree) else {
            throw GeomagneticCoefficientError.incompleteCoefficients
        }
        var values = [Double](repeating: 0,count: 2*(index(degree,degree)+1)), rate = values
        for (i,v) in records { values[2*i] = v.0; values[2*i+1] = v.1; rate[2*i] = v.2; rate[2*i+1] = v.3 }
        return Self(epochs: [epoch],values: [values],secularVariation: rate,endYear: epoch+validityYears,
                    degree: degree,name: String(header[1]))
    }
    static func igrf(_ text: String) throws -> Self {
        let rows = lines(text)
        guard let headerIndex = rows.firstIndex(where: { $0.1.first == "g/h" }) else {
            throw GeomagneticCoefficientError.invalidHeader
        }
        let header = rows[headerIndex].1
        guard header.count >= 5 else { throw GeomagneticCoefficientError.invalidHeader }
        let epochs = header.dropFirst(3).dropLast().compactMap { Double($0) }
        let interval = header.last!.split(separator: "-")
        guard epochs.count == header.count-4, epochs.count >= 2, epochs.allSatisfy(\.isFinite),
              zip(epochs,epochs.dropFirst()).allSatisfy({ $0 < $1 }), interval.count == 2,
              let start = Double(interval[0]), start == epochs.last!, let suffix = Double(interval[1]) else {
            throw GeomagneticCoefficientError.invalidHeader
        }
        var end = interval[1].count == 2 ? floor(start/100)*100+suffix : suffix
        if interval[1].count == 2, end <= start { end += 100 }
        guard end.isFinite, end > start else { throw GeomagneticCoefficientError.invalidHeader }
        var records: [Int:[Double]] = [:], degree = 0
        for (line,t) in rows.dropFirst(headerIndex+1) {
            guard t.count == epochs.count+4, t[0] == "g" || t[0] == "h",
                  let n = Int(t[1]), let m = Int(t[2]) else { throw GeomagneticCoefficientError.malformedLine(line) }
            try validDegree(n,m)
            guard t[0] != "h" || m > 0 else { throw GeomagneticCoefficientError.malformedLine(line) }
            let key = 2*index(n,m)+(t[0] == "h" ? 1 : 0)
            guard records[key] == nil else { throw GeomagneticCoefficientError.duplicateCoefficient }
            let numbers = t.dropFirst(3).compactMap { Double($0) }
            guard numbers.count == epochs.count+1, numbers.allSatisfy(\.isFinite) else {
                throw GeomagneticCoefficientError.malformedLine(line)
            }
            records[key] = numbers; degree = max(degree,n)
        }
        guard degree > 0, records.count == degree*(degree+2) else { throw GeomagneticCoefficientError.incompleteCoefficients }
        let count = 2*(index(degree,degree)+1)
        var values = epochs.map { _ in [Double](repeating: 0,count: count) }, rate = [Double](repeating: 0,count: count)
        for (key,numbers) in records {
            for i in epochs.indices { values[i][key] = numbers[i] }
            rate[key] = numbers.last!
        }
        return Self(epochs: epochs,values: values,secularVariation: rate,endYear: end,degree: degree,name: "IGRF epoch table")
    }
}
