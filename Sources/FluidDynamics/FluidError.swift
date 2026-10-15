// SPDX-License-Identifier: Apache-2.0
public enum FluidError: Error, Equatable {
    case invalidParameter, invalidGeometry, invalidFlow, invalidCoefficients
    case malformedData, duplicateCondition, incompleteGrid, outsideCoverage, nonFiniteResult
    case missingFlowProperty, unsupportedQuery, nonConvergence, workBudgetExceeded
}

