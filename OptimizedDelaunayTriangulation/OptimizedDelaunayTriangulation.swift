//
//  OptimizedDelaunayTriangulation.swift
//  OptimizedDelaunayTriangulation
//
//  Created by Paul Ollivier on 27/10/2024.
//

import Foundation
import simd

// MARK: - Data Structures

public struct DPoint: Hashable {
    public let x: Double
    public let y: Double
    let index: Int
    
    public init(x: Double, y: Double, index: Int = 0) {
        self.x = x
        self.y = y
        self.index = index
    }
}

public struct DTriangle: Hashable {
    public let point1: DPoint
    public let point2: DPoint
    public let point3: DPoint
    
    public init(point1: DPoint, point2: DPoint, point3: DPoint) {
        self.point1 = point1
        self.point2 = point2
        self.point3 = point3
    }
}

internal struct DCircumcircle: Hashable {
    let point1: DPoint
    let point2: DPoint
    let point3: DPoint
    let x: Double
    let y: Double
    let rsqr: Double
}

internal struct DEdge: Hashable {
    let point1: DPoint
    let point2: DPoint
    
    init(point1: DPoint, point2: DPoint) {
        if point1.index <= point2.index {
            self.point1 = point1
            self.point2 = point2
        } else {
            self.point1 = point2
            self.point2 = point1
        }
    }
}

// MARK: - Memory Management

final class TriangulationContext {
    var points: ContiguousArray<DPoint>
    var open: Set<DCircumcircle>
    var completed: Set<DCircumcircle>
    var edges: [DEdge: Int]
    let pointCount: Int
    
    init(capacity: Int) {
        self.pointCount = capacity
        self.points = ContiguousArray()
        self.points.reserveCapacity(capacity + 3)
        self.open = Set(minimumCapacity: capacity)
        self.completed = Set(minimumCapacity: capacity * 2)
        self.edges = Dictionary(minimumCapacity: capacity * 3)
    }
}

// MARK: - Geometric Utilities

@inline(__always)
private func fastCircumcircle(_ i: DPoint, _ j: DPoint, _ k: DPoint) -> DCircumcircle {
    #if arch(arm64) || arch(x86_64)
    let p1 = SIMD2<Double>(i.x, i.y)
    let p2 = SIMD2<Double>(j.x, j.y)
    let p3 = SIMD2<Double>(k.x, k.y)
    
    let d = 2 * (p1.x * (p2.y - p3.y) + p2.x * (p3.y - p1.y) + p3.x * (p1.y - p2.y))
    
    if abs(d) < Double.ulpOfOne {
        // Collinear case
        let center = SIMD2<Double>((p1.x + p3.x) * 0.5, (p1.y + p3.y) * 0.5)
        let diff = p1 - center
        let rsqr = diff.x * diff.x + diff.y * diff.y
        return DCircumcircle(point1: i, point2: j, point3: k,
                           x: center.x, y: center.y, rsqr: rsqr)
    }
    
    let ux = ((p1.x * p1.x + p1.y * p1.y) * (p2.y - p3.y) +
              (p2.x * p2.x + p2.y * p2.y) * (p3.y - p1.y) +
              (p3.x * p3.x + p3.y * p3.y) * (p1.y - p2.y)) / d
    
    let uy = ((p1.x * p1.x + p1.y * p1.y) * (p3.x - p2.x) +
              (p2.x * p2.x + p2.y * p2.y) * (p1.x - p3.x) +
              (p3.x * p3.x + p3.y * p3.y) * (p2.x - p1.x)) / d
    
    let center = SIMD2<Double>(ux, uy)
    let diff = p1 - center
    let rsqr = diff.x * diff.x + diff.y * diff.y
    
    return DCircumcircle(point1: i, point2: j, point3: k,
                        x: center.x, y: center.y, rsqr: rsqr)
    #else
    return circumcircleScalar(i, j, k)
    #endif
}

private func circumcircleScalar(_ i: DPoint, _ j: DPoint, _ k: DPoint) -> DCircumcircle {
    let x1 = i.x, y1 = i.y
    let x2 = j.x, y2 = j.y
    let x3 = k.x, y3 = k.y
    
    let fabsy1y2 = abs(y1 - y2)
    let fabsy2y3 = abs(y2 - y3)
    
    var xc: Double = 0.0
    var yc: Double = 0.0
    
    if fabsy1y2 < Double.ulpOfOne {
        let m2 = -((x3 - x2) / (y3 - y2))
        let mx2 = (x2 + x3) * 0.5
        let my2 = (y2 + y3) * 0.5
        xc = (x2 + x1) * 0.5
        yc = m2 * (xc - mx2) + my2
    } else if fabsy2y3 < Double.ulpOfOne {
        let m1 = -((x2 - x1) / (y2 - y1))
        let mx1 = (x1 + x2) * 0.5
        let my1 = (y1 + y2) * 0.5
        xc = (x3 + x2) * 0.5
        yc = m1 * (xc - mx1) + my1
    } else {
        let m1 = -((x2 - x1) / (y2 - y1))
        let m2 = -((x3 - x2) / (y3 - y2))
        let mx1 = (x1 + x2) * 0.5
        let my1 = (y1 + y2) * 0.5
        let mx2 = (x2 + x3) * 0.5
        let my2 = (y2 + y3) * 0.5
        
        if abs(m1 - m2) < Double.ulpOfOne {
            xc = (x1 + x3) * 0.5
            yc = (y1 + y3) * 0.5
        } else {
            xc = (m1 * mx1 - m2 * mx2 + my2 - my1) / (m1 - m2)
            yc = m1 * (xc - mx1) + my1
        }
    }
    
    let dx = x2 - xc
    let dy = y2 - yc
    let rsqr = dx * dx + dy * dy
    
    return DCircumcircle(point1: i, point2: j, point3: k, x: xc, y: yc, rsqr: rsqr)
}

@inline(__always)
private func optimizedSupertriangle(_ points: ContiguousArray<DPoint>) -> [DPoint] {
    var xmin = Double.greatestFiniteMagnitude
    var ymin = Double.greatestFiniteMagnitude
    var xmax = -Double.greatestFiniteMagnitude
    var ymax = -Double.greatestFiniteMagnitude
    
    points.withUnsafeBufferPointer { buffer in
        for point in buffer {
            xmin = min(xmin, point.x)
            ymin = min(ymin, point.y)
            xmax = max(xmax, point.x)
            ymax = max(ymax, point.y)
        }
    }
    
    let dx = xmax - xmin
    let dy = ymax - ymin
    let dmax = max(dx, dy)
    let xmid = xmin + dx * 0.5
    let ymid = ymin + dy * 0.5
    
    let margin = dmax * 20
    return [
        DPoint(x: xmid - margin, y: ymid - dmax, index: -1),
        DPoint(x: xmid, y: ymid + margin, index: -2),
        DPoint(x: xmid + margin, y: ymid - dmax, index: -3)
    ]
}

// MARK: - Sorting

private func swapValues(_ indices: UnsafeMutablePointer<Int>, _ i: Int, _ j: Int) {
    let temp = indices[i]
    indices[i] = indices[j]
    indices[j] = temp
}

private func quicksortPoints(_ points: ContiguousArray<DPoint>, indices: UnsafeMutablePointer<Int>, low: Int, high: Int) {
    guard low < high else { return }
    
    let pivot = partitionPoints(points, indices: indices, low: low, high: high)
    quicksortPoints(points, indices: indices, low: low, high: pivot - 1)
    quicksortPoints(points, indices: indices, low: pivot + 1, high: high)
}

private func partitionPoints(_ points: ContiguousArray<DPoint>, indices: UnsafeMutablePointer<Int>, low: Int, high: Int) -> Int {
    let pivot = points[indices[high]].x
    var i = low - 1
    
    for j in low..<high {
        if points[indices[j]].x <= pivot {
            i += 1
            swapValues(indices, i, j)
        }
    }
    
    swapValues(indices, i + 1, high)
    return i + 1
}

// MARK: - Point Processing

private func processPoint(context: TriangulationContext, currentPoint: DPoint) {
    var trianglesToRemove = Set<DCircumcircle>()
    context.edges.removeAll(keepingCapacity: true)
    
    for circle in context.open {
        let dx = currentPoint.x - circle.x
        
        if dx > 0 && dx * dx > circle.rsqr {
            context.completed.insert(circle)
            trianglesToRemove.insert(circle)
            continue
        }
        
        let dy = currentPoint.y - circle.y
        if dx * dx + dy * dy - circle.rsqr > Double.ulpOfOne {
            continue
        }
        
        trianglesToRemove.insert(circle)
        
        let edges = [
            DEdge(point1: circle.point1, point2: circle.point2),
            DEdge(point1: circle.point2, point2: circle.point3),
            DEdge(point1: circle.point3, point2: circle.point1)
        ]
        
        for edge in edges {
            context.edges[edge, default: 0] += 1
        }
    }
    
    context.open.subtract(trianglesToRemove)
    
    for (edge, count) in context.edges where count == 1 {
        let newCircle = fastCircumcircle(edge.point1, edge.point2, currentPoint)
        context.open.insert(newCircle)
    }
}

// MARK: - Main Triangulation Function

public func triangulate(_ points: [DPoint]) -> [DTriangle] {
    guard points.count >= 3 else { return [] }
    
    let context = TriangulationContext(capacity: points.count)
    
    // Remove duplicates while preserving order
    var seen = Set<DPoint>()
    context.points = ContiguousArray(points.enumerated().compactMap { index, point in
        let pointWithIndex = DPoint(x: point.x, y: point.y, index: index)
        return seen.insert(pointWithIndex).inserted ? pointWithIndex : nil
    })
    
    let n = context.points.count
    guard n >= 3 else { return [] }
    
    // Sort points
    let indices = UnsafeMutablePointer<Int>.allocate(capacity: n)
    defer { indices.deallocate() }
    
    for i in 0..<n {
        indices[i] = i
    }
    
    quicksortPoints(context.points, indices: indices, low: 0, high: n - 1)
    
    // Add supertriangle points
    let superPoints = optimizedSupertriangle(context.points)
    context.points.append(contentsOf: superPoints)
    
    // Initialize with supertriangle
    let superTriangle = fastCircumcircle(context.points[n],
                                         context.points[n + 1],
                                         context.points[n + 2])
    context.open.insert(superTriangle)
    
    // Process points
    let threshold = 1000
    if n > threshold {
        // Create a regular array for concurrent access
        let pointIndices = (0..<n).map { indices[$0] }
        
        DispatchQueue.concurrentPerform(iterations: n) { i in
            let currentPoint = context.points[pointIndices[i]]
            processPoint(context: context, currentPoint: currentPoint)
        }
    } else {
        for i in 0..<n {
            let currentPoint = context.points[indices[i]]
            processPoint(context: context, currentPoint: currentPoint)
        }
    }
    
    // Final processing
    context.completed.formUnion(context.open)
    
    // Remove triangles connected to supertriangle
    return context.completed.compactMap { circle in
        if circle.point1.index < 0 || circle.point2.index < 0 || circle.point3.index < 0 {
            return nil
        }
        return DTriangle(point1: circle.point1,
                         point2: circle.point2,
                         point3: circle.point3)
    }
}
