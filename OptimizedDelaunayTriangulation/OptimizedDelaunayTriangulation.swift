import Foundation
import simd

// MARK: - Data Structures

public struct DPoint: Hashable {
    public let x: Double
    public let y: Double
    public let index: Int
    
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

struct DCircumcircle: Hashable {
    let point1: DPoint
    let point2: DPoint
    let point3: DPoint
    let x: Double
    let y: Double
    let rsqr: Double
    
    // Pre-compute hash for better performance in Sets
    private let precomputedHash: Int
    
    init(point1: DPoint, point2: DPoint, point3: DPoint, x: Double, y: Double, rsqr: Double) {
        self.point1 = point1
        self.point2 = point2
        self.point3 = point3
        self.x = x
        self.y = y
        self.rsqr = rsqr
        
        // Pre-compute hash combining all values
        var hasher = Hasher()
        hasher.combine(x)
        hasher.combine(y)
        hasher.combine(point1.index)
        hasher.combine(point2.index)
        hasher.combine(point3.index)
        self.precomputedHash = hasher.finalize()
    }
    
    func hash(into hasher: inout Hasher) {
        hasher.combine(precomputedHash)
    }
}

struct DEdge: Hashable {
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
    private(set) var points: UnsafeMutableBufferPointer<DPoint>
    var open: Set<DCircumcircle>
    var completed: Set<DCircumcircle>
    var edges: [DEdge: Int]
    let pointCount: Int
    
    init(capacity: Int) {
        self.pointCount = capacity
        let allocatedPoints = UnsafeMutablePointer<DPoint>.allocate(capacity: capacity + 3)
        self.points = UnsafeMutableBufferPointer(start: allocatedPoints, count: capacity + 3)
        self.open = Set(minimumCapacity: capacity)
        self.completed = Set(minimumCapacity: capacity * 2)
        self.edges = Dictionary(minimumCapacity: capacity * 3)
    }
    
    deinit {
        points.deallocate()
    }
}

// MARK: - SIMD Optimized Geometric Calculations

private func fastCircumcircle(_ i: DPoint, _ j: DPoint, _ k: DPoint) -> DCircumcircle {
    let p1 = SIMD2<Double>(i.x, i.y)
    let p2 = SIMD2<Double>(j.x, j.y)
    let p3 = SIMD2<Double>(k.x, k.y)
    
    let d = 2 * simd_cross(p2 - p1, p3 - p1).z
    
    if abs(d) < Double.ulpOfOne {
        let center = (p1 + p3) * 0.5
        let diff = p1 - center
        let rsqr = simd_dot(diff, diff)
        return DCircumcircle(point1: i, point2: j, point3: k,
                           x: center.x, y: center.y, rsqr: rsqr)
    }
    
    let sq1 = simd_dot(p1, p1)
    let sq2 = simd_dot(p2, p2)
    let sq3 = simd_dot(p3, p3)
    
    let center = SIMD2<Double>(
        (sq1 * (p2.y - p3.y) + sq2 * (p3.y - p1.y) + sq3 * (p1.y - p2.y)) / d,
        (sq1 * (p3.x - p2.x) + sq2 * (p1.x - p3.x) + sq3 * (p2.x - p1.x)) / d
    )
    
    let diff = p1 - center
    let rsqr = simd_dot(diff, diff)
    
    return DCircumcircle(point1: i, point2: j, point3: k,
                        x: center.x, y: center.y, rsqr: rsqr)
}

// MARK: - Optimized Sorting with Block-Based Quicksort

private func blockQuicksort(_ points: UnsafeMutableBufferPointer<DPoint>, low: Int, high: Int) {
    let blockSize = 64
    var stack = [(low: Int, high: Int)](repeating: (0,0), count: 32)
    var stackTop = 0
    
    stack[stackTop] = (low, high)
    stackTop += 1
    
    while stackTop > 0 {
        stackTop -= 1
        let (left, right) = stack[stackTop]
        
        if right - left <= blockSize {
            // Use insertion sort for small arrays
            for i in (left + 1)...right {
                let temp = points[i]
                var j = i - 1
                while j >= left && points[j].x > temp.x {
                    points[j + 1] = points[j]
                    j -= 1
                }
                points[j + 1] = temp
            }
            continue
        }
        
        // Median-of-three pivot selection
        let mid = left + (right - left) / 2
        if points[left].x > points[mid].x { points.swapAt(left, mid) }
        if points[left].x > points[right].x { points.swapAt(left, right) }
        if points[mid].x > points[right].x { points.swapAt(mid, right) }
        
        let pivot = points[mid].x
        var i = left
        var j = right
        
        while i <= j {
            while points[i].x < pivot { i += 1 }
            while points[j].x > pivot { j -= 1 }
            if i <= j {
                points.swapAt(i, j)
                i += 1
                j -= 1
            }
        }
        
        if j - left > right - i {
            if left < j {
                stack[stackTop] = (left, j)
                stackTop += 1
            }
            if i < right {
                stack[stackTop] = (i, right)
                stackTop += 1
            }
        } else {
            if i < right {
                stack[stackTop] = (i, right)
                stackTop += 1
            }
            if left < j {
                stack[stackTop] = (left, j)
                stackTop += 1
            }
        }
    }
}

// MARK: - Parallel Point Processing

private func processPointBatch(context: TriangulationContext, startIndex: Int, endIndex: Int) {
    var localEdges = [DEdge: Int](minimumCapacity: (endIndex - startIndex) * 3)
    var trianglesToRemove = Set<DCircumcircle>()
    
    for i in startIndex..<endIndex {
        let currentPoint = context.points[i]
        trianglesToRemove.removeAll(keepingCapacity: true)
        localEdges.removeAll(keepingCapacity: true)
        
        // Process circles and collect edges
        for circle in context.open {
            let dx = currentPoint.x - circle.x
            
            if dx > 0 && dx * dx > circle.rsqr {
                context.completed.insert(circle)
                trianglesToRemove.insert(circle)
                continue
            }
            
            let dy = currentPoint.y - circle.y
            if dx * dx + dy * dy - circle.rsqr <= Double.ulpOfOne {
                trianglesToRemove.insert(circle)
                
                localEdges[DEdge(point1: circle.point1, point2: circle.point2), default: 0] += 1
                localEdges[DEdge(point1: circle.point2, point2: circle.point3), default: 0] += 1
                localEdges[DEdge(point1: circle.point3, point2: circle.point1), default: 0] += 1
            }
        }
        
        context.open.subtract(trianglesToRemove)
        
        // Create new triangles
        for (edge, count) in localEdges where count == 1 {
            let newCircle = fastCircumcircle(edge.point1, edge.point2, currentPoint)
            context.open.insert(newCircle)
        }
    }
}

// MARK: - Main Triangulation Function

public func triangulate(_ points: [DPoint]) -> [DTriangle] {
    guard points.count >= 3 else { return [] }
    
    let context = TriangulationContext(capacity: points.count)
    
    // Remove duplicates while preserving order using a more efficient approach
    var seen = Set<DPoint>(minimumCapacity: points.count)
    var uniqueCount = 0
    
    for (index, point) in points.enumerated() {
        let pointWithIndex = DPoint(x: point.x, y: point.y, index: index)
        if seen.insert(pointWithIndex).inserted {
            context.points[uniqueCount] = pointWithIndex
            uniqueCount += 1
        }
    }
    
    guard uniqueCount >= 3 else { return [] }
    
    // Sort points using block quicksort
    blockQuicksort(UnsafeMutableBufferPointer(start: context.points.baseAddress!, count: uniqueCount), low: 0, high: uniqueCount - 1)
    
    // Add supertriangle points efficiently
    let (minX, maxX, minY, maxY) = points.reduce((Double.infinity, -Double.infinity, Double.infinity, -Double.infinity)) { result, point in
        (min(result.0, point.x), max(result.1, point.x), min(result.2, point.y), max(result.3, point.y))
    }
    
    let dx = maxX - minX
    let dy = maxY - minY
    let dmax = max(dx, dy)
    let xmid = minX + dx * 0.5
    let ymid = minY + dy * 0.5
    let margin = dmax * 20
    
    let superPoints = [
        DPoint(x: xmid - margin, y: ymid - dmax, index: -1),
        DPoint(x: xmid, y: ymid + margin, index: -2),
        DPoint(x: xmid + margin, y: ymid - dmax, index: -3)
    ]
    
    for (i, point) in superPoints.enumerated() {
        context.points[uniqueCount + i] = point
    }
    
    // Initialize with supertriangle
    context.open.insert(fastCircumcircle(context.points[uniqueCount],
                                       context.points[uniqueCount + 1],
                                       context.points[uniqueCount + 2]))
    
    // Process points in parallel batches
    let batchSize = max(1000, uniqueCount / ProcessInfo.processInfo.activeProcessorCount)
    DispatchQueue.concurrentPerform(iterations: (uniqueCount + batchSize - 1) / batchSize) { batch in
        let start = batch * batchSize
        let end = min(start + batchSize, uniqueCount)
        processPointBatch(context: context, startIndex: start, endIndex: end)
    }
    
    // Final processing
    context.completed.formUnion(context.open)
    
    // Efficiently filter and transform results
    return context.completed.compactMap { circle in
        guard circle.point1.index >= 0,
              circle.point2.index >= 0,
              circle.point3.index >= 0 else {
            return nil
        }
        return DTriangle(point1: circle.point1,
                        point2: circle.point2,
                        point3: circle.point3)
    }
}
