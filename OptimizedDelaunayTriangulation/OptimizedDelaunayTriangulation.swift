import Foundation
import simd

// MARK: - Public Types

public struct DPoint: Hashable {
    public let x: Double
    public let y: Double
    @usableFromInline let index: Int
    
    @usableFromInline
    init(x: Double, y: Double, index: Int) {
        self.x = x
        self.y = y
        self.index = index
    }
    
    public init(x: Double, y: Double) {
        self.init(x: x, y: y, index: 0)
    }
}

public struct DTriangle: Hashable {
    public let point1: DPoint
    public let point2: DPoint
    public let point3: DPoint
    
    @usableFromInline
    init(point1: DPoint, point2: DPoint, point3: DPoint) {
        self.point1 = point1
        self.point2 = point2
        self.point3 = point3
    }
}

// MARK: - Internal Types

@usableFromInline
struct DCircumcircle: Hashable {
    @usableFromInline let point1: DPoint
    @usableFromInline let point2: DPoint
    @usableFromInline let point3: DPoint
    @usableFromInline let x: Double
    @usableFromInline let y: Double
    @usableFromInline let rsqr: Double
    
    @usableFromInline
    init(point1: DPoint, point2: DPoint, point3: DPoint, x: Double, y: Double, rsqr: Double) {
        self.point1 = point1
        self.point2 = point2
        self.point3 = point3
        self.x = x
        self.y = y
        self.rsqr = rsqr
    }
}

@usableFromInline
struct DEdge: Hashable {
    @usableFromInline let point1: DPoint
    @usableFromInline let point2: DPoint
    
    @usableFromInline
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

@usableFromInline
final class TriangulationContext {
    let points: UnsafeMutableBufferPointer<DPoint>
    var open: Set<DCircumcircle>
    var completed: Set<DCircumcircle>
    private var edgePool: [DEdge: Int]
    private var removePool: Set<DCircumcircle>
    
    @usableFromInline
    init(capacity: Int) {
        let allocatedPoints = UnsafeMutablePointer<DPoint>.allocate(capacity: capacity + 3)
        self.points = UnsafeMutableBufferPointer(start: allocatedPoints, count: capacity + 3)
        // Use power of 2 for better hash table performance
        let openCapacity = max(16, (capacity / 4).nextPowerOf2())
        let completedCapacity = max(16, (capacity / 2).nextPowerOf2())
        self.open = Set(minimumCapacity: openCapacity)
        self.completed = Set(minimumCapacity: completedCapacity)
        self.edgePool = Dictionary(minimumCapacity: max(32, capacity))
        self.removePool = Set(minimumCapacity: max(32, capacity / 4))
    }
    
    deinit {
        points.deallocate()
    }
    
    @inline(__always)
    func getReusableCollections() -> (edges: [DEdge: Int], toRemove: Set<DCircumcircle>) {
        edgePool.removeAll(keepingCapacity: true)
        removePool.removeAll(keepingCapacity: true)
        return (edgePool, removePool)
    }
    
    @inline(__always)
    func updateEdges(_ edges: [DEdge: Int]) {
        edgePool = edges
    }
    
    @inline(__always)
    func updateRemoveSet(_ toRemove: Set<DCircumcircle>) {
        removePool = toRemove
    }
}

// MARK: - Geometric Calculations

@usableFromInline
@inline(__always)
func fastCircumcircle(_ i: DPoint, _ j: DPoint, _ k: DPoint) -> DCircumcircle {
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
    
    let invD = 1.0 / d
    let sq1 = simd_dot(p1, p1)
    let sq2 = simd_dot(p2, p2)
    let sq3 = simd_dot(p3, p3)
    
    let center = SIMD2<Double>(
        (sq1 * (p2.y - p3.y) + sq2 * (p3.y - p1.y) + sq3 * (p1.y - p2.y)) * invD,
        (sq1 * (p3.x - p2.x) + sq2 * (p1.x - p3.x) + sq3 * (p2.x - p1.x)) * invD
    )
    
    let diff = p1 - center
    let rsqr = simd_dot(diff, diff)
    
    return DCircumcircle(point1: i, point2: j, point3: k,
                        x: center.x, y: center.y, rsqr: rsqr)
}

@inline(__always)
private func processEdges(_ edges: inout [DEdge: Int], circle: DCircumcircle) {
    edges[DEdge(point1: circle.point1, point2: circle.point2), default: 0] += 1
    edges[DEdge(point1: circle.point2, point2: circle.point3), default: 0] += 1
    edges[DEdge(point1: circle.point3, point2: circle.point1), default: 0] += 1
}

@usableFromInline
@inline(__always)
func processPoint(context: TriangulationContext, currentPoint: DPoint) {
    let collections = context.getReusableCollections()
    var edges = collections.edges
    var trianglesToRemove = collections.toRemove
    
    for circle in context.open {
        let dx = currentPoint.x - circle.x
        let dxSquared = dx * dx
        
        if dx > 0 && dxSquared > circle.rsqr {
            context.completed.insert(circle)
            trianglesToRemove.insert(circle)
            continue
        }
        
        let dy = currentPoint.y - circle.y
        if dxSquared + dy * dy - circle.rsqr <= Double.ulpOfOne {
            trianglesToRemove.insert(circle)
            processEdges(&edges, circle: circle)
        }
    }
    
    context.open.subtract(trianglesToRemove)
    
    for (edge, count) in edges where count == 1 {
        context.open.insert(fastCircumcircle(edge.point1, edge.point2, currentPoint))
    }
    
    context.updateEdges(edges)
    context.updateRemoveSet(trianglesToRemove)
}

// MARK: - Main Triangulation

public func triangulate(_ points: [DPoint]) -> [DTriangle] {
    guard points.count >= 3 else { return [] }
    
    let context = TriangulationContext(capacity: points.count)
    var uniqueCount = 0
    
    // Remove duplicates and collect points
    var seen = Set<DPoint>(minimumCapacity: min(points.count, 1024))
    var sortedPoints = ContiguousArray<DPoint>()
    sortedPoints.reserveCapacity(points.count)
    
    for (index, point) in points.enumerated() {
        let pointWithIndex = DPoint(x: point.x, y: point.y, index: index)
        if seen.insert(pointWithIndex).inserted {
            sortedPoints.append(pointWithIndex)
            uniqueCount += 1
        }
    }
    
    guard uniqueCount >= 3 else { return [] }
    
    // Sort points by x-coordinate
    sortedPoints.sort { $0.x < $1.x }
    
    // Copy sorted points to context
    _ = UnsafeMutableBufferPointer(start: context.points.baseAddress!,
                                  count: uniqueCount)
        .initialize(from: sortedPoints)
    
    // Calculate supertriangle bounds
    let bounds = points.reduce(into: (min: SIMD2<Double>(Double.infinity, Double.infinity),
                                    max: SIMD2<Double>(-Double.infinity, -Double.infinity))) { result, point in
        result.min = simd_min(result.min, SIMD2(point.x, point.y))
        result.max = simd_max(result.max, SIMD2(point.x, point.y))
    }
    
    let delta = bounds.max - bounds.min
    let dmax = max(delta.x, delta.y)
    let mid = (bounds.max + bounds.min) * 0.5
    let margin = dmax * 20
    
    // Add supertriangle points
    context.points[uniqueCount] = DPoint(x: mid.x - margin, y: mid.y - dmax, index: -1)
    context.points[uniqueCount + 1] = DPoint(x: mid.x, y: mid.y + margin, index: -2)
    context.points[uniqueCount + 2] = DPoint(x: mid.x + margin, y: mid.y - dmax, index: -3)
    
    // Initialize with supertriangle
    context.open.insert(fastCircumcircle(context.points[uniqueCount],
                                       context.points[uniqueCount + 1],
                                       context.points[uniqueCount + 2]))
    
    // Process points
    for point in sortedPoints {
        processPoint(context: context, currentPoint: point)
    }
    
    // Add remaining open triangles to completed set
    context.completed.formUnion(context.open)
    
    // Prepare final result
    var result = [DTriangle]()
    result.reserveCapacity(context.completed.count)
    
    // Filter out supertriangle
    for circle in context.completed where circle.point1.index >= 0 && circle.point2.index >= 0 && circle.point3.index >= 0 {
        result.append(DTriangle(point1: circle.point1,
                              point2: circle.point2,
                              point3: circle.point3))
    }
    
    return result
}

private extension Int {
    func nextPowerOf2() -> Int {
        var n = self - 1
        n |= n >> 1
        n |= n >> 2
        n |= n >> 4
        n |= n >> 8
        n |= n >> 16
        return n + 1
    }
}
