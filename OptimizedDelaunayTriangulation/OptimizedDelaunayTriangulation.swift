import Foundation
import simd

// MARK: - Public Types

public struct DPoint {
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

extension DPoint: Hashable {
    public func hash(into hasher: inout Hasher) {
        hasher.combine(index)
    }

    public static func == (lhs: DPoint, rhs: DPoint) -> Bool {
        return lhs.index == rhs.index
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
struct DCircumcircle {
    @usableFromInline let point1: DPoint
    @usableFromInline let point2: DPoint
    @usableFromInline let point3: DPoint
    @usableFromInline let x: Double
    @usableFromInline let y: Double
    @usableFromInline let rsqr: Double
}

extension DCircumcircle: Hashable {
    @usableFromInline func hash(into hasher: inout Hasher) {
        hasher.combine(point1.index)
        hasher.combine(point2.index)
        hasher.combine(point3.index)
    }

    @usableFromInline static func == (lhs: DCircumcircle, rhs: DCircumcircle) -> Bool {
        return lhs.point1.index == rhs.point1.index &&
               lhs.point2.index == rhs.point2.index &&
               lhs.point3.index == rhs.point3.index
    }
}

@usableFromInline
struct DEdge {
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

extension DEdge: Hashable {
    @usableFromInline func hash(into hasher: inout Hasher) {
        hasher.combine(point1.index)
        hasher.combine(point2.index)
    }

    @usableFromInline static func == (lhs: DEdge, rhs: DEdge) -> Bool {
        return lhs.point1.index == rhs.point1.index &&
               lhs.point2.index == rhs.point2.index
    }
}

// MARK: - Triangulation Context

@usableFromInline
final class TriangulationContext {
    var points: [DPoint]
    var open: Set<DCircumcircle>
    var completed: Set<DCircumcircle>
    var edgePool: Set<DEdge>
    var removePool: Set<DCircumcircle>

    @usableFromInline
    init(capacity: Int) {
        self.points = []
        self.points.reserveCapacity(capacity + 3)
        self.open = Set(minimumCapacity: capacity)
        self.completed = Set(minimumCapacity: capacity)
        self.edgePool = Set(minimumCapacity: capacity)
        self.removePool = Set(minimumCapacity: capacity)
    }
}

// MARK: - Geometric Calculations

@usableFromInline
func fastCircumcircle(_ i: DPoint, _ j: DPoint, _ k: DPoint) -> DCircumcircle {
    let p1 = SIMD2<Double>(i.x, i.y)
    let p2 = SIMD2<Double>(j.x, j.y)
    let p3 = SIMD2<Double>(k.x, k.y)

    let d = 2 * ((p1.x * (p2.y - p3.y)) +
                 (p2.x * (p3.y - p1.y)) +
                 (p3.x * (p1.y - p2.y)))

    if abs(d) < Double.ulpOfOne {
        let center = (p1 + p3) * 0.5
        let diff = p1 - center
        let rsqr = simd_length_squared(diff)
        return DCircumcircle(point1: i, point2: j, point3: k,
                             x: center.x, y: center.y, rsqr: rsqr)
    }

    let sq1 = simd_length_squared(p1)
    let sq2 = simd_length_squared(p2)
    let sq3 = simd_length_squared(p3)

    let x = ((sq1 * (p2.y - p3.y)) +
             (sq2 * (p3.y - p1.y)) +
             (sq3 * (p1.y - p2.y))) / d

    let y = ((sq1 * (p3.x - p2.x)) +
             (sq2 * (p1.x - p3.x)) +
             (sq3 * (p2.x - p1.x))) / d

    let center = SIMD2<Double>(x, y)
    let diff = p1 - center
    let rsqr = simd_length_squared(diff)

    return DCircumcircle(point1: i, point2: j, point3: k,
                         x: center.x, y: center.y, rsqr: rsqr)
}

@inline(__always)
private func processEdges(_ edges: inout Set<DEdge>, circle: DCircumcircle) {
    toggleEdge(&edges, DEdge(point1: circle.point1, point2: circle.point2))
    toggleEdge(&edges, DEdge(point1: circle.point2, point2: circle.point3))
    toggleEdge(&edges, DEdge(point1: circle.point3, point2: circle.point1))
}

@inline(__always)
private func toggleEdge(_ edges: inout Set<DEdge>, _ edge: DEdge) {
    if !edges.insert(edge).inserted {
        edges.remove(edge)
    }
}

@usableFromInline
func processPoint(context: TriangulationContext, currentPoint: DPoint) {
    // Reset the reusable collections
    context.edgePool.removeAll(keepingCapacity: true)
    context.removePool.removeAll(keepingCapacity: true)

    // Use references to the collections for convenience
    var edges = context.edgePool
    var trianglesToRemove = context.removePool

    for circle in context.open {
        let dx = currentPoint.x - circle.x
        let dxSquared = dx * dx

        // Early exit if the circle cannot contain the point
        if dx > 0 && dxSquared > circle.rsqr {
            context.completed.insert(circle)
            trianglesToRemove.insert(circle)
            continue
        }

        let dy = currentPoint.y - circle.y
        let distanceSquared = dxSquared + dy * dy

        if distanceSquared <= circle.rsqr {
            trianglesToRemove.insert(circle)
            processEdges(&edges, circle: circle)
        }
    }

    context.open.subtract(trianglesToRemove)

    for edge in edges {
        context.open.insert(fastCircumcircle(edge.point1, edge.point2, currentPoint))
    }

    // Update the context's edge pool
    context.edgePool = edges
    context.removePool = trianglesToRemove
}

// MARK: - Main Triangulation

public func triangulate(_ inputPoints: [DPoint]) -> [DTriangle] {
    guard inputPoints.count >= 3 else { return [] }

    let context = TriangulationContext(capacity: inputPoints.count)
    var uniquePoints = [DPoint]()
    var index = 0

    // Remove duplicates and collect points
    var seen = Set<DPoint>(minimumCapacity: inputPoints.count)
    for point in inputPoints {
        let pointWithIndex = DPoint(x: point.x, y: point.y, index: index)
        if seen.insert(pointWithIndex).inserted {
            uniquePoints.append(pointWithIndex)
            index += 1
        }
    }

    guard uniquePoints.count >= 3 else { return [] }

    // Sort points by x-coordinate
    uniquePoints.sort { $0.x < $1.x }

    // Calculate supertriangle bounds
    let xs = uniquePoints.map { $0.x }
    let ys = uniquePoints.map { $0.y }
    let xmin = xs.min()!
    let xmax = xs.max()!
    let ymin = ys.min()!
    let ymax = ys.max()!

    let dx = xmax - xmin
    let dy = ymax - ymin
    let dmax = max(dx, dy)
    let midx = (xmax + xmin) / 2
    let midy = (ymax + ymin) / 2
    let margin = dmax * 20

    // Create supertriangle points
    let superTrianglePoints = [
        DPoint(x: midx - margin, y: midy - dmax, index: -1),
        DPoint(x: midx, y: midy + margin, index: -2),
        DPoint(x: midx + margin, y: midy - dmax, index: -3)
    ]

    // Combine uniquePoints and superTrianglePoints to form context.points
    context.points = uniquePoints + superTrianglePoints

    // Initialize with supertriangle
    let superTriangleStartIndex = context.points.count - 3
    context.open.insert(fastCircumcircle(
        context.points[superTriangleStartIndex],
        context.points[superTriangleStartIndex + 1],
        context.points[superTriangleStartIndex + 2]
    ))

    // Process points
    for point in uniquePoints {
        processPoint(context: context, currentPoint: point)
    }

    // Add remaining open triangles to completed set
    context.completed.formUnion(context.open)

    // Prepare final result
    var result = [DTriangle]()
    result.reserveCapacity(context.completed.count)

    // Filter out triangles that include supertriangle points
    for circle in context.completed {
        if circle.point1.index >= 0 && circle.point2.index >= 0 && circle.point3.index >= 0 {
            result.append(DTriangle(
                point1: circle.point1,
                point2: circle.point2,
                point3: circle.point3
            ))
        }
    }

    return result
}
