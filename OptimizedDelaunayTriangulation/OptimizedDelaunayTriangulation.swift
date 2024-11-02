import Foundation
import simd

// Public types for the API
struct Point {
    var x: Double
    var y: Double
}

struct Delaunator {
    let triangles: [UInt32]
    
    static func from(points: [Point]) -> Delaunator {
        // Convert Point array to DPoint array
        let dpoints = points.enumerated().map { index, point in
            DPoint(x: point.x, y: point.y, index: index)
        }
        
        // Get triangulation
        let triangulation = triangulate(dpoints)
        
        // Convert triangulation to array of indices
        var triangleIndices: [UInt32] = []
        triangleIndices.reserveCapacity(triangulation.count * 3)
        
        for triangle in triangulation {
            triangleIndices.append(UInt32(triangle.point1.index))
            triangleIndices.append(UInt32(triangle.point2.index))
            triangleIndices.append(UInt32(triangle.point3.index))
        }
        
        return Delaunator(triangles: triangleIndices)
    }
}

// Internal types for triangulation
private struct DPoint {
    let x: Double
    let y: Double
    let index: Int
    
    init(x: Double, y: Double, index: Int) {
        self.x = x
        self.y = y
        self.index = index
    }
}

private struct DTriangle {
    let point1: DPoint
    let point2: DPoint
    let point3: DPoint
}

private struct DCircumcircle {
    let point1: DPoint
    let point2: DPoint
    let point3: DPoint
    let x: Double
    let y: Double
    let rsqr: Double
}

private struct DEdge {
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

// Hashable conformance
extension DPoint: Hashable {
    func hash(into hasher: inout Hasher) {
        hasher.combine(index)
    }
    
    static func == (lhs: DPoint, rhs: DPoint) -> Bool {
        return lhs.index == rhs.index
    }
}

extension DCircumcircle: Hashable {
    func hash(into hasher: inout Hasher) {
        hasher.combine(point1.index)
        hasher.combine(point2.index)
        hasher.combine(point3.index)
    }
    
    static func == (lhs: DCircumcircle, rhs: DCircumcircle) -> Bool {
        return lhs.point1.index == rhs.point1.index &&
               lhs.point2.index == rhs.point2.index &&
               lhs.point3.index == rhs.point3.index
    }
}

extension DEdge: Hashable {
    func hash(into hasher: inout Hasher) {
        hasher.combine(point1.index)
        hasher.combine(point2.index)
    }
    
    static func == (lhs: DEdge, rhs: DEdge) -> Bool {
        return lhs.point1.index == rhs.point1.index &&
               lhs.point2.index == rhs.point2.index
    }
}

// Triangulation context
private final class TriangulationContext {
    var points: [DPoint]
    var open: [DCircumcircle]
    var completed: Set<DCircumcircle>
    var edgePool: Set<DEdge>
    
    init(capacity: Int) {
        self.points = []
        self.points.reserveCapacity(capacity + 3)
        self.open = []
        self.completed = Set(minimumCapacity: capacity)
        self.edgePool = Set(minimumCapacity: capacity)
    }
}

// Geometric calculations
private func fastCircumcircle(_ i: DPoint, _ j: DPoint, _ k: DPoint) -> DCircumcircle {
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

private func processEdges(_ edges: inout Set<DEdge>, circle: DCircumcircle) {
    toggleEdge(&edges, DEdge(point1: circle.point1, point2: circle.point2))
    toggleEdge(&edges, DEdge(point1: circle.point2, point2: circle.point3))
    toggleEdge(&edges, DEdge(point1: circle.point3, point2: circle.point1))
}

private func toggleEdge(_ edges: inout Set<DEdge>, _ edge: DEdge) {
    if !edges.insert(edge).inserted {
        edges.remove(edge)
    }
}

private func insertCircleIntoOpen(_ open: inout [DCircumcircle], _ circle: DCircumcircle) {
    var low = 0
    var high = open.count
    
    while low < high {
        let mid = (low + high) / 2
        if open[mid].x < circle.x {
            low = mid + 1
        } else {
            high = mid
        }
    }
    
    open.insert(circle, at: low)
}

private func processPoint(context: TriangulationContext, currentPoint: DPoint) {
    context.edgePool.removeAll(keepingCapacity: true)
    var edges = context.edgePool
    var removeIndices = [Int]()
    var index = 0
    
    while index < context.open.count {
        let circle = context.open[index]
        let dx = currentPoint.x - circle.x
        let dxSquared = dx * dx
        
        if dx > 0 && dxSquared > circle.rsqr {
            context.completed.insert(circle)
            removeIndices.append(index)
            index += 1
            continue
        }
        
        if dx < 0 && dxSquared > circle.rsqr {
            break
        }
        
        let dy = currentPoint.y - circle.y
        let distanceSquared = dxSquared + dy * dy
        
        if distanceSquared <= circle.rsqr {
            removeIndices.append(index)
            processEdges(&edges, circle: circle)
        }
        
        index += 1
    }
    
    for i in removeIndices.reversed() {
        context.open.remove(at: i)
    }
    
    for edge in edges {
        let newCircle = fastCircumcircle(edge.point1, edge.point2, currentPoint)
        insertCircleIntoOpen(&context.open, newCircle)
    }
    
    context.edgePool = edges
}

private func triangulate(_ inputPoints: [DPoint]) -> [DTriangle] {
    guard inputPoints.count >= 3 else { return [] }
    
    let context = TriangulationContext(capacity: inputPoints.count)
    var uniquePoints = [DPoint]()
    var seen = Set<DPoint>(minimumCapacity: inputPoints.count)
    
    for point in inputPoints {
        if seen.insert(point).inserted {
            uniquePoints.append(point)
        }
    }
    
    guard uniquePoints.count >= 3 else { return [] }
    
    uniquePoints.sort { $0.x < $1.x }
    
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
    
    let superTrianglePoints = [
        DPoint(x: midx - margin, y: midy - dmax, index: -1),
        DPoint(x: midx, y: midy + margin, index: -2),
        DPoint(x: midx + margin, y: midy - dmax, index: -3)
    ]
    
    context.points = uniquePoints + superTrianglePoints
    
    let superTriangleStartIndex = context.points.count - 3
    let initialCircle = fastCircumcircle(
        context.points[superTriangleStartIndex],
        context.points[superTriangleStartIndex + 1],
        context.points[superTriangleStartIndex + 2]
    )
    context.open.append(initialCircle)
    
    for point in uniquePoints {
        processPoint(context: context, currentPoint: point)
    }
    
    context.completed.formUnion(context.open)
    
    var result = [DTriangle]()
    result.reserveCapacity(context.completed.count)
    
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
