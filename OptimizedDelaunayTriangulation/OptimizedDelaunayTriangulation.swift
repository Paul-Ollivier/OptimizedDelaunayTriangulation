import Foundation
import simd

// MARK: - Data Types
public struct DPoint {
    public let x: Double
    public let y: Double
    let index: Int
    
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
    
    init(point1: DPoint, point2: DPoint, point3: DPoint) {
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
    
    func containsPoint(_ point: DPoint) -> Bool {
        let dx = point.x - x
        let dy = point.y - y
        return (dx * dx + dy * dy) <= rsqr
    }
    
    static func == (lhs: DCircumcircle, rhs: DCircumcircle) -> Bool {
        return lhs.point1.index == rhs.point1.index &&
               lhs.point2.index == rhs.point2.index &&
               lhs.point3.index == rhs.point3.index
    }
    
    func hash(into hasher: inout Hasher) {
        hasher.combine(point1.index)
        hasher.combine(point2.index)
        hasher.combine(point3.index)
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

// MARK: - Geometric Calculations
func fastCircumcircle(_ p1: DPoint, _ p2: DPoint, _ p3: DPoint) -> DCircumcircle {
    let v1 = SIMD2<Double>(p1.x, p1.y)
    let v2 = SIMD2<Double>(p2.x, p2.y)
    let v3 = SIMD2<Double>(p3.x, p3.y)
    
    let d = 2 * ((v1.x * (v2.y - v3.y)) +
                 (v2.x * (v3.y - v1.y)) +
                 (v3.x * (v1.y - v2.y)))
    
    if abs(d) < Double.ulpOfOne {
        let center = (v1 + v3) * 0.5
        let diff = v1 - center
        let rsqr = simd_length_squared(diff)
        return DCircumcircle(point1: p1, point2: p2, point3: p3,
                           x: center.x, y: center.y, rsqr: rsqr)
    }
    
    let sq1 = simd_length_squared(v1)
    let sq2 = simd_length_squared(v2)
    let sq3 = simd_length_squared(v3)
    
    let x = ((sq1 * (v2.y - v3.y)) +
             (sq2 * (v3.y - v1.y)) +
             (sq3 * (v1.y - v2.y))) / d
    
    let y = ((sq1 * (v3.x - v2.x)) +
             (sq2 * (p1.x - v3.x)) +
             (sq3 * (v2.x - p1.x))) / d
    
    let center = SIMD2<Double>(x, y)
    let diff = v1 - center
    let rsqr = simd_length_squared(diff)
    
    return DCircumcircle(point1: p1, point2: p2, point3: p3,
                       x: x, y: y, rsqr: rsqr)
}

// MARK: - Worker Actor
actor TriangulationWorker {
    private var open: [DCircumcircle]
    private var completed: Set<DCircumcircle>
    private var edgePool: Set<DEdge>
    private let workerId: Int
    
    init(workerId: Int, capacity: Int) {
        self.workerId = workerId
        self.open = []
        self.completed = Set(minimumCapacity: capacity)
        self.edgePool = Set(minimumCapacity: capacity)
    }
    
    func process(points: [DPoint], initialOpen: [DCircumcircle]) async -> (completed: Set<DCircumcircle>, remaining: [DCircumcircle]) {
        self.open = initialOpen
        
        for point in points {
            processPoint(point)
        }
        
        return (self.completed, self.open)
    }
    
    private func processPoint(_ point: DPoint) {
        edgePool.removeAll(keepingCapacity: true)
        var removeIndices = [Int]()
        var index = 0
        
        while index < open.count {
            let circle = open[index]
            let dx = point.x - circle.x
            let dxSquared = dx * dx
            
            if dx > 0 && dxSquared > circle.rsqr {
                completed.insert(circle)
                removeIndices.append(index)
                index += 1
                continue
            }
            
            if dx < 0 && dxSquared > circle.rsqr {
                break
            }
            
            if circle.containsPoint(point) {
                removeIndices.append(index)
                processEdges(circle)
            }
            
            index += 1
        }
        
        for i in removeIndices.reversed() {
            open.remove(at: i)
        }
        
        for edge in edgePool {
            let newCircle = fastCircumcircle(edge.point1, edge.point2, point)
            insertSorted(circle: newCircle)
        }
    }
    
    private func processEdges(_ circle: DCircumcircle) {
        let edges = [
            DEdge(point1: circle.point1, point2: circle.point2),
            DEdge(point1: circle.point2, point2: circle.point3),
            DEdge(point1: circle.point3, point2: circle.point1)
        ]
        
        for edge in edges {
            if !edgePool.insert(edge).inserted {
                edgePool.remove(edge)
            }
        }
    }
    
    private func insertSorted(circle: DCircumcircle) {
        let index = open.firstIndex { $0.x > circle.x } ?? open.count
        open.insert(circle, at: index)
    }
}

// MARK: - Main Triangulation Function
public func parallelTriangulate(_ inputPoints: [DPoint]) async -> [DTriangle] {
    let startTime = CFAbsoluteTimeGetCurrent()
    guard inputPoints.count >= 3 else { return [] }
    
    // Process points
    var uniquePoints = [DPoint]()
    uniquePoints.reserveCapacity(inputPoints.count)
    var seen = Set<DPoint>(minimumCapacity: inputPoints.count)
    
    for (index, point) in inputPoints.enumerated() {
        let pointWithIndex = DPoint(x: point.x, y: point.y, index: index)
        if seen.insert(pointWithIndex).inserted {
            uniquePoints.append(pointWithIndex)
        }
    }
    
    // Sort points by x-coordinate
    uniquePoints.sort { $0.x < $1.x }
    
    // Calculate bounds for super triangle
    let xmin = uniquePoints.map(\.x).min()!
    let xmax = uniquePoints.map(\.x).max()!
    let ymin = uniquePoints.map(\.y).min()!
    let ymax = uniquePoints.map(\.y).max()!
    
    let dx = xmax - xmin
    let dy = ymax - ymin
    let dmax = max(dx, dy)
    let midx = (xmax + xmin) / 2
    let midy = (ymax + ymin) / 2
    let margin = dmax * 20
    
    // Create super triangle
    let superTrianglePoints = [
        DPoint(x: midx - margin, y: midy - dmax, index: -1),
        DPoint(x: midx, y: midy + margin, index: -2),
        DPoint(x: midx + margin, y: midy - dmax, index: -3)
    ]
    
    let initialCircle = fastCircumcircle(
        superTrianglePoints[0],
        superTrianglePoints[1],
        superTrianglePoints[2]
    )
    
    // Parallel processing setup using TaskGroup
    let availableCores = ProcessInfo.processInfo.activeProcessorCount
    let optimalChunkSize = max(50, uniquePoints.count / (availableCores * 2))
    
    var completedCircles = Set<DCircumcircle>()
    var openCircles = [initialCircle]
    
    await withTaskGroup(of: (Set<DCircumcircle>, [DCircumcircle]).self) { group in
        for chunk in uniquePoints.chunked(into: optimalChunkSize) {
            group.addTask {
                let worker = TriangulationWorker(workerId: 0, capacity: chunk.count)
                return await worker.process(points: chunk, initialOpen: openCircles)
            }
        }
        
        for await (completed, remaining) in group {
            completedCircles.formUnion(completed)
            openCircles.append(contentsOf: remaining)
        }
    }
    
    // Add remaining circles and convert to triangles
    completedCircles.formUnion(openCircles)
    
    let triangles = completedCircles.compactMap { circle -> DTriangle? in
        guard circle.point1.index >= 0 && circle.point2.index >= 0 && circle.point3.index >= 0 else {
            return nil
        }
        return DTriangle(point1: circle.point1, point2: circle.point2, point3: circle.point3)
    }
    
    let endTime = CFAbsoluteTimeGetCurrent()
    print("Triangulation completed in \((endTime - startTime) * 1000)ms")
    print("Generated \(triangles.count) triangles from \(uniquePoints.count) points")
    
    return triangles
}

// Extension to chunk the array into smaller parts
extension Array {
    func chunked(into size: Int) -> [[Element]] {
        stride(from: 0, to: count, by: size).map {
            Array(self[$0..<Swift.min($0 + size, count)])
        }
    }
}
