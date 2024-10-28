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
    var open: [DCircumcircle] // Changed from Set to Array
    var completed: Set<DCircumcircle>
    var edgePool: Set<DEdge>
    var totalProcessPointTime: Double

    @usableFromInline
    init(capacity: Int) {
        self.points = []
        self.points.reserveCapacity(capacity + 3)
        self.open = []
        self.completed = Set(minimumCapacity: capacity)
        self.edgePool = Set(minimumCapacity: capacity)
        self.totalProcessPointTime = 0
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

// Binary search extension
extension Array {
    func binarySearch(comparator: (Element) -> ComparisonResult) -> Int {
        var low = 0
        var high = self.count
        while low < high {
            let mid = (low + high) / 2
            switch comparator(self[mid]) {
            case .orderedAscending:
                low = mid + 1
            case .orderedDescending:
                high = mid
            case .orderedSame:
                return mid
            }
        }
        return low
    }
}

@usableFromInline
func insertCircleIntoOpen(_ open: inout [DCircumcircle], _ circle: DCircumcircle) {
    let insertionIndex = open.binarySearch { existingCircle in
        if existingCircle.x < circle.x {
            return .orderedAscending
        } else if existingCircle.x > circle.x {
            return .orderedDescending
        } else {
            return .orderedSame
        }
    }
    open.insert(circle, at: insertionIndex)
}

@usableFromInline
func processPoint(context: TriangulationContext, currentPoint: DPoint) {
    let startTime = DispatchTime.now()
    // Reset the reusable collections
    context.edgePool.removeAll(keepingCapacity: true)

    var edges = context.edgePool
    var removeIndices = [Int]()
    var index = 0

    while index < context.open.count {
        let circle = context.open[index]
        let dx = currentPoint.x - circle.x
        let dxSquared = dx * dx

        // Early exit if the circle cannot contain the point
        if dx > 0 && dxSquared > circle.rsqr {
            context.completed.insert(circle)
            removeIndices.append(index)
            index += 1
            continue
        }

        // Early break if circle is too far to the left
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

    // Remove circles that are no longer needed
    for i in removeIndices.reversed() {
        context.open.remove(at: i)
    }

    // Insert new circumcircles, maintaining the sorted order
    for edge in edges {
        let newCircle = fastCircumcircle(edge.point1, edge.point2, currentPoint)
        insertCircleIntoOpen(&context.open, newCircle)
    }

    // Update the context's edge pool
    context.edgePool = edges

    let endTime = DispatchTime.now()
    // Accumulate time per processPoint call in milliseconds
    context.totalProcessPointTime += Double(endTime.uptimeNanoseconds - startTime.uptimeNanoseconds) * 1e-6
}

// MARK: - Main Triangulation

public func triangulate(_ inputPoints: [DPoint]) -> [DTriangle] {
    var processTimes = [String: Double]()
    let totalStartTime = DispatchTime.now()

    guard inputPoints.count >= 3 else {
        processTimes["Total"] = 0
        print(processTimes)
        return []
    }

    // Measure time for removing duplicates and collecting points
    let startTime1 = DispatchTime.now()
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
    let endTime1 = DispatchTime.now()
    processTimes["Removing Duplicates"] = Double(endTime1.uptimeNanoseconds - startTime1.uptimeNanoseconds) * 1e-6

    guard uniquePoints.count >= 3 else {
        processTimes["Total"] = Double(DispatchTime.now().uptimeNanoseconds - totalStartTime.uptimeNanoseconds) * 1e-6
        print(processTimes)
        return []
    }

    // Measure time for sorting points
    let startTime2 = DispatchTime.now()
    // Sort points by x-coordinate
    uniquePoints.sort { $0.x < $1.x }
    let endTime2 = DispatchTime.now()
    processTimes["Sorting Points"] = Double(endTime2.uptimeNanoseconds - startTime2.uptimeNanoseconds) * 1e-6

    // Measure time for calculating supertriangle bounds
    let startTime3 = DispatchTime.now()
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
    let endTime3 = DispatchTime.now()
    processTimes["Calculating Supertriangle Bounds"] = Double(endTime3.uptimeNanoseconds - startTime3.uptimeNanoseconds) * 1e-6

    // Measure time for creating supertriangle points
    let startTime4 = DispatchTime.now()
    // Create supertriangle points
    let superTrianglePoints = [
        DPoint(x: midx - margin, y: midy - dmax, index: -1),
        DPoint(x: midx, y: midy + margin, index: -2),
        DPoint(x: midx + margin, y: midy - dmax, index: -3)
    ]
    let endTime4 = DispatchTime.now()
    processTimes["Creating Supertriangle Points"] = Double(endTime4.uptimeNanoseconds - startTime4.uptimeNanoseconds) * 1e-6

    // Measure time for combining points
    let startTime5 = DispatchTime.now()
    // Combine uniquePoints and superTrianglePoints to form context.points
    context.points = uniquePoints + superTrianglePoints
    let endTime5 = DispatchTime.now()
    processTimes["Combining Points"] = Double(endTime5.uptimeNanoseconds - startTime5.uptimeNanoseconds) * 1e-6

    // Measure time for initializing with supertriangle
    let startTime6 = DispatchTime.now()
    // Initialize with supertriangle
    let superTriangleStartIndex = context.points.count - 3
    let initialCircle = fastCircumcircle(
        context.points[superTriangleStartIndex],
        context.points[superTriangleStartIndex + 1],
        context.points[superTriangleStartIndex + 2]
    )
    context.open.append(initialCircle)
    let endTime6 = DispatchTime.now()
    processTimes["Initializing Supertriangle"] = Double(endTime6.uptimeNanoseconds - startTime6.uptimeNanoseconds) * 1e-6

    // Measure time for processing points
    let startTime7 = DispatchTime.now()
    // Process points
    for point in uniquePoints {
        processPoint(context: context, currentPoint: point)
    }
    let endTime7 = DispatchTime.now()
    processTimes["Processing Points"] = Double(endTime7.uptimeNanoseconds - startTime7.uptimeNanoseconds) * 1e-6

    // Add the total time spent in processPoint function
    processTimes["ProcessPoint Total Time"] = context.totalProcessPointTime

    // Measure time for collecting open and completed triangles
    let startTime8 = DispatchTime.now()
    // Add remaining open triangles to completed set
    context.completed.formUnion(context.open)
    let endTime8 = DispatchTime.now()
    processTimes["Collecting Triangles"] = Double(endTime8.uptimeNanoseconds - startTime8.uptimeNanoseconds) * 1e-6

    // Measure time for preparing final result
    let startTime9 = DispatchTime.now()
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
    let endTime9 = DispatchTime.now()
    processTimes["Preparing Final Result"] = Double(endTime9.uptimeNanoseconds - startTime9.uptimeNanoseconds) * 1e-6

    let totalEndTime = DispatchTime.now()
    processTimes["Total"] = Double(totalEndTime.uptimeNanoseconds - totalStartTime.uptimeNanoseconds) * 1e-6

    // Print the process times in milliseconds
    print("Process Times (ms):", processTimes)

    return result
}
