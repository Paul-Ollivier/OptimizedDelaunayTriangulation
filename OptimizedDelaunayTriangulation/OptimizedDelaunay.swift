import SwiftUI
import simd

// MARK: - Point Struct

struct Point: Hashable, Codable {
    var vector: SIMD2<Double>
    
    var x: Double {
        get { vector.x }
        set { vector.x = newValue }
    }
    
    var y: Double {
        get { vector.y }
        set { vector.y = newValue }
    }
    
    init(x: Double, y: Double) {
        self.vector = SIMD2<Double>(x, y)
    }
}

// MARK: - Triangulation Struct

struct Triangulation: Hashable, Codable, Identifiable {
    var points: [Point]
    var triangles: [Int]
    var halfEdges: [Int]
    var hull: [Int]
    var numberOfEdges: Int
    var id: UUID = UUID()
    
    init(using delaunay: OptimizedDelaunay, with points: [Point]) {
        self.triangles = delaunay.triangles
        self.halfEdges = delaunay.halfEdges
        self.hull = delaunay.hull
        self.numberOfEdges = delaunay.numberOfEdges
        self.points = points
    }
    
    init() {
        self.triangles = []
        self.halfEdges = []
        self.hull = []
        self.points = []
        self.numberOfEdges = 0
    }
}

// MARK: - JCVDelaunay Struct

final class OptimizedDelaunay {
    static let epsilon: Double = pow(2.0, -53)
    static let orientThreshold: Double = (3.0 + 16.0 * epsilon) * epsilon
    static let circleThreshold: Double = (10.0 + 96.0 * epsilon) * epsilon
    
    var triangles: [Int]
    var halfEdges: [Int]
    var hull: [Int]
    var numberOfEdges: Int
    private var maxPoints: Int
    private var numberOfPoints: Int
    private var maxTriangles: Int
    private var hashSize: Int
    private var hullStartIndex: Int
    private var hullSize: Int
    private var hashFactor: Double
    private var hullTriangles: [Int]
    private var hullPrevious: [Int]
    private var coordinates: [SIMD2<Double>]
    private var center: SIMD2<Double>
    private var edgeStack: [Int]
    
    init(maxPoints: Int) {
        self.maxPoints = maxPoints
        self.numberOfPoints = 0
        self.maxTriangles = max(2 * maxPoints - 5, 0)
        self.numberOfEdges = 0
        
        // Preallocate arrays with required capacity
        self.triangles = [Int]()
        self.triangles.reserveCapacity(3 * maxTriangles)
        
        self.halfEdges = [Int]()
        self.halfEdges.reserveCapacity(3 * maxTriangles)
        
        self.hull = [Int]()
        self.hull.reserveCapacity(maxPoints)
        
        self.hullTriangles = [Int](repeating: -1, count: maxPoints)
        self.hullPrevious = [Int](repeating: -1, count: maxPoints)
        self.coordinates = [SIMD2<Double>](repeating: SIMD2<Double>(0.0, 0.0), count: maxPoints)
        
        // Hash array parameters
        self.hashFactor = 0.0
        self.hashSize = 0
        self.hullStartIndex = 0
        self.hullSize = 0
        self.center = SIMD2<Double>(0.0, 0.0)
        
        self.edgeStack = [Int]()
        self.edgeStack.reserveCapacity(512)
    }
    
    func triangulate(points: [Point]) {
        self.numberOfPoints = points.count
        assert(numberOfPoints <= maxPoints, "Number of points exceeds maxPoints")
        
        self.numberOfEdges = 0
        self.hull.removeAll(keepingCapacity: true)
        self.hullSize = 0
        self.hullStartIndex = 0
        self.hashFactor = Double(numberOfPoints).squareRoot().rounded(.up)
        self.hashSize = Int(hashFactor)
        self.hashFactor *= 0.25
        self.center = SIMD2<Double>(0.0, 0.0)
        self.edgeStack.removeAll(keepingCapacity: true)
        
        for i in 0..<numberOfPoints {
            hullTriangles[i] = -1
            hullPrevious[i] = -1
        }
        
        triangles.removeAll(keepingCapacity: true)
        halfEdges.removeAll(keepingCapacity: true)
        
        // Populate coordinates array
        for i in 0..<numberOfPoints {
            coordinates[i] = points[i].vector
        }
        
        if numberOfPoints == 0 {
            return
        }
        
        var hullNext = [Int](repeating: -1, count: numberOfPoints)
        var hullHash = [Int](repeating: -1, count: hashSize)
        var distances = [Double](repeating: 0.0, count: numberOfPoints)
        var pointIndices = [Int](repeating: 0, count: numberOfPoints)
        
        // Compute bounds
        var minX = Double.infinity
        var maxX = -Double.infinity
        var minY = Double.infinity
        var maxY = -Double.infinity
        
        for i in 0..<numberOfPoints {
            let point = coordinates[i]
            minX = min(minX, point.x)
            minY = min(minY, point.y)
            maxX = max(maxX, point.x)
            maxY = max(maxY, point.y)
            pointIndices[i] = i
        }
        
        // Compute the center point of the bounding box
        let centerPoint = SIMD2<Double>(0.5 * (minX + maxX), 0.5 * (minY + maxY))
        
        // Find the seed point closest to the center of the bounding box
        var seedIndex0 = 0
        var minimumDistance = Double.infinity
        for i in 0..<numberOfPoints {
            let distance = distanceSquared(coordinates[i], centerPoint)
            if distance < minimumDistance {
                seedIndex0 = i
                minimumDistance = distance
            }
        }
        
        let seedPoint0 = coordinates[seedIndex0]
        var seedIndex1 = 0, seedIndex2 = 0
        minimumDistance = Double.infinity
        
        // Find the second seed point
        for i in 0..<numberOfPoints {
            if i == seedIndex0 { continue }
            let distance = distanceSquared(seedPoint0, coordinates[i])
            if distance < minimumDistance && distance > OptimizedDelaunay.epsilon {
                seedIndex1 = i
                minimumDistance = distance
            }
        }
        
        var seedPoint1 = coordinates[seedIndex1]
        var minimumRadius = Double.infinity
        
        // Find the third seed point
        for i in 0..<numberOfPoints {
            if i == seedIndex0 || i == seedIndex1 { continue }
            let radius = circumRadius(seedPoint0, seedPoint1, coordinates[i])
            if radius < minimumRadius {
                seedIndex2 = i
                minimumRadius = radius
            }
        }
        
        var seedPoint2 = coordinates[seedIndex2]
        
        if minimumRadius == Double.infinity {
            for i in 0..<numberOfPoints {
                let deltaX = coordinates[i].x - coordinates[0].x
                distances[i] = isNearZero(x: deltaX) ? coordinates[i].y - coordinates[0].y : deltaX
            }
            pointIndices.sort { distances[$0] < distances[$1] }
            hull.removeAll(keepingCapacity: true)
            hull.reserveCapacity(numberOfPoints)
            hull.append(contentsOf: pointIndices)
            triangles.removeAll(keepingCapacity: true)
            halfEdges.removeAll(keepingCapacity: true)
            return
        }
        
        if !orient(seedPoint0, seedPoint1, seedPoint2) {
            (seedIndex1, seedIndex2) = (seedIndex2, seedIndex1)
            (seedPoint1, seedPoint2) = (seedPoint2, seedPoint1)
        }
        
        center = circumCenter(seedPoint0, seedPoint1, seedPoint2)
        
        for i in 0..<numberOfPoints {
            distances[i] = distanceSquared(coordinates[i], center)
        }
        
        pointIndices.sort { distances[$0] < distances[$1] }
        
        hullStartIndex = seedIndex0
        hullSize = 3
        
        hullNext[seedIndex0] = seedIndex1; hullPrevious[seedIndex2] = seedIndex1
        hullNext[seedIndex1] = seedIndex2; hullPrevious[seedIndex0] = seedIndex2
        hullNext[seedIndex2] = seedIndex0; hullPrevious[seedIndex1] = seedIndex0
        
        hullTriangles[seedIndex0] = 0
        hullTriangles[seedIndex1] = 1
        hullTriangles[seedIndex2] = 2
        
        hullHash[hashKey(point: seedPoint0)] = seedIndex0
        hullHash[hashKey(point: seedPoint1)] = seedIndex1
        hullHash[hashKey(point: seedPoint2)] = seedIndex2
        
        _ = addTriangle(seedIndex0, seedIndex1, seedIndex2, -1, -1, -1)
        
        // Optimized Hull Loop
        var previousX: Double = 0
        var previousY: Double = 0
        var cachedDifferenceX: Double = 0
        var cachedDifferenceY: Double = 0
        
        for (order, currentIndex) in pointIndices.enumerated() {
            let currentPoint = coordinates[currentIndex]
            let currentX = currentPoint.x
            let currentY = currentPoint.y
            
            // Compute and cache the differences only once per iteration
            if order > 0 {
                cachedDifferenceX = currentX - previousX
                cachedDifferenceY = currentY - previousY
            }
            
            // Skip this point if it’s too close to the previous one
            if order > 0 && isNearZero(x: cachedDifferenceX) && isNearZero(x: cachedDifferenceY) {
                continue
            }
            
            // Update previous point values
            previousX = currentX
            previousY = currentY
            
            if currentIndex == seedIndex0 || currentIndex == seedIndex1 || currentIndex == seedIndex2 {
                continue
            }
            
            var startIndex = 0
            let key = hashKey(point: currentPoint)
            
            for j in 0..<hashSize {
                startIndex = hullHash[(key + j) % hashSize]
                if startIndex != -1 && startIndex != hullNext[startIndex] { break }
            }
            
            var edgeIndex = hullPrevious[startIndex]
            var hullIndex = startIndex
            startIndex = edgeIndex
            
            while !orient(currentPoint, coordinates[hullIndex], coordinates[edgeIndex]) {
                if hullIndex == startIndex { break }
                edgeIndex = hullIndex
                hullIndex = hullNext[edgeIndex]
            }
            
            var triangleIndex = addTriangle(edgeIndex, currentIndex, hullIndex, -1, -1, hullTriangles[edgeIndex])
            
            hullTriangles[currentIndex] = legalize(edge: triangleIndex + 2)
            hullTriangles[edgeIndex] = triangleIndex
            
            hullSize += 1
            
            var nextHullIndex = hullIndex
            hullIndex = hullNext[hullIndex]
            
            while orient(currentPoint, coordinates[hullIndex], coordinates[nextHullIndex]) {
                triangleIndex = addTriangle(nextHullIndex, currentIndex, hullIndex, hullTriangles[currentIndex], -1, hullTriangles[nextHullIndex])
                hullTriangles[currentIndex] = legalize(edge: triangleIndex + 2)
                hullNext[nextHullIndex] = nextHullIndex
                hullSize -= 1
                nextHullIndex = hullIndex
                hullIndex = hullNext[nextHullIndex]
            }
            
            if edgeIndex == startIndex {
                hullIndex = hullPrevious[edgeIndex]
                while orient(currentPoint, coordinates[edgeIndex], coordinates[hullIndex]) {
                    triangleIndex = addTriangle(hullIndex, currentIndex, edgeIndex, -1, hullTriangles[edgeIndex], hullTriangles[hullIndex])
                    _ = legalize(edge: triangleIndex + 2)
                    hullTriangles[hullIndex] = triangleIndex
                    hullNext[edgeIndex] = edgeIndex
                    hullSize -= 1
                    edgeIndex = hullIndex
                    hullIndex = hullPrevious[edgeIndex]
                }
            }
            
            hullStartIndex = edgeIndex
            hullPrevious[currentIndex] = edgeIndex; hullNext[edgeIndex] = currentIndex
            hullPrevious[nextHullIndex] = currentIndex; hullNext[currentIndex] = nextHullIndex
            
            hullHash[hashKey(point: currentPoint)] = currentIndex
            hullHash[hashKey(point: coordinates[edgeIndex])] = edgeIndex
        }
        
        // Final hull creation
        hull.removeAll(keepingCapacity: true)
        hull.reserveCapacity(hullSize)
        var edgeIndex = hullStartIndex
        for _ in 0..<hullSize {
            hull.append(edgeIndex)
            edgeIndex = hullNext[edgeIndex]
        }
    }

    
    private func legalize(edge edgeIndex: Int) -> Int {
        var a = edgeIndex, a2 = 0
        
        flipEdge: while true {
            let b = halfEdges[a]
            
            if b == -1 {
                a = edgeStack.popLast() ?? -1
                if a == -1 { break flipEdge }
                continue flipEdge
            }
            
            let a0 = a - a % 3
            a2 = a0 + (a + 2) % 3
            
            let a1 = a0 + (a + 1) % 3
            let b0 = b - b % 3
            let b2 = b0 + (b + 2) % 3
            
            let pointN = triangles[a1]
            let pointI = triangles[a2]
            let pointQ = triangles[a]
            let pointP = triangles[b2]
            
            let isIllegal = inCircumCircle(coordinates[pointN], coordinates[pointI], coordinates[pointQ], coordinates[pointP])
            
            if isIllegal {
                triangles[a] = pointP
                triangles[b] = pointI
                
                let hb2 = halfEdges[b2]
                if hb2 == -1 {
                    hullTriangles[pointP] = a
                }
                
                linkEdges(a, hb2)
                linkEdges(b, halfEdges[a2])
                linkEdges(a2, b2)
                
                edgeStack.append(b0 + (b + 1) % 3)
            } else {
                a = edgeStack.popLast() ?? -1
                if a == -1 { break flipEdge }
            }
        }
        return a2
    }
    
    private func hashKey(point: SIMD2<Double>) -> Int {
        return Int(hashFactor * pseudoAngle(point - center).rounded(.down)) % hashSize
    }
    
    private func linkEdges(_ a: Int, _ b: Int) {
        if a < halfEdges.count {
            halfEdges[a] = b
        } else {
            halfEdges.append(b)
        }
        if b != -1 {
            if b < halfEdges.count {
                halfEdges[b] = a
            } else {
                halfEdges.append(a)
            }
        }
    }
    
    private func addTriangle(_ index0: Int, _ index1: Int, _ index2: Int, _ edgeA: Int, _ edgeB: Int, _ edgeC: Int) -> Int {
        let triangleIndex = triangles.count
        // Ensure capacity before appending
        if triangles.count + 3 > triangles.capacity {
            triangles.reserveCapacity(triangles.capacity + 3 * maxTriangles)
            halfEdges.reserveCapacity(halfEdges.capacity + 3 * maxTriangles)
        }
        triangles.append(contentsOf: [index0, index1, index2])
        
        linkEdges(triangleIndex, edgeA)
        linkEdges(triangleIndex + 1, edgeB)
        linkEdges(triangleIndex + 2, edgeC)
        
        numberOfEdges += 3
        return triangleIndex
    }
    
    // MARK: - Geometry Functions
    
    @inline(__always)
    func distanceSquared(_ a: SIMD2<Double>, _ b: SIMD2<Double>) -> Double {
        return simd_distance_squared(a, b)
    }

    func inCircumCircle(_ a: SIMD2<Double>, _ b: SIMD2<Double>,
                        _ c: SIMD2<Double>, _ p: SIMD2<Double>) -> Bool {
        let ad = a - p
        let bd = b - p
        let cd = c - p

        let ab = ad.x * bd.y - ad.y * bd.x
        let bc = bd.x * cd.y - bd.y * cd.x
        let ca = cd.x * ad.y - cd.y * ad.x

        let alift = simd_length_squared(ad)
        let blift = simd_length_squared(bd)
        let clift = simd_length_squared(cd)

        let determinant = alift * bc + blift * ca + clift * ab

        let permanent = (abs(bd.x * cd.y - bd.y * cd.x) + abs(cd.x * ad.y - cd.y * ad.x)) * alift
                      + (abs(cd.x * ad.y - cd.y * ad.x) + abs(ad.x * bd.y - ad.y * bd.x)) * blift
                      + (abs(ad.x * bd.y - ad.y * bd.x) + abs(bd.x * cd.y - bd.y * cd.x)) * clift
        let errorBound = OptimizedDelaunay.circleThreshold * permanent
        return determinant < -errorBound
    }
    
    func circumRadius(_ a: SIMD2<Double>, _ b: SIMD2<Double>, _ c: SIMD2<Double>) -> Double {
        let ab = b - a
        let ac = c - a

        let abLengthSquared = simd_length_squared(ab)
        let acLengthSquared = simd_length_squared(ac)
        let cross = ab.x * ac.y - ab.y * ac.x
        let crossLength = abs(cross)

        if isNearZero(x: crossLength) {
            return Double.infinity
        }

        let denominator = 0.5 / crossLength

        let x = (ac.y * abLengthSquared - ab.y * acLengthSquared) * denominator
        let y = (ab.x * acLengthSquared - ac.x * abLengthSquared) * denominator

        return x * x + y * y
    }
    
    func circumCenter(_ a: SIMD2<Double>, _ b: SIMD2<Double>, _ c: SIMD2<Double>) -> SIMD2<Double> {
        let d = 2 * (a.x * (b.y - c.y) + b.x * (c.y - a.y) + c.x * (a.y - b.y))
        if isNearZero(x: d) {
            return SIMD2<Double>(Double.nan, Double.nan)
        }

        let asq = simd_length_squared(a)
        let bsq = simd_length_squared(b)
        let csq = simd_length_squared(c)

        let x = (asq * (b.y - c.y) + bsq * (c.y - a.y) + csq * (a.y - b.y)) / d
        let y = (asq * (c.x - b.x) + bsq * (a.x - c.x) + csq * (b.x - a.x)) / d
        return SIMD2<Double>(x, y)
    }
    
    @inline(__always)
    func pseudoAngle(_ vector: SIMD2<Double>) -> Double {
        let p = vector.x / (abs(vector.x) + abs(vector.y))
        return vector.y > 0 ? 3 - p : 1 + p // [0..4]
    }
    
    func orientIfSure(_ a: SIMD2<Double>, _ b: SIMD2<Double>, _ c: SIMD2<Double>) -> Double {
        let detLeft = (a.y - c.y) * (b.x - c.x)
        let detRight = (a.x - c.x) * (b.y - c.y)
        let determinant = detLeft - detRight
        return abs(determinant) >= OptimizedDelaunay.orientThreshold * abs(detLeft + detRight) ? determinant : 0
    }
    
    func orient(_ a: SIMD2<Double>, _ b: SIMD2<Double>, _ c: SIMD2<Double>) -> Bool {
        var determinant = orientIfSure(a, b, c)
        if determinant != 0 { return determinant > 0 }
        determinant = orientIfSure(b, c, a)
        if determinant != 0 { return determinant > 0 }
        determinant = orientIfSure(c, a, b)
        return determinant >= 0
    }

    @inline(__always)
    func isNearZero(x: Double) -> Bool {
        return abs(x) <= OptimizedDelaunay.epsilon
    }
}
