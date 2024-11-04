import SwiftUI
import simd

// MARK: - Static Constants

struct Static {
    static let Epsilon: Double = pow(2.0, -53)
    static let orientThreshold: Double = (3.0 + 16.0 * Epsilon) * Epsilon
    static let circleThreshold: Double = (10.0 + 96.0 * Epsilon) * Epsilon
}

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
    var numberEdges: Int
    var id: UUID? = UUID()

    init(using delaunay: JCVDelaunay, with points: [Point]) {
        self.triangles = delaunay.triangles
        self.halfEdges = delaunay.halfEdges
        self.hull = delaunay.hull
        self.numberEdges = delaunay.numberEdges
        self.points = points
    }

    init() {
        self.triangles = []
        self.halfEdges = []
        self.hull = []
        self.points = []
        self.numberEdges = 0
    }
}

// MARK: - JCVDelaunay Struct

struct JCVDelaunay {
    var triangles: [Int]
    var halfEdges: [Int]
    var hull: [Int]
    var numberEdges: Int
    private var maxPoints: Int
    private var numberPoints: Int
    private var maxTriangles: Int
    private var hashSize: Int
    private var hullStart: Int
    private var hullSize: Int
    private var hashFactor: Double
    private var hullTri: [Int]
    private var hullPrev: [Int]
    private var centre: SIMD2<Double>
    private var coords: [SIMD2<Double>]
    var edgeStack: [Int]
    
    init(maxPoints: Int) {
        self.maxPoints = maxPoints
        self.numberPoints = 0
        self.maxTriangles = max(2 * maxPoints - 5, 0)
        self.numberEdges = 0
        
        // Preallocate capacity for triangles and halfEdges
        self.triangles = [Int](repeating: 0, count: 3 * maxTriangles)
        self.halfEdges = [Int](repeating: -1, count: 3 * maxTriangles)
        
        // Preallocate hullTri and hullPrev
        self.hullTri = [Int](repeating: -1, count: maxPoints)
        self.hullPrev = [Int](repeating: -1, count: maxPoints)
        self.coords = [SIMD2<Double>](repeating: SIMD2<Double>(0.0, 0.0), count: maxPoints)
        
        // Hash array parameters
        self.hashFactor = 0.0
        self.hashSize = 0
        self.hullStart = 0
        self.hullSize = 0
        self.hull = []
        self.centre = SIMD2<Double>(0.0, 0.0)
        
        // Initialize edgeStack and preallocate capacity
        self.edgeStack = [Int]()
        self.edgeStack.reserveCapacity(512) // Initial capacity, adjust as needed
    }
    
    mutating func triangulate(points: [Point]) {
        self.numberPoints = points.count
        assert(numberPoints <= maxPoints, "Number of points exceeds maxPoints")
        
        // Re-initialize variables that depend on numberPoints
        self.numberEdges = 0
        self.hull.removeAll(keepingCapacity: true)
        self.hullSize = 0
        self.hullStart = 0
        self.hashFactor = Double(numberPoints).squareRoot().rounded(.up)
        self.hashSize = Int(hashFactor)
        self.hashFactor *= 0.25
        self.centre = SIMD2<Double>(0.0, 0.0)
        self.edgeStack.removeAll(keepingCapacity: true)
        
        // Reset hullTri and hullPrev arrays up to maxPoints
        for i in 0..<maxPoints {
            hullTri[i] = -1
            hullPrev[i] = -1
        }
        
        // Reset triangles and halfEdges arrays
        triangles.removeAll(keepingCapacity: true)
        halfEdges.removeAll(keepingCapacity: true)
        
        // Populate coords array
        for i in 0..<numberPoints {
            coords[i] = points[i].vector
        }
        
        if numberPoints == 0 {
            return
        }
        
        var hullNext = [Int](repeating: -1, count: numberPoints)
        var hullHash = [Int](repeating: -1, count: hashSize)
        
        var dists = [Double](repeating: 0.0, count: numberPoints)
        var ids = [Int](repeating: 0, count: numberPoints)
        
        // Compute bounds and populate coords array
        var minX = Double.infinity
        var maxX = -Double.infinity
        var minY = Double.infinity
        var maxY = -Double.infinity
        
        for (i, _) in points.enumerated() {
            let p = coords[i]
            if p.x < minX { minX = p.x }
            if p.y < minY { minY = p.y }
            if p.x > maxX { maxX = p.x }
            if p.y > maxY { maxY = p.y }
            ids[i] = i
        }
        
        // Find seed points
        let c = SIMD2<Double>(0.5 * (minX + maxX), 0.5 * (minY + maxY))
        
        var i0 = 0, i1 = 0, i2 = 0
        
        var minDist = Double.infinity
        for i in 0..<numberPoints {
            let d = distanceSquared(coords[i], c)
            if d < minDist {
                i0 = i
                minDist = d
            }
        }
        let i0v = coords[i0]
        
        minDist = Double.infinity
        for i in 0..<numberPoints {
            if i == i0 { continue }
            let d = distanceSquared(coords[i0], coords[i])
            if d < minDist && d > Static.Epsilon {
                i1 = i
                minDist = d
            }
        }
        var i1v = coords[i1]
        
        var minRadius = Double.infinity
        for i in 0..<numberPoints {
            if i == i0 || i == i1 { continue }
            let r = circumRadius(coords[i0], coords[i1], coords[i])
            if r < minRadius {
                i2 = i
                minRadius = r
            }
        }
        var i2v = coords[i2]
        
        if minRadius == Double.infinity {
            for i in 0..<numberPoints {
                let deltaX = coords[i].x - coords[0].x
                dists[i] = isNearZero(x: deltaX) ? coords[i].y - coords[0].y : deltaX
            }
            ids.sort { dists[$0] < dists[$1] }
            hull = ids
            hull.removeLast(numberPoints - hull.count)
            triangles.removeAll(keepingCapacity: true)
            halfEdges.removeAll(keepingCapacity: true)
            return
        }
        
        if !orient(i0v, i1v, i2v) {
            (i1, i2) = (i2, i1)
            (i1v, i2v) = (i2v, i1v)
        }
        
        centre = circumCentre(coords[i0], coords[i1], coords[i2])
        
        for i in 0..<numberPoints {
            dists[i] = distanceSquared(coords[i], centre)
        }
        
        ids.sort { dists[$0] < dists[$1] }
        
        hullStart = i0
        hullSize = 3
        
        hullNext[i0] = i1; hullPrev[i2] = i1
        hullNext[i1] = i2; hullPrev[i0] = i2
        hullNext[i2] = i0; hullPrev[i1] = i0
        
        hullTri[i0] = 0
        hullTri[i1] = 1
        hullTri[i2] = 2
        
        hullHash[hashKey(p: coords[i0])] = i0
        hullHash[hashKey(p: coords[i1])] = i1
        hullHash[hashKey(p: coords[i2])] = i2
        
        _ = addTriangle(i0, i1, i2, -1, -1, -1)
        
        var xp: Double = 0, yp: Double = 0
        
        for (k, i) in ids.enumerated() {
            let v = coords[i]
            let x = v.x
            let y = v.y
            
            if k > 0 && isNearZero(x: x - xp) && isNearZero(x: y - yp) { continue }
            
            xp = x
            yp = y
            
            if i == i0 || i == i1 || i == i2 { continue }
            
            var start = 0
            let key = hashKey(p: v)
            
            for j in 0..<hashSize {
                start = hullHash[(key + j) % hashSize]
                if start != -1 && start != hullNext[start] { break }
            }
            
            var e = hullPrev[start]
            var q = start
            start = e
            
            while !orient(v, coords[q], coords[e]) {
                if q == start { break }
                e = q
                q = hullNext[e]
            }
            
            var t = addTriangle(e, i, q, -1, -1, hullTri[e])
            
            hullTri[i] = legalize(edge: t + 2)
            hullTri[e] = t
            
            hullSize += 1
            
            var n = q
            q = hullNext[q]
            
            while orient(v, coords[q], coords[n]) {
                
                t = addTriangle(n, i, q, hullTri[i], -1, hullTri[n])
                hullTri[i] = legalize(edge: t + 2)
                hullNext[n] = n
                hullSize -= 1
                n = q
                q = hullNext[n]
            }
            
            if e == start {
                q = hullPrev[e]
                while orient(v, coords[e], coords[q]) {
                    
                    t = addTriangle(q, i, e, -1, hullTri[e], hullTri[q])
                    _ = legalize(edge: t + 2)
                    hullTri[q] = t
                    hullNext[e] = e
                    hullSize -= 1
                    e = q
                    q = hullPrev[e]
                }
            }
            
            hullStart = e
            hullPrev[i] = e; hullNext[e] = i
            hullPrev[n] = i; hullNext[i] = n
            
            hullHash[hashKey(p: v)] = i
            hullHash[hashKey(p: coords[e])] = e
        }
        
        hull = [Int](repeating: 0, count: hullSize)
        var e = hullStart
        for i in 0..<hullSize {
            hull[i] = e
            e = hullNext[e]
        }
        
        // Triangles and halfEdges are already correct
    }
    
    private mutating func legalize(edge e: Int) -> Int {
        var a = e, a2 = 0
        
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
            
            let n = triangles[a1]
            let i = triangles[a2]
            let q = triangles[a]
            let p = triangles[b2]
            
            let illegal = inCircumCircle(coords[n], coords[i], coords[q], coords[p])
            
            if illegal {
                triangles[a] = p
                triangles[b] = i
                
                let hb2 = halfEdges[b2]
                if hb2 == -1 {
                    hullTri[p] = a
                }
                
                link(a, hb2)
                link(b, halfEdges[a2])
                link(a2, b2)
                
                edgeStack.append(b0 + (b + 1) % 3)
            } else {
                a = edgeStack.popLast() ?? -1
                if a == -1 { break flipEdge }
            }
        }
        return a2
    }
    
    private func hashKey(p: SIMD2<Double>) -> Int {
        return Int(hashFactor * pseudoAngle(p - centre).rounded(.down)) % hashSize
    }
    
    private mutating func link(_ a: Int, _ b: Int) {
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
    
    private mutating func addTriangle(_ i0: Int, _ i1: Int, _ i2: Int,
                                      _ a: Int, _ b: Int, _ c: Int) -> Int {
        let t = triangles.count
        triangles.append(contentsOf: [i0, i1, i2])
        
        link(t, a)
        link(t + 1, b)
        link(t + 2, c)
        
        numberEdges += 3
        return t
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

        let det = alift * bc + blift * ca + clift * ab

        let permanent = (abs(bd.x * cd.y - bd.y * cd.x) + abs(cd.x * ad.y - cd.y * ad.x)) * alift
                      + (abs(cd.x * ad.y - cd.y * ad.x) + abs(ad.x * bd.y - ad.y * bd.x)) * blift
                      + (abs(ad.x * bd.y - ad.y * bd.x) + abs(bd.x * cd.y - bd.y * cd.x)) * clift
        let errbound = Static.circleThreshold * permanent
        return det < -errbound
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
    
    func circumCentre(_ a: SIMD2<Double>, _ b: SIMD2<Double>, _ c: SIMD2<Double>) -> SIMD2<Double> {
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
    func pseudoAngle(_ v: SIMD2<Double>) -> Double {
        let p = v.x / (abs(v.x) + abs(v.y))
        return v.y > 0 ? 3 - p : 1 + p // [0..4]
    }
    
    func orientIfSure(_ a: SIMD2<Double>, _ b: SIMD2<Double>, _ c: SIMD2<Double>) -> Double {
        let detLeft = (a.y - c.y) * (b.x - c.x)
        let detRight = (a.x - c.x) * (b.y - c.y)
        let det = detLeft - detRight
        return abs(det) >= Static.orientThreshold * abs(detLeft + detRight) ? det : 0
    }
    
    func orient(_ a: SIMD2<Double>, _ b: SIMD2<Double>, _ c: SIMD2<Double>) -> Bool {
        var det = orientIfSure(a, b, c)
        if det != 0 { return det > 0 }
        det = orientIfSure(b, c, a)
        if det != 0 { return det > 0 }
        det = orientIfSure(c, a, b)
        return det >= 0
    }

    @inline(__always)
    func isNearZero(x: Double) -> Bool {
        return abs(x) <= Static.Epsilon
    }
}
