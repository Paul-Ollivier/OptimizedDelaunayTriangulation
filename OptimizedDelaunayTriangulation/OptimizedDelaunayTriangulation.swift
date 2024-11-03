import Foundation
import Darwin
import simd

struct Point {
    let coords: SIMD2<Double>
    
    init(_ coords: SIMD2<Double>)
    {
        self.coords = coords
    }
}

class Delaunator {
    private let EPSILON: Double = Darwin.pow(2.0, -52)
    private var EDGE_STACK = [UInt](repeating: 0, count: 512)
    
    private var coords: [Double]
    public private(set) var triangles: [UInt]
    public private(set) var halfedges: [Int]
    private var trianglesLen: Int = 0
    public private(set) var hull: [UInt] = []
    
    // Temporary arrays
    private var _hashSize: Int
    private var _hullPrev: [UInt]
    private var _hullNext: [UInt]
    private var _hullTri: [UInt]
    private var _hullHash: [Int]
    private var _ids: [UInt]
    private var _dists: [Double]
    private var _center = SIMD2<Double>(0, 0)
    private var _hullStart: UInt = 0
    
    static func from(points: [Point]) -> Delaunator {
        let n = points.count
        var coords = [Double](repeating: 0, count: n * 2)
        
        for i in 0..<n {
            coords[2 * i] = points[i].coords.x
            coords[2 * i + 1] = points[i].coords.y
        }
        
        return Delaunator(coords: coords)
    }
    
    init(coords: [Double]) {
        let n = coords.count >> 1
        self.coords = coords
        
        let maxTriangles = max(2 * n - 5, 0)
        self.triangles = [UInt](repeating: 0, count: maxTriangles * 3)
        self.halfedges = [Int](repeating: 0, count: maxTriangles * 3)
        
        self._hashSize = Int(ceil(sqrt(Double(n))))
        self._hullPrev = [UInt](repeating: 0, count: n)
        self._hullNext = [UInt](repeating: 0, count: n)
        self._hullTri = [UInt](repeating: 0, count: n)
        self._hullHash = [Int](repeating: -1, count: self._hashSize)
        self._ids = [UInt](repeating: 0, count: n)
        self._dists = [Double](repeating: 0, count: n)
        
        update()
    }
    
    private func update() {
        let n = coords.count >> 1
        
        var minX = Double.infinity
        var minY = Double.infinity
        var maxX = -Double.infinity
        var maxY = -Double.infinity
        
        for i in 0..<n {
            let x = coords[2 * i]
            let y = coords[2 * i + 1]
            if x < minX { minX = x }
            if y < minY { minY = y }
            if x > maxX { maxX = x }
            if y > maxY { maxY = y }
            _ids[i] = UInt(i)
        }
        
        let c = SIMD2<Double>((minX + maxX) / 2, (minY + maxY) / 2)
        
        var i0: Int = 0
        var minDist = Double.infinity
        
        for i in 0..<n {
            let d = squaredDistance(c, SIMD2<Double>(coords[2 * i], coords[2 * i + 1]))
            if d < minDist {
                i0 = i
                minDist = d
            }
        }
        
        let i0p = SIMD2<Double>(coords[2 * i0], coords[2 * i0 + 1])
        
        var i1: Int = 0
        minDist = Double.infinity
        
        for i in 0..<n {
            if i == i0 { continue }
            let d = squaredDistance(i0p, SIMD2<Double>(coords[2 * i], coords[2 * i + 1]))
            if d < minDist && d > 0 {
                i1 = i
                minDist = d
            }
        }
        
        var i1p = SIMD2<Double>(coords[2 * i1], coords[2 * i1 + 1])
        
        var minRadius = Double.infinity
        var i2: Int = 0
        
        for i in 0..<n {
            if i == i0 || i == i1 { continue }
            let r = circumradius(i0p, i1p, SIMD2<Double>(coords[2 * i], coords[2 * i + 1]))
            if r < minRadius {
                i2 = i
                minRadius = r
            }
        }
        
        var i2p = SIMD2<Double>(coords[2 * i2], coords[2 * i2 + 1])
        
        if minRadius == Double.infinity {
            // Handle collinear points
            for i in 0..<n {
                _dists[i] = coords[2 * i] - coords[0] != 0 ?
                coords[2 * i] - coords[0] : coords[2 * i + 1] - coords[1]
            }
            quicksort(ids: &_ids, dists: _dists, left: 0, right: n - 1)
            var hull = [UInt](repeating: 0, count: n)
            var j = 0
            var d0 = -Double.infinity
            
            for i in 0..<n {
                let id = _ids[i]
                if _dists[Int(id)] > d0 {
                    hull[j] = id
                    j += 1
                    d0 = _dists[Int(id)]
                }
            }
            
            self.hull = Array(hull[0..<j])
            self.triangles = []
            self.halfedges = []
            trianglesLen = 0
            return
        }
        
        if orient2d(i0p, i1p, i2p) < 0 {
            let i = i1
            let p = i1p
            i1 = i2
            i1p = i2p
            i2 = i
            i2p = p
        }
        
        let center = circumcenter(i0p, i1p, i2p)
        self._center = center.coords
        
        for i in 0..<n {
            _dists[i] = squaredDistance(SIMD2<Double>(coords[2 * i], coords[2 * i + 1]), center.coords)
        }
        
        quicksort(ids: &_ids, dists: _dists, left: 0, right: n - 1)
        
        self._hullStart = UInt(i0)
        var hullSize = 3
        
        _hullNext[i0] = UInt(i1)
        _hullPrev[i2] = UInt(i1)
        _hullNext[i1] = UInt(i2)
        _hullPrev[i0] = UInt(i2)
        _hullNext[i2] = UInt(i0)
        _hullPrev[i1] = UInt(i0)
        
        _hullTri[i0] = 0
        _hullTri[i1] = 1
        _hullTri[i2] = 2
        
        _hullHash = [Int](repeating: -1, count: _hashSize)
        _hullHash[hashKey(i0p)] = i0
        _hullHash[hashKey(i1p)] = i1
        _hullHash[hashKey(i2p)] = i2
        
        trianglesLen = 0
        _ = addTriangle(i0: UInt(i0), i1: UInt(i1), i2: UInt(i2),
                        a: -1, b: -1, c: -1)
        
        var xp: Double = 0
        var yp: Double = 0
        
        for k in 0..<_ids.count {
            let i = Int(_ids[k])
            let p = SIMD2<Double>(coords[2 * i], coords[2 * i + 1])
            
            if k > 0 && abs(p.x - xp) <= EPSILON && abs(p.y - yp) <= EPSILON { continue }
            xp = p.x
            yp = p.y
            
            if i == i0 || i == i1 || i == i2 { continue }
            
            var start: Int = 0
            let key = hashKey(p)
            
            for j in 0..<_hashSize {
                start = Int(_hullHash[(key + j) % _hashSize])
                if start != -1 && start != Int(_hullNext[start]) { break }
            }
            
            start = Int(_hullPrev[start])
            var e = start
            
            while true {
                let q = Int(_hullNext[e])
                if orient2d(p,
                            SIMD2<Double>(coords[2 * e], coords[2 * e + 1]),
                            SIMD2<Double>(coords[2 * q], coords[2 * q + 1])) < 0 { break }
                e = q
                if e == start {
                    e = -1
                    break
                }
            }
            
            if e == -1 { continue }
            
            var t = addTriangle(i0: UInt(e), i1: UInt(i),
                                i2: _hullNext[e],
                                a: -1, b: -1, c: Int32(_hullTri[e]))
            
            _hullTri[i] = legalize(a: t + 2)
            _hullTri[e] = UInt(t)
            hullSize += 1
            
            var n = Int(_hullNext[e])
            while true {
                let q = Int(_hullNext[n])
                if orient2d(p,
                            SIMD2<Double>(coords[2 * n], coords[2 * n + 1]),
                            SIMD2<Double>(coords[2 * q], coords[2 * q + 1])) >= 0 { break }
                
                t = addTriangle(i0: UInt(n), i1: UInt(i), i2: UInt(q),
                                a: Int32(_hullTri[i]), b: -1, c: Int32(_hullTri[n]))
                _hullTri[i] = legalize(a: t + 2)
                _hullNext[n] = UInt(n)
                hullSize -= 1
                n = q
            }
            
            if e == start {
                while true {
                    let q = Int(_hullPrev[e])
                    if orient2d(p,
                                SIMD2<Double>(coords[2 * q], coords[2 * q + 1]),
                                SIMD2<Double>(coords[2 * e], coords[2 * e + 1])) >= 0 { break }
                    
                    t = addTriangle(i0: UInt(q), i1: UInt(i), i2: UInt(e),
                                    a: -1, b: Int32(_hullTri[e]), c: Int32(_hullTri[q]))
                    _ = legalize(a: t + 2)
                    _hullTri[q] = UInt(t)
                    _hullNext[e] = UInt(e)
                    hullSize -= 1
                    e = q
                }
            }
            
            _hullPrev[i] = UInt(e)
            _hullStart = UInt(e)
            _hullNext[e] = UInt(i)
            _hullPrev[n] = UInt(i)
            _hullNext[i] = UInt(n)
            
            _hullHash[hashKey(p)] = i
            _hullHash[hashKey(SIMD2<Double>(coords[2 * e], coords[2 * e + 1]))] = e
        }
        
        hull = Array(repeating: 0, count: hullSize)
        var e = Int(_hullStart)
        for i in 0..<hullSize {
            hull[i] = UInt(e)
            e = Int(_hullNext[e])
        }
        
        if trianglesLen > 0 {
            triangles = Array(triangles[0..<trianglesLen])
            halfedges = Array(halfedges[0..<trianglesLen])
        } else {
            triangles = []
            halfedges = []
        }
    }
    
    @inline(__always)
    private func hashKey(_ coords: SIMD2<Double>) -> Int {
        let angle = pseudoAngle(d: coords - _center)
        return Int(floor(Double(angle) * Double(_hashSize))) % _hashSize
    }
    
    private func legalize(a: Int) -> UInt {
        var i: Int = 0
        var ar: Int = 0
        var a = a
        
        while true {
            let b = Int(halfedges[a])
            
            if b == -1 {
                if i == 0 { break }
                a = Int(EDGE_STACK[i - 1])
                i -= 1
                continue
            }
            
            let a0 = a - a % 3
            ar = a0 + (a + 2) % 3
            
            let b0 = b - b % 3
            let al = a0 + (a + 1) % 3
            let bl = b0 + (b + 2) % 3
            
            let p0 = Int(triangles[ar])
            let pr = Int(triangles[a])
            let pl = Int(triangles[al])
            let p1 = Int(triangles[bl])
            
            let illegal = inCircle(SIMD2<Double>(coords[2 * p0], coords[2 * p0 + 1]),
                                   SIMD2<Double>(coords[2 * pr], coords[2 * pr + 1]),
                                   SIMD2<Double>(coords[2 * pl], coords[2 * pl + 1]),
                                   SIMD2<Double>(coords[2 * p1], coords[2 * p1 + 1]))
            
            if illegal {
                triangles[a] = UInt(p1)
                triangles[b] = UInt(p0)
                
                let hbl = halfedges[bl]
                
                // Continuing legalize function...
                if hbl == -1 {
                    var e = Int(_hullStart)
                    repeat {
                        if _hullTri[e] == UInt(bl) {
                            _hullTri[e] = UInt(a)
                            break
                        }
                        e = Int(_hullPrev[e])
                    } while e != Int(_hullStart)
                }
                
                link(a: a, b: Int(hbl))
                link(a: b, b: Int(halfedges[ar]))
                link(a: ar, b: bl)
                
                let br = b0 + (b + 1) % 3
                
                if i < EDGE_STACK.count {
                    EDGE_STACK[i] = UInt(br)
                    i += 1
                }
            } else {
                if i == 0 { break }
                i -= 1
                a = Int(EDGE_STACK[i])
            }
        }
        
        return UInt(ar)
    }
    
    private func link(a: Int, b: Int) {
        halfedges[a] = b
        if b != -1 {
            halfedges[b] = a
        }
    }
    
    private func addTriangle(i0: UInt, i1: UInt, i2: UInt,
                             a: Int32, b: Int32, c: Int32) -> Int {
        let t = trianglesLen
        
        triangles[t] = i0
        triangles[t + 1] = i1
        triangles[t + 2] = i2
        
        link(a: t, b: Int(a))
        link(a: t + 1, b: Int(b))
        link(a: t + 2, b: Int(c))
        
        trianglesLen += 3
        return t
    }
    
    @inline(__always)
    private func pseudoAngle(d: SIMD2<Double>) -> Double {
        let abs_d = abs(d)
        let p = d.x / (abs_d.x + abs_d.y)
        return (d.y > 0 ? 3 - p : 1 + p) / 4
    }
    
    @inline(__always)
    func squaredDistance(_ a: SIMD2<Double>, _ b: SIMD2<Double>) -> Double {
        let diff = a - b
        return dot(diff, diff)
    }
    
    @inline(__always)
    private func orient2d(_ a: SIMD2<Double>, _ b: SIMD2<Double>, _ c: SIMD2<Double>) -> Double {
        let ab = b - a
        let bc = c - b
        return ab.y * bc.x - ab.x * bc.y
    }
    
    @inline(__always)
    private func inCircle(_ a: SIMD2<Double>, _ b: SIMD2<Double>, _ c: SIMD2<Double>, _ p: SIMD2<Double>) -> Bool {
        // Compute vectors from p to each point
        let ap = a - p
        let bp = b - p
        let cp = c - p
        
        // Compute squared distances
        let ap_sq = dot(ap, ap)
        let bp_sq = dot(bp, bp)
        let cp_sq = dot(cp, cp)
        
        // Original determinant calculation which is numerically stable
        return (ap.x * (bp.y * cp_sq - bp_sq * cp.y) -
                ap.y * (bp.x * cp_sq - bp_sq * cp.x) +
                ap_sq * (bp.x * cp.y - bp.y * cp.x)) < 0
    }
    
    @inline(__always)
    private func circumradius(_ a: SIMD2<Double>, _ b: SIMD2<Double>, _ c: SIMD2<Double>) -> Double {
        let ab = b - a
        let ac = c - a
        
        let bl = dot(ab, ab)
        let cl = dot(ac, ac)
        let d = 0.5 / (ab.x * ac.y - ab.y * ac.x)
        
        let center = SIMD2<Double>(ac.y * bl - ab.y * cl,
                                  ab.x * cl - ac.x * bl) * d
        
        return dot(center, center)
    }
    
    @inline(__always)
    private func circumcenter(_ a: SIMD2<Double>, _ b: SIMD2<Double>, _ c: SIMD2<Double>) -> Point {
        let ab = b - a
        let ac = c - a
        
        let bl = dot(ab, ab)
        let cl = dot(ac, ac)
        let d = 0.5 / (ab.x * ac.y - ab.y * ac.x)
        
        let coords = a + SIMD2<Double>(ac.y * bl - ab.y * cl,
                                      ab.x * cl - ac.x * bl) * d
        
        return Point(coords)
    }
    
    
    private func quicksort(ids: inout [UInt], dists: [Double],
                           left: Int, right: Int) {
        if right <= left {
            return
        }
        
        if right - left <= 20 {
            // Insertion sort for small subarrays
            for i in (left + 1)...right {
                let temp = ids[i]
                let tempDist = dists[Int(temp)]
                var j = i - 1
                while j >= left && dists[Int(ids[j])] > tempDist {
                    ids[j + 1] = ids[j]
                    j -= 1
                }
                ids[j + 1] = temp
            }
            return
        }
        
        // Regular quicksort for larger arrays
        let median = (left + right) >> 1
        var i = left + 1
        var j = right
        
        swap(arr: &ids, i: median, j: i)
        
        if dists[Int(ids[left])] > dists[Int(ids[right])] {
            swap(arr: &ids, i: left, j: right)
        }
        if dists[Int(ids[i])] > dists[Int(ids[right])] {
            swap(arr: &ids, i: i, j: right)
        }
        if dists[Int(ids[left])] > dists[Int(ids[i])] {
            swap(arr: &ids, i: left, j: i)
        }
        
        let temp = ids[i]
        let tempDist = dists[Int(temp)]
        
        while true {
            repeat { i += 1 } while i <= right && dists[Int(ids[i])] < tempDist
            repeat { j -= 1 } while j >= left && dists[Int(ids[j])] > tempDist
            if j < i { break }
            swap(arr: &ids, i: i, j: j)
        }
        
        ids[left + 1] = ids[j]
        ids[j] = temp
        
        // Ensure the recursive calls have valid ranges
        if j > left {
            quicksort(ids: &ids, dists: dists, left: left, right: j - 1)
        }
        if j + 1 < right {
            quicksort(ids: &ids, dists: dists, left: j + 1, right: right)
        }
    }
    
    @inline(__always)
    private func swap(arr: inout [UInt], i: Int, j: Int) {
        let tmp = arr[i]
        arr[i] = arr[j]
        arr[j] = tmp
    }
}
