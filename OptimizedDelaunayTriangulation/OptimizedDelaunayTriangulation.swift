import Foundation
import Darwin

struct Point {
    let x: Double
    let y: Double
}

class Delaunator {
    private let EPSILON: Double = Darwin.pow(2.0, -52)
    private var EDGE_STACK = [UInt32](repeating: 0, count: 512)
    
    private var coords: [Double]
    public private(set) var triangles: [UInt32]
    public private(set) var halfedges: [Int32]
    private var trianglesLen: Int = 0
    public private(set) var hull: [UInt32] = []
    
    // Temporary arrays
    private var _hashSize: Int
    private var _hullPrev: [UInt32]
    private var _hullNext: [UInt32]
    private var _hullTri: [UInt32]
    private var _hullHash: [Int32]
    private var _ids: [UInt32]
    private var _dists: [Double]
    private var _cx: Double = 0
    private var _cy: Double = 0
    private var _hullStart: UInt32 = 0
    
    static func from(points: [Point]) -> Delaunator {
        let n = points.count
        var coords = [Double](repeating: 0, count: n * 2)
        
        for i in 0..<n {
            coords[2 * i] = points[i].x
            coords[2 * i + 1] = points[i].y
        }
        
        return Delaunator(coords: coords)
    }
    
    init(coords: [Double]) {
        let n = coords.count >> 1
        self.coords = coords
        
        let maxTriangles = max(2 * n - 5, 0)
        self.triangles = [UInt32](repeating: 0, count: maxTriangles * 3)
        self.halfedges = [Int32](repeating: 0, count: maxTriangles * 3)
        
        self._hashSize = Int(ceil(sqrt(Double(n))))
        self._hullPrev = [UInt32](repeating: 0, count: n)
        self._hullNext = [UInt32](repeating: 0, count: n)
        self._hullTri = [UInt32](repeating: 0, count: n)
        self._hullHash = [Int32](repeating: -1, count: self._hashSize)
        self._ids = [UInt32](repeating: 0, count: n)
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
            _ids[i] = UInt32(i)
        }
        
        let cx = (minX + maxX) / 2
        let cy = (minY + maxY) / 2
        
        var i0: Int = 0
        var minDist = Double.infinity
        
        for i in 0..<n {
            let d = dist(ax: cx, ay: cy, bx: coords[2 * i], by: coords[2 * i + 1])
            if d < minDist {
                i0 = i
                minDist = d
            }
        }
        
        let i0x = coords[2 * i0]
        let i0y = coords[2 * i0 + 1]
        
        var i1: Int = 0
        minDist = Double.infinity
        
        for i in 0..<n {
            if i == i0 { continue }
            let d = dist(ax: i0x, ay: i0y, bx: coords[2 * i], by: coords[2 * i + 1])
            if d < minDist && d > 0 {
                i1 = i
                minDist = d
            }
        }
        
        var i1x = coords[2 * i1]
        var i1y = coords[2 * i1 + 1]
        
        var minRadius = Double.infinity
        var i2: Int = 0
        
        for i in 0..<n {
            if i == i0 || i == i1 { continue }
            let r = circumradius(ax: i0x, ay: i0y, bx: i1x, by: i1y,
                                 cx: coords[2 * i], cy: coords[2 * i + 1])
            if r < minRadius {
                i2 = i
                minRadius = r
            }
        }
        
        var i2x = coords[2 * i2]
        var i2y = coords[2 * i2 + 1]
        
        if minRadius == Double.infinity {
            // Handle collinear points
            for i in 0..<n {
                _dists[i] = coords[2 * i] - coords[0] != 0 ?
                coords[2 * i] - coords[0] : coords[2 * i + 1] - coords[1]
            }
            quicksort(ids: &_ids, dists: _dists, left: 0, right: n - 1)
            var hull = [UInt32](repeating: 0, count: n)
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
        
        if orient2d(ax: i0x, ay: i0y, bx: i1x, by: i1y, cx: i2x, cy: i2y) < 0 {
            let i = i1
            let x = i1x
            let y = i1y
            i1 = i2
            i1x = i2x
            i1y = i2y
            i2 = i
            i2x = x
            i2y = y
        }
        
        let center = circumcenter(ax: i0x, ay: i0y, bx: i1x, by: i1y, cx: i2x, cy: i2y)
        self._cx = center.x
        self._cy = center.y
        
        for i in 0..<n {
            _dists[i] = dist(ax: coords[2 * i], ay: coords[2 * i + 1],
                             bx: center.x, by: center.y)
        }
        
        quicksort(ids: &_ids, dists: _dists, left: 0, right: n - 1)
        
        self._hullStart = UInt32(i0)
        var hullSize = 3
        
        _hullNext[i0] = UInt32(i1)
        _hullPrev[i2] = UInt32(i1)
        _hullNext[i1] = UInt32(i2)
        _hullPrev[i0] = UInt32(i2)
        _hullNext[i2] = UInt32(i0)
        _hullPrev[i1] = UInt32(i0)
        
        _hullTri[i0] = 0
        _hullTri[i1] = 1
        _hullTri[i2] = 2
        
        _hullHash = [Int32](repeating: -1, count: _hashSize)
        _hullHash[hashKey(x: i0x, y: i0y)] = Int32(i0)
        _hullHash[hashKey(x: i1x, y: i1y)] = Int32(i1)
        _hullHash[hashKey(x: i2x, y: i2y)] = Int32(i2)
        
        trianglesLen = 0
        _ = addTriangle(i0: UInt32(i0), i1: UInt32(i1), i2: UInt32(i2),
                        a: -1, b: -1, c: -1)
        
        var xp: Double = 0
        var yp: Double = 0
        
        for k in 0..<_ids.count {
            let i = Int(_ids[k])
            let x = coords[2 * i]
            let y = coords[2 * i + 1]
            
            if k > 0 && abs(x - xp) <= EPSILON && abs(y - yp) <= EPSILON { continue }
            xp = x
            yp = y
            
            if i == i0 || i == i1 || i == i2 { continue }
            
            var start: Int = 0
            let key = hashKey(x: x, y: y)
            
            for j in 0..<_hashSize {
                start = Int(_hullHash[(key + j) % _hashSize])
                if start != -1 && start != Int(_hullNext[start]) { break }
            }
            
            start = Int(_hullPrev[start])
            var e = start
            
            while true {
                let q = Int(_hullNext[e])
                if orient2d(ax: x, ay: y,
                            bx: coords[2 * e], by: coords[2 * e + 1],
                            cx: coords[2 * q], cy: coords[2 * q + 1]) < 0 { break }
                e = q
                if e == start {
                    e = -1
                    break
                }
            }
            
            if e == -1 { continue }
            
            var t = addTriangle(i0: UInt32(e), i1: UInt32(i),
                                i2: _hullNext[e],
                                a: -1, b: -1, c: Int32(_hullTri[e]))
            
            _hullTri[i] = legalize(a: t + 2)
            _hullTri[e] = UInt32(t)
            hullSize += 1
            
            var n = Int(_hullNext[e])
            while true {
                let q = Int(_hullNext[n])
                if orient2d(ax: x, ay: y,
                            bx: coords[2 * n], by: coords[2 * n + 1],
                            cx: coords[2 * q], cy: coords[2 * q + 1]) >= 0 { break }
                
                t = addTriangle(i0: UInt32(n), i1: UInt32(i), i2: UInt32(q),
                                a: Int32(_hullTri[i]), b: -1, c: Int32(_hullTri[n]))
                _hullTri[i] = legalize(a: t + 2)
                _hullNext[n] = UInt32(n)
                hullSize -= 1
                n = q
            }
            
            if e == start {
                while true {
                    let q = Int(_hullPrev[e])
                    if orient2d(ax: x, ay: y,
                                bx: coords[2 * q], by: coords[2 * q + 1],
                                cx: coords[2 * e], cy: coords[2 * e + 1]) >= 0 { break }
                    
                    t = addTriangle(i0: UInt32(q), i1: UInt32(i), i2: UInt32(e),
                                    a: -1, b: Int32(_hullTri[e]), c: Int32(_hullTri[q]))
                    _ = legalize(a: t + 2)
                    _hullTri[q] = UInt32(t)
                    _hullNext[e] = UInt32(e)
                    hullSize -= 1
                    e = q
                }
            }
            
            _hullPrev[i] = UInt32(e)
            _hullStart = UInt32(e)
            _hullNext[e] = UInt32(i)
            _hullPrev[n] = UInt32(i)
            _hullNext[i] = UInt32(n)
            
            _hullHash[hashKey(x: x, y: y)] = Int32(i)
            _hullHash[hashKey(x: coords[2 * e], y: coords[2 * e + 1])] = Int32(e)
        }
        
        hull = Array(repeating: 0, count: hullSize)
        var e = Int(_hullStart)
        for i in 0..<hullSize {
            hull[i] = UInt32(e)
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
    
    // Helper functions remain the same...
    private func hashKey(x: Double, y: Double) -> Int {
        let angle = pseudoAngle(dx: x - _cx, dy: y - _cy)
        return Int(floor(Double(angle) * Double(_hashSize))) % _hashSize
    }
    
    private func legalize(a: Int) -> UInt32 {
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
            
            let illegal = inCircle(ax: coords[2 * p0], ay: coords[2 * p0 + 1],
                                   bx: coords[2 * pr], by: coords[2 * pr + 1],
                                   cx: coords[2 * pl], cy: coords[2 * pl + 1],
                                   px: coords[2 * p1], py: coords[2 * p1 + 1])
            
            if illegal {
                triangles[a] = UInt32(p1)
                triangles[b] = UInt32(p0)
                
                let hbl = halfedges[bl]
                
                // Continuing legalize function...
                if hbl == -1 {
                    var e = Int(_hullStart)
                    repeat {
                        if _hullTri[e] == UInt32(bl) {
                            _hullTri[e] = UInt32(a)
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
                    EDGE_STACK[i] = UInt32(br)
                    i += 1
                }
            } else {
                if i == 0 { break }
                i -= 1
                a = Int(EDGE_STACK[i])
            }
        }
        
        return UInt32(ar)
    }
    
    private func link(a: Int, b: Int) {
        halfedges[a] = Int32(b)
        if b != -1 {
            halfedges[b] = Int32(a)
        }
    }
    
    private func addTriangle(i0: UInt32, i1: UInt32, i2: UInt32,
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
    
    private func pseudoAngle(dx: Double, dy: Double) -> Double {
        let p = dx / (abs(dx) + abs(dy))
        return (dy > 0 ? 3 - p : 1 + p) / 4
    }
    
    private func dist(ax: Double, ay: Double, bx: Double, by: Double) -> Double {
        let dx = ax - bx
        let dy = ay - by
        return dx * dx + dy * dy
    }
    
    private func orient2d(ax: Double, ay: Double,
                          bx: Double, by: Double,
                          cx: Double, cy: Double) -> Double {
        return (by - ay) * (cx - bx) - (bx - ax) * (cy - by)
    }
    
    private func inCircle(ax: Double, ay: Double,
                          bx: Double, by: Double,
                          cx: Double, cy: Double,
                          px: Double, py: Double) -> Bool {
        let dx = ax - px
        let dy = ay - py
        let ex = bx - px
        let ey = by - py
        let fx = cx - px
        let fy = cy - py
        
        let ap = dx * dx + dy * dy
        let bp = ex * ex + ey * ey
        let cp = fx * fx + fy * fy
        
        return dx * (ey * cp - bp * fy) -
        dy * (ex * cp - bp * fx) +
        ap * (ex * fy - ey * fx) < 0
    }
    
    private func circumradius(ax: Double, ay: Double,
                              bx: Double, by: Double,
                              cx: Double, cy: Double) -> Double {
        let dx = bx - ax
        let dy = by - ay
        let ex = cx - ax
        let ey = cy - ay
        
        let bl = dx * dx + dy * dy
        let cl = ex * ex + ey * ey
        let d = 0.5 / (dx * ey - dy * ex)
        
        let x = (ey * bl - dy * cl) * d
        let y = (dx * cl - ex * bl) * d
        
        return x * x + y * y
    }
    
    private func circumcenter(ax: Double, ay: Double,
                              bx: Double, by: Double,
                              cx: Double, cy: Double) -> Point {
        let dx = bx - ax
        let dy = by - ay
        let ex = cx - ax
        let ey = cy - ay
        
        let bl = dx * dx + dy * dy
        let cl = ex * ex + ey * ey
        let d = 0.5 / (dx * ey - dy * ex)
        
        let x = ax + (ey * bl - dy * cl) * d
        let y = ay + (dx * cl - ex * bl) * d
        
        return Point(x: x, y: y)
    }
    
    private func quicksort(ids: inout [UInt32], dists: [Double],
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
    
    private func swap(arr: inout [UInt32], i: Int, j: Int) {
        let tmp = arr[i]
        arr[i] = arr[j]
        arr[j] = tmp
    }
}
