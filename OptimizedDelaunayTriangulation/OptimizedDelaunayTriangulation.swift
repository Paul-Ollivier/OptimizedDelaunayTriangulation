import Foundation
import Darwin
import simd
import Dispatch

typealias DPoint = SIMD2<Double>

struct Point {
    let coords: DPoint
    
    init(_ coords: DPoint) {
        self.coords = coords
    }
}

private struct DTriangle {
    var i0, i1, i2: Int
    var p0, p1, p2: DPoint
}

class DelaunatorCache {
    var triangles: [UInt]
    var halfedges: [Int]
    var hull: [UInt]
    var hullPrev: [UInt]
    var hullNext: [UInt]
    var hullTri: [UInt]
    var hullHash: [Int]
    var ids: [UInt]
    var dists: [Double]
    var coords: [Double]
    var hullStart: UInt = 0  // Added this property
    
    init(maxPoints: Int) {
        let maxTriangles = max(2 * maxPoints - 5, 0)
        let hashSize = Int(ceil(sqrt(Double(maxPoints))))
        
        self.triangles = [UInt](repeating: 0, count: maxTriangles * 3)
        self.halfedges = [Int](repeating: 0, count: maxTriangles * 3)
        self.hull = [UInt](repeating: 0, count: maxPoints)
        self.hullPrev = [UInt](repeating: 0, count: maxPoints)
        self.hullNext = [UInt](repeating: 0, count: maxPoints)
        self.hullTri = [UInt](repeating: 0, count: maxPoints)
        self.hullHash = [Int](repeating: -1, count: hashSize)
        self.ids = [UInt](repeating: 0, count: maxPoints)
        self.dists = [Double](repeating: 0, count: maxPoints)
        self.coords = [Double](repeating: 0, count: maxPoints * 2)
    }
    
    func reset() {
        hullHash.withUnsafeMutableBufferPointer { ptr in
            ptr.baseAddress?.initialize(repeating: -1, count: ptr.count)
        }
        hullStart = 0
    }
}

class CachedDelaunator {
    private let EPSILON: Double = Darwin.pow(2.0, -52)
    private var EDGE_STACK = [UInt](repeating: 0, count: 512)
    
    private let cache: DelaunatorCache
    private var trianglesLen: Int = 0
    private var center: DPoint = DPoint(0, 0)
    private var currentPointCount: Int = 0
    
    private let processingQueue: DispatchQueue
    
    public var triangles: [UInt] { Array(cache.triangles[0..<trianglesLen]) }
    public var halfedges: [Int] { Array(cache.halfedges[0..<trianglesLen]) }
    public var hull: [UInt] { Array(cache.hull[0..<currentPointCount]) }
    
    init(maxPoints: Int) {
        self.cache = DelaunatorCache(maxPoints: maxPoints)
        self.processingQueue = DispatchQueue(label: "com.delaunator.processing",
                                             qos: .userInitiated,
                                             attributes: .concurrent)
    }
    
    func update(from points: UnsafeBufferPointer<Point>) {
        let n = points.count
        currentPointCount = n
        
        cache.coords.withUnsafeMutableBufferPointer { coordsPtr in
            for i in 0..<n {
                let point = points[i]
                let baseIndex = 2 * i
                coordsPtr[baseIndex] = point.coords.x
                coordsPtr[baseIndex + 1] = point.coords.y
            }
        }
        
        triangulate()
    }
    
    func update(from points: [Point]) {
        points.withUnsafeBufferPointer { pointsPtr in
            update(from: pointsPtr)
        }
    }
    
    private func triangulate() {
        let n = currentPointCount
        cache.reset()
        trianglesLen = 0
        
        // Initialize parallel processing groups
        let boundsGroup = DispatchGroup()
        var minX = Double.infinity
        var minY = Double.infinity
        var maxX = -Double.infinity
        var maxY = -Double.infinity
        let boundsLock = NSLock()
        
        // Compute bounds in parallel
        let chunkSize = max(1, n / ProcessInfo.processInfo.activeProcessorCount)
        for start in stride(from: 0, to: n, by: chunkSize) {
            let end = min(start + chunkSize, n)
            processingQueue.async(group: boundsGroup) {
                var localMinX = Double.infinity
                var localMinY = Double.infinity
                var localMaxX = -Double.infinity
                var localMaxY = -Double.infinity
                
                for i in start..<end {
                    let x = self.cache.coords[2 * i]
                    let y = self.cache.coords[2 * i + 1]
                    localMinX = min(localMinX, x)
                    localMinY = min(localMinY, y)
                    localMaxX = max(localMaxX, x)
                    localMaxY = max(localMaxY, y)
                }
                
                boundsLock.lock()
                minX = min(minX, localMinX)
                minY = min(minY, localMinY)
                maxX = max(maxX, localMaxX)
                maxY = max(maxY, localMaxY)
                boundsLock.unlock()
            }
        }
        boundsGroup.wait()
        
        let c = DPoint((minX + maxX) / 2, (minY + maxY) / 2)
        
        // Find initial point
        var i0: Int = 0
        var minDist = Double.infinity
        let distanceGroup = DispatchGroup()
        let distanceLock = NSLock()
        
        for start in stride(from: 0, to: n, by: chunkSize) {
            let end = min(start + chunkSize, n)
            processingQueue.async(group: distanceGroup) {
                var localMinDist = Double.infinity
                var localI0 = 0
                
                for i in start..<end {
                    let d = self.squaredDistance(c, DPoint(self.cache.coords[2 * i],
                                                           self.cache.coords[2 * i + 1]))
                    if d < localMinDist {
                        localI0 = i
                        localMinDist = d
                    }
                }
                
                distanceLock.lock()
                if localMinDist < minDist {
                    minDist = localMinDist
                    i0 = localI0
                }
                distanceLock.unlock()
            }
        }
        distanceGroup.wait()
        
        let i0p = DPoint(cache.coords[2 * i0], cache.coords[2 * i0 + 1])
        
        // Initialize ids array
        DispatchQueue.concurrentPerform(iterations: n) { i in
            cache.ids[i] = UInt(i)
        }
        
        // Find second point
        var i1: Int = 0
        minDist = Double.infinity
        let secondPointGroup = DispatchGroup()
        
        for start in stride(from: 0, to: n, by: chunkSize) {
            let end = min(start + chunkSize, n)
            processingQueue.async(group: secondPointGroup) {
                var localMinDist = Double.infinity
                var localI1 = 0
                
                for i in start..<end {
                    if i == i0 { continue }
                    let d = self.squaredDistance(i0p, DPoint(self.cache.coords[2 * i],
                                                             self.cache.coords[2 * i + 1]))
                    if d < localMinDist && d > 0 {
                        localI1 = i
                        localMinDist = d
                    }
                }
                
                distanceLock.lock()
                if localMinDist < minDist {
                    minDist = localMinDist
                    i1 = localI1
                }
                distanceLock.unlock()
            }
        }
        secondPointGroup.wait()
        
        var triangle = DTriangle(
            i0: i0,
            i1: i1,
            i2: 0,
            p0: i0p,
            p1: DPoint(cache.coords[2 * i1], cache.coords[2 * i1 + 1]),
            p2: DPoint(0, 0)
        )
        
        // Find third point
        var minRadius = Double.infinity
        var i2: Int = 0
        let thirdPointGroup = DispatchGroup()
        let radiusLock = NSLock()
        
        for start in stride(from: 0, to: n, by: chunkSize) {
            let end = min(start + chunkSize, n)
            processingQueue.async(group: thirdPointGroup) {
                var localMinRadius = Double.infinity
                var localI2 = 0
                
                for i in start..<end {
                    if i == triangle.i0 || i == triangle.i1 { continue }
                    let p = DPoint(self.cache.coords[2 * i], self.cache.coords[2 * i + 1])
                    let r = self.circumradius(triangle.p0, triangle.p1, p)
                    if r < localMinRadius {
                        localI2 = i
                        localMinRadius = r
                    }
                }
                
                radiusLock.lock()
                if localMinRadius < minRadius {
                    minRadius = localMinRadius
                    i2 = localI2
                }
                radiusLock.unlock()
            }
        }
        thirdPointGroup.wait()
        
        if minRadius == Double.infinity {
            handleCollinearPoints(n)
            return
        }
        
        triangle.i2 = i2
        triangle.p2 = DPoint(cache.coords[2 * i2], cache.coords[2 * i2 + 1])
        
        if orient2d(triangle.p0, triangle.p1, triangle.p2) < 0 {
            let tempI = triangle.i1
            let tempP = triangle.p1
            triangle.i1 = triangle.i2
            triangle.p1 = triangle.p2
            triangle.i2 = tempI
            triangle.p2 = tempP
        }
        
        center = circumcenter(triangle.p0, triangle.p1, triangle.p2).coords
        
        // Parallel distance computation
        DispatchQueue.concurrentPerform(iterations: n) { i in
            cache.dists[i] = squaredDistance(DPoint(cache.coords[2 * i],
                                                    cache.coords[2 * i + 1]), center)
        }
        
        quicksort(ids: &cache.ids, dists: cache.dists, left: 0, right: n - 1)
        
        cache.hullStart = UInt(triangle.i0)
        var hullSize = 3
        
        cache.hullNext[triangle.i0] = UInt(triangle.i1)
        cache.hullPrev[triangle.i2] = UInt(triangle.i1)
        cache.hullNext[triangle.i1] = UInt(triangle.i2)
        cache.hullPrev[triangle.i0] = UInt(triangle.i2)
        cache.hullNext[triangle.i2] = UInt(triangle.i0)
        cache.hullPrev[triangle.i1] = UInt(triangle.i0)
        
        cache.hullTri[triangle.i0] = 0
        cache.hullTri[triangle.i1] = 1
        cache.hullTri[triangle.i2] = 2
        
        cache.hullHash[hashKey(triangle.p0)] = triangle.i0
        cache.hullHash[hashKey(triangle.p1)] = triangle.i1
        cache.hullHash[hashKey(triangle.p2)] = triangle.i2
        
        trianglesLen = 0
        _ = addTriangle(i0: UInt(triangle.i0), i1: UInt(triangle.i1), i2: UInt(triangle.i2),
                        a: -1, b: -1, c: -1)
        
        var xp = Double.infinity
        var yp = Double.infinity
        
        for k in 0..<n {
            let i = Int(cache.ids[k])
            let x = cache.coords[2 * i]
            let y = cache.coords[2 * i + 1]
            
            if k > 0 && abs(x - xp) <= EPSILON && abs(y - yp) <= EPSILON { continue }
            xp = x
            yp = y
            
            if i == triangle.i0 || i == triangle.i1 || i == triangle.i2 { continue }
            
            var start = 0
            let key = hashKey(DPoint(x, y))
            
            for j in 0..<cache.hullHash.count {
                start = Int(cache.hullHash[(key + j) % cache.hullHash.count])
                if start != -1 && start != Int(cache.hullNext[start]) { break }
            }
            
            start = Int(cache.hullPrev[start])
            var e = start
            
            while true {
                let q = Int(cache.hullNext[e])
                if orient2d(DPoint(x, y),
                            DPoint(cache.coords[2 * e], cache.coords[2 * e + 1]),
                            DPoint(cache.coords[2 * q], cache.coords[2 * q + 1])) < 0 { break }
                e = q
                if e == start {
                    e = -1
                    break
                }
            }
            
            if e == -1 { continue }
            
            var t = addTriangle(i0: UInt(e), i1: UInt(i), i2: cache.hullNext[e],
                                a: -1, b: -1, c: Int32(cache.hullTri[e]))
            
            cache.hullTri[i] = legalize(a: t + 2)
            cache.hullTri[e] = UInt(t)
            hullSize += 1
            
            var n = Int(cache.hullNext[e])
            while true {
                let q = Int(cache.hullNext[n])
                if orient2d(DPoint(x, y),
                            DPoint(cache.coords[2 * n], cache.coords[2 * n + 1]),
                            DPoint(cache.coords[2 * q], cache.coords[2 * q + 1])) >= 0 { break }
                
                t = addTriangle(i0: UInt(n), i1: UInt(i), i2: UInt(q),
                                a: Int32(cache.hullTri[i]), b: -1, c: Int32(cache.hullTri[n]))
                
                cache.hullTri[i] = legalize(a: t + 2)
                cache.hullNext[n] = UInt(n)
                hullSize -= 1
                n = q
            }
            
            if e == start {
                while true {
                    let q = Int(cache.hullPrev[e])
                    if orient2d(DPoint(x, y),
                                DPoint(cache.coords[2 * q], cache.coords[2 * q + 1]),
                                DPoint(cache.coords[2 * e], cache.coords[2 * e + 1])) >= 0 { break }
                    
                    t = addTriangle(i0: UInt(q), i1: UInt(i), i2: UInt(e),
                                    a: -1, b: Int32(cache.hullTri[e]), c: Int32(cache.hullTri[q]))
                    
                    _ = legalize(a: t + 2)
                    cache.hullTri[q] = UInt(t)
                    cache.hullNext[e] = UInt(e)
                    hullSize -= 1
                    e = q
                }
            }
            
            cache.hullPrev[i] = UInt(e)
            cache.hullStart = UInt(e)
            cache.hullNext[e] = UInt(i)
            cache.hullPrev[n] = UInt(i)
            cache.hullNext[i] = UInt(n)
            
            cache.hullHash[hashKey(DPoint(x, y))] = i
            cache.hullHash[hashKey(DPoint(cache.coords[2 * e], cache.coords[2 * e + 1]))] = e
        }
        
        // Update hull array
        cache.hull.withUnsafeMutableBufferPointer { hullPtr in
            var e = Int(cache.hullStart)
            for i in 0..<hullSize {
                hullPtr[i] = UInt(e)
                e = Int(cache.hullNext[e])
            }
        }
    }
    
    private func handleCollinearPoints(_ n: Int) {
        DispatchQueue.concurrentPerform(iterations: n) { i in
            cache.dists[i] = cache.coords[2 * i] - cache.coords[0] != 0 ?
            cache.coords[2 * i] - cache.coords[0] :
            cache.coords[2 * i + 1] - cache.coords[1]
        }
        
        quicksort(ids: &cache.ids, dists: cache.dists, left: 0, right: n - 1)
        
        var j = 0
        var d0 = -Double.infinity
        
        for i in 0..<n {
            let id = cache.ids[i]
            if cache.dists[Int(id)] > d0 {
                cache.hull[j] = id
                j += 1
                d0 = cache.dists[Int(id)]
            }
        }
        
        trianglesLen = 0
    }
    
    private func quicksort(ids: inout [UInt], dists: [Double], left: Int, right: Int) {
        if left >= right { return }  // Added this check first
        
        if right - left <= 20 {
            // Use insertion sort for small arrays
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
            repeat { i += 1 } while i <= right && dists[Int(ids[i])] < tempDist  // Added bounds check
            repeat { j -= 1 } while j >= left && dists[Int(ids[j])] > tempDist   // Added bounds check
            if j < i { break }
            swap(arr: &ids, i: i, j: j)
        }
        
        ids[left + 1] = ids[j]
        ids[j] = temp
        
        // Recursively sort the smaller partition first
        if right - j < j - left {
            quicksort(ids: &ids, dists: dists, left: j + 1, right: right)
            quicksort(ids: &ids, dists: dists, left: left, right: j - 1)
        } else {
            quicksort(ids: &ids, dists: dists, left: left, right: j - 1)
            quicksort(ids: &ids, dists: dists, left: j + 1, right: right)
        }
    }
    
    @inline(__always)
    private func swap(arr: inout [UInt], i: Int, j: Int) {
        let tmp = arr[i]
        arr[i] = arr[j]
        arr[j] = tmp
    }
    
    private func legalize(a: Int) -> UInt {
        var i = 0
        var ar = 0
        var a = a
        
        while true {
            let b = cache.halfedges[a]
            
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
            
            let p0 = Int(cache.triangles[ar])
            let pr = Int(cache.triangles[a])
            let pl = Int(cache.triangles[al])
            let p1 = Int(cache.triangles[bl])
            
            let illegal = inCircle(
                DPoint(cache.coords[2 * p0], cache.coords[2 * p0 + 1]),
                DPoint(cache.coords[2 * pr], cache.coords[2 * pr + 1]),
                DPoint(cache.coords[2 * pl], cache.coords[2 * pl + 1]),
                DPoint(cache.coords[2 * p1], cache.coords[2 * p1 + 1])
            )
            
            if illegal {
                cache.triangles[a] = UInt(p1)
                cache.triangles[b] = UInt(p0)
                
                let hbl = cache.halfedges[bl]
                
                if hbl == -1 {
                    var e = Int(cache.hullStart)  // Updated reference
                    repeat {
                        if cache.hullTri[e] == UInt(bl) {
                            cache.hullTri[e] = UInt(a)
                            break
                        }
                        e = Int(cache.hullPrev[e])
                    } while e != Int(cache.hullStart)  // Updated reference
                }
                
                link(a: a, b: Int(hbl))
                link(a: b, b: cache.halfedges[ar])
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
        cache.halfedges[a] = b
        if b != -1 {
            cache.halfedges[b] = a
        }
    }
    
    private func addTriangle(i0: UInt, i1: UInt, i2: UInt, a: Int32, b: Int32, c: Int32) -> Int {
        let t = trianglesLen
        
        cache.triangles[t] = i0
        cache.triangles[t + 1] = i1
        cache.triangles[t + 2] = i2
        
        link(a: t, b: Int(a))
        link(a: t + 1, b: Int(b))
        link(a: t + 2, b: Int(c))
        
        trianglesLen += 3
        return t
    }
    
    @inline(__always)
    private func hashKey(_ coords: DPoint) -> Int {
        let angle = pseudoAngle(d: coords - center)
        return Int(floor(Double(angle) * Double(cache.hullHash.count))) % cache.hullHash.count
    }
    
    @inline(__always)
    private func pseudoAngle(d: DPoint) -> Double {
        let abs_d = abs(d)
        let sum = abs_d.x + abs_d.y
        guard sum != 0 else { return 0 }
        let p = d.x / sum
        return (d.y > 0 ? 3 - p : 1 + p) / 4
    }
    
    @inline(__always)
    private func squaredDistance(_ a: DPoint, _ b: DPoint) -> Double {
        return simd_distance_squared(a, b)
    }
    
    @inline(__always)
    private func orient2d(_ a: DPoint, _ b: DPoint, _ c: DPoint) -> Double {
        let ab = b - a
        let bc = c - b
        return ab.y * bc.x - ab.x * bc.y
    }
    
    @inline(__always)
    private func inCircle(_ a: DPoint, _ b: DPoint, _ c: DPoint, _ p: DPoint) -> Bool {
        let ap = a - p
        let bp = b - p
        let cp = c - p
        
        let ap_sq = simd_dot(ap, ap)
        let bp_sq = simd_dot(bp, bp)
        let cp_sq = simd_dot(cp, cp)
        
        let m = SIMD3<Double>(
            ap_sq * simd_cross(bp, cp).z,
            bp_sq * simd_cross(cp, ap).z,
            cp_sq * simd_cross(ap, bp).z
        )
        return simd_reduce_add(m) < 0
    }
    
    @inline(__always)
    private func circumcenter(_ a: DPoint, _ b: DPoint, _ c: DPoint) -> Point {
        let ab = b - a
        let ac = c - a
        
        let bl = simd_dot(ab, ab)
        let cl = simd_dot(ac, ac)
        
        let d = 0.5 / simd_cross(ab, ac).z
        
        let coords = a + DPoint(
            ac.y * bl - ab.y * cl,
            ab.x * cl - ac.x * bl
        ) * d
        
        return Point(coords)
    }
    
    @inline(__always)
    private func circumradius(_ a: DPoint, _ b: DPoint, _ c: DPoint) -> Double {
        let ab = b - a
        let ac = c - a
        
        let bl = simd_dot(ab, ab)
        let cl = simd_dot(ac, ac)
        
        let d = 0.5 / simd_cross(ab, ac).z
        
        let center = DPoint(
            ac.y * bl - ab.y * cl,
            ab.x * cl - ac.x * bl
        ) * d
        
        return simd_dot(center, center)
    }
}
