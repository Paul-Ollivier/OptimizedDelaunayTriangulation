import Foundation
import Darwin
import simd

typealias DPoint = SIMD2<Double>

struct Point {
    public let coords: DPoint
    
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
    var hullStart: UInt = 0
    
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
    
    public var triangles: [UInt] { Array(cache.triangles[0..<trianglesLen]) }
    public var halfedges: [Int] { Array(cache.halfedges[0..<trianglesLen]) }
    public var hull: [UInt] { Array(cache.hull[0..<currentPointCount]) }
    
    init(maxPoints: Int) {
        self.cache = DelaunatorCache(maxPoints: maxPoints)
    }
    
    func update(from points: [Point]) {
        points.withUnsafeBufferPointer { pointsPtr in
            let n = pointsPtr.count
            currentPointCount = n
            
            cache.coords.withUnsafeMutableBufferPointer { coordsPtr in
                for i in 0..<n {
                    let point = pointsPtr[i]
                    let baseIndex = 2 * i
                    coordsPtr[baseIndex] = point.coords.x
                    coordsPtr[baseIndex + 1] = point.coords.y
                }
            }
            
            triangulate()
        }
    }
    
    private func triangulate() {
        let n = currentPointCount
        cache.reset()
        trianglesLen = 0
        
        // Find bounds
        var minX = Double.infinity
        var minY = Double.infinity
        var maxX = -Double.infinity
        var maxY = -Double.infinity
        
        for i in 0..<n {
            let x = cache.coords[2 * i]
            let y = cache.coords[2 * i + 1]
            minX = min(minX, x)
            minY = min(minY, y)
            maxX = max(maxX, x)
            maxY = max(maxY, y)
        }
        
        var center = DPoint((minX + maxX) / 2, (minY + maxY) / 2)
        
        // Find first vertex
        var firstIdx = 0
        var minDist = Double.infinity
        
        for i in 0..<n {
            let dist = squaredDistance(center, DPoint(cache.coords[2 * i], cache.coords[2 * i + 1]))
            if dist < minDist {
                firstIdx = i
                minDist = dist
            }
        }
        
        let first = DPoint(cache.coords[2 * firstIdx], cache.coords[2 * firstIdx + 1])
        
        // Initialize ids array
        for i in 0..<n {
            cache.ids[i] = UInt(i)
        }
        
        // Find second vertex
        var secondIdx = 0
        minDist = Double.infinity
        
        for i in 0..<n {
            if i == firstIdx { continue }
            let dist = squaredDistance(first, DPoint(cache.coords[2 * i], cache.coords[2 * i + 1]))
            if dist < minDist && dist > 0 {
                secondIdx = i
                minDist = dist
            }
        }
        
        var tri = DTriangle(
            i0: firstIdx,
            i1: secondIdx,
            i2: 0,
            p0: first,
            p1: DPoint(cache.coords[2 * secondIdx], cache.coords[2 * secondIdx + 1]),
            p2: DPoint(0, 0)
        )
        
        // Find third vertex
        var minRadius = Double.infinity
        var thirdIdx = 0
        
        for i in 0..<n {
            if i == tri.i0 || i == tri.i1 { continue }
            let point = DPoint(cache.coords[2 * i], cache.coords[2 * i + 1])
            let radius = circumradius(tri.p0, tri.p1, point)
            if radius < minRadius {
                thirdIdx = i
                minRadius = radius
            }
        }
        
        if minRadius == Double.infinity {
            handleCollinearPoints(n)
            return
        }
        
        tri.i2 = thirdIdx
        tri.p2 = DPoint(cache.coords[2 * thirdIdx], cache.coords[2 * thirdIdx + 1])
        
        if orient2d(tri.p0, tri.p1, tri.p2) < 0 {
            let temp = tri.i1
            let tempPoint = tri.p1
            tri.i1 = tri.i2
            tri.p1 = tri.p2
            tri.i2 = temp
            tri.p2 = tempPoint
        }
        
        center = circumcenter(tri.p0, tri.p1, tri.p2).coords
        
        // Calculate distances
        for i in 0..<n {
            cache.dists[i] = squaredDistance(
                DPoint(cache.coords[2 * i], cache.coords[2 * i + 1]),
                center
            )
        }
        
        quicksort(ids: &cache.ids, dists: cache.dists, left: 0, right: n - 1)
        
        // Hull initialization
        cache.hullStart = UInt(tri.i0)
        var hullSize = 3
        
        cache.hullNext[tri.i0] = UInt(tri.i1)
        cache.hullPrev[tri.i2] = UInt(tri.i1)
        cache.hullNext[tri.i1] = UInt(tri.i2)
        cache.hullPrev[tri.i0] = UInt(tri.i2)
        cache.hullNext[tri.i2] = UInt(tri.i0)
        cache.hullPrev[tri.i1] = UInt(tri.i0)
        
        cache.hullTri[tri.i0] = 0
        cache.hullTri[tri.i1] = 1
        cache.hullTri[tri.i2] = 2
        
        cache.hullHash[hashKey(tri.p0)] = tri.i0
        cache.hullHash[hashKey(tri.p1)] = tri.i1
        cache.hullHash[hashKey(tri.p2)] = tri.i2
        
        trianglesLen = 0
        _ = addTriangle(i0: UInt(tri.i0), i1: UInt(tri.i1), i2: UInt(tri.i2),
                        a: -1, b: -1, c: -1)
        
        var prevX = Double.infinity
        var prevY = Double.infinity
        
        for k in 0..<n {
            let i = Int(cache.ids[k])
            let x = cache.coords[2 * i]
            let y = cache.coords[2 * i + 1]
            
            if k > 0 && abs(x - prevX) <= EPSILON && abs(y - prevY) <= EPSILON { continue }
            prevX = x
            prevY = y
            
            if i == tri.i0 || i == tri.i1 || i == tri.i2 { continue }
            
            var start = 0
            let key = hashKey(DPoint(x, y))
            
            for j in 0..<cache.hullHash.count {
                start = Int(cache.hullHash[(key + j) % cache.hullHash.count])
                if start != -1 && start != Int(cache.hullNext[start]) { break }
            }
            
            start = Int(cache.hullPrev[start])
            var edge = start
            
            while true {
                let next = Int(cache.hullNext[edge])
                if orient2d(DPoint(x, y),
                            DPoint(cache.coords[2 * edge], cache.coords[2 * edge + 1]),
                            DPoint(cache.coords[2 * next], cache.coords[2 * next + 1])) < 0 { break }
                edge = next
                if edge == start {
                    edge = -1
                    break
                }
            }
            
            if edge == -1 { continue }
            
            var t = addTriangle(i0: UInt(edge), i1: UInt(i), i2: cache.hullNext[edge],
                                a: -1, b: -1, c: Int32(cache.hullTri[edge]))
            
            cache.hullTri[i] = legalize(edge: t + 2)
            cache.hullTri[edge] = UInt(t)
            hullSize += 1
            
            var next = Int(cache.hullNext[edge])
            while true {
                let nextNext = Int(cache.hullNext[next])
                if orient2d(DPoint(x, y),
                            DPoint(cache.coords[2 * next], cache.coords[2 * next + 1]),
                            DPoint(cache.coords[2 * nextNext], cache.coords[2 * nextNext + 1])) >= 0 { break }
                
                t = addTriangle(i0: UInt(next), i1: UInt(i), i2: UInt(nextNext),
                                a: Int32(cache.hullTri[i]), b: -1, c: Int32(cache.hullTri[next]))
                
                cache.hullTri[i] = legalize(edge: t + 2)
                cache.hullNext[next] = UInt(next)
                hullSize -= 1
                next = nextNext
            }
            
            if edge == start {
                while true {
                    let prev = Int(cache.hullPrev[edge])
                    if orient2d(DPoint(x, y),
                                DPoint(cache.coords[2 * prev], cache.coords[2 * prev + 1]),
                                DPoint(cache.coords[2 * edge], cache.coords[2 * edge + 1])) >= 0 { break }
                    
                    t = addTriangle(i0: UInt(prev), i1: UInt(i), i2: UInt(edge),
                                    a: -1, b: Int32(cache.hullTri[edge]), c: Int32(cache.hullTri[prev]))
                    
                    _ = legalize(edge: t + 2)
                    cache.hullTri[prev] = UInt(t)
                    cache.hullNext[edge] = UInt(edge)
                    hullSize -= 1
                    edge = prev
                }
            }
            
            cache.hullPrev[i] = UInt(edge)
            cache.hullStart = UInt(edge)
            cache.hullNext[edge] = UInt(i)
            cache.hullPrev[next] = UInt(i)
            cache.hullNext[i] = UInt(next)
            
            cache.hullHash[hashKey(DPoint(x, y))] = i
            cache.hullHash[hashKey(DPoint(cache.coords[2 * edge], cache.coords[2 * edge + 1]))] = edge
        }
        
        cache.hull.withUnsafeMutableBufferPointer { hullPtr in
            var e = Int(cache.hullStart)
            for i in 0..<hullSize {
                hullPtr[i] = UInt(e)
                e = Int(cache.hullNext[e])
            }
        }
    }
    
    private func handleCollinearPoints(_ n: Int) {
        for i in 0..<n {
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
        if left >= right { return }
        
        if right - left <= 20 {
            for curr in (left + 1)...right {
                let tmp = ids[curr]
                let tmpDist = dists[Int(tmp)]
                var prev = curr - 1
                while prev >= left && dists[Int(ids[prev])] > tmpDist {
                    ids[prev + 1] = ids[prev]
                    prev -= 1
                }
                ids[prev + 1] = tmp
            }
            return
        }
        
        let mid = (left + right) >> 1
        var i = left + 1
        var j = right
        swap(arr: &ids, index1: mid, index2: i)
        
        if dists[Int(ids[left])] > dists[Int(ids[right])] {
            swap(arr: &ids, index1: left, index2: right)
        }
        if dists[Int(ids[i])] > dists[Int(ids[right])] {
            swap(arr: &ids, index1: i, index2: right)
        }
        if dists[Int(ids[left])] > dists[Int(ids[i])] {
            swap(arr: &ids, index1: left, index2: i)
        }
        
        let tmp = ids[i]
        let tmpDist = dists[Int(tmp)]
        
        while true {
            repeat { i += 1 } while i <= right && dists[Int(ids[i])] < tmpDist
            repeat { j -= 1 } while j >= left && dists[Int(ids[j])] > tmpDist
            if j < i { break }
            swap(arr: &ids, index1: i, index2: j)
        }
        
        ids[left + 1] = ids[j]
        ids[j] = tmp
        
        if right - j < j - left {
            quicksort(ids: &ids, dists: dists, left: j + 1, right: right)
            quicksort(ids: &ids, dists: dists, left: left, right: j - 1)
        } else {
            quicksort(ids: &ids, dists: dists, left: left, right: j - 1)
            quicksort(ids: &ids, dists: dists, left: j + 1, right: right)
        }
    }
    
    private func legalize(edge: Int) -> UInt {
        var stackPos = 0
        var result = 0
        var curr = edge
        
        while true {
            let opposite = cache.halfedges[curr]
            
            if opposite == -1 {
                if stackPos == 0 { break }
                curr = Int(EDGE_STACK[stackPos - 1])
                stackPos -= 1
                continue
            }
            
            let triStart = curr - curr % 3
            result = triStart + (curr + 2) % 3
            
            let oppTriStart = opposite - opposite % 3
            let adjLeft = triStart + (curr + 1) % 3
            let oppLeft = oppTriStart + (opposite + 2) % 3
            
            let v0 = Int(cache.triangles[result])
            let v1 = Int(cache.triangles[curr])
            let v2 = Int(cache.triangles[adjLeft])
            let v3 = Int(cache.triangles[oppLeft])
            
            let illegal = inCircle(
                DPoint(cache.coords[2 * v0], cache.coords[2 * v0 + 1]),
                DPoint(cache.coords[2 * v1], cache.coords[2 * v1 + 1]),
                DPoint(cache.coords[2 * v2], cache.coords[2 * v2 + 1]),
                DPoint(cache.coords[2 * v3], cache.coords[2 * v3 + 1])
            )
            
            if illegal {
                cache.triangles[curr] = UInt(v3)
                cache.triangles[opposite] = UInt(v0)
                
                let oppHull = cache.halfedges[oppLeft]
                
                if oppHull == -1 {
                    var e = Int(cache.hullStart)
                    repeat {
                        if cache.hullTri[e] == UInt(oppLeft) {
                            cache.hullTri[e] = UInt(curr)
                            break
                        }
                        e = Int(cache.hullPrev[e])
                    } while e != Int(cache.hullStart)
                }
                
                link(a: curr, b: Int(oppHull))
                link(a: opposite, b: cache.halfedges[result])
                link(a: result, b: oppLeft)
                
                let nextEdge = oppTriStart + (opposite + 1) % 3
                
                if stackPos < EDGE_STACK.count {
                    EDGE_STACK[stackPos] = UInt(nextEdge)
                    stackPos += 1
                }
            } else {
                if stackPos == 0 { break }
                stackPos -= 1
                curr = Int(EDGE_STACK[stackPos])
            }
        }
        
        return UInt(result)
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
    
    private func hashKey(_ p: DPoint) -> Int {
        let angle = pseudoAngle(delta: p - center)
        return Int(floor(Double(angle) * Double(cache.hullHash.count))) % cache.hullHash.count
    }
    
    private func pseudoAngle(delta: DPoint) -> Double {
        let magnitude = abs(delta)
        let manhattanDist = magnitude.x + magnitude.y
        guard manhattanDist != 0 else { return 0 }
        let normalizedX = delta.x / manhattanDist
        return (delta.y > 0 ? 3 - normalizedX : 1 + normalizedX) / 4
    }

    @inline(__always)
    private func squaredDistance(_ point1: DPoint, _ point2: DPoint) -> Double {
        return simd_distance_squared(point1, point2)
    }

    @inline(__always)
    private func orient2d(_ vertex1: DPoint, _ vertex2: DPoint, _ vertex3: DPoint) -> Double {
        let edge1 = vertex2 - vertex1  // First edge vector
        let edge2 = vertex3 - vertex2  // Second edge vector
        return edge1.y * edge2.x - edge1.x * edge2.y  // Cross product z component
    }

    @inline(__always)
    private func inCircle(_ vertex1: DPoint, _ vertex2: DPoint, _ vertex3: DPoint, _ test: DPoint) -> Bool {
        let vec1 = vertex1 - test
        let vec2 = vertex2 - test
        let vec3 = vertex3 - test
        
        let sqDist1 = simd_dot(vec1, vec1)
        let sqDist2 = simd_dot(vec2, vec2)
        let sqDist3 = simd_dot(vec3, vec3)
        
        let determinants = SIMD3<Double>(
            sqDist1 * simd_cross(vec2, vec3).z,
            sqDist2 * simd_cross(vec3, vec1).z,
            sqDist3 * simd_cross(vec1, vec2).z
        )
        return simd_reduce_add(determinants) < 0
    }

    @inline(__always)
    private func circumcenter(_ vertex1: DPoint, _ vertex2: DPoint, _ vertex3: DPoint) -> Point {
        let edge1 = vertex2 - vertex1
        let edge2 = vertex3 - vertex1
        
        let sqLen1 = simd_dot(edge1, edge1)
        let sqLen2 = simd_dot(edge2, edge2)
        
        let perpDivisor = 0.5 / simd_cross(edge1, edge2).z
        
        let center = vertex1 + DPoint(
            edge2.y * sqLen1 - edge1.y * sqLen2,
            edge1.x * sqLen2 - edge2.x * sqLen1
        ) * perpDivisor
        
        return Point(center)
    }

    @inline(__always)
    private func circumradius(_ vertex1: DPoint, _ vertex2: DPoint, _ vertex3: DPoint) -> Double {
        let edge1 = vertex2 - vertex1
        let edge2 = vertex3 - vertex1
        
        let sqLen1 = simd_dot(edge1, edge1)
        let sqLen2 = simd_dot(edge2, edge2)
        
        let perpDivisor = 0.5 / simd_cross(edge1, edge2).z
        
        let centerVec = DPoint(
            edge2.y * sqLen1 - edge1.y * sqLen2,
            edge1.x * sqLen2 - edge2.x * sqLen1
        ) * perpDivisor
        
        return simd_dot(centerVec, centerVec)
    }

    private func swap(arr: inout [UInt], index1: Int, index2: Int) {
        let temp = arr[index1]
        arr[index1] = arr[index2]
        arr[index2] = temp
    }
}
