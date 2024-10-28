import SwiftUI
import MetalKit

// Main SwiftUI view for the app
@main
struct OptimizedDelaunayTriangulationApp: App {
    var body: some Scene {
        WindowGroup {
            ContentView()
        }
    }
}

// ContentView for displaying the Metal rendering
struct ContentView: View {
    var body: some View {
        MetalView()
            .frame(maxWidth: .infinity, maxHeight: .infinity)
    }
}

// MetalView for Metal rendering, conforming to NSViewRepresentable
struct MetalView: NSViewRepresentable {

    func makeCoordinator() -> Coordinator {
        Coordinator(self)
    }

    func makeNSView(context: Context) -> MTKView {
        let metalView = MTKView()
        guard let device = MTLCreateSystemDefaultDevice() else {
            fatalError("Metal is not supported on this device")
        }
        
        metalView.device = device
        metalView.delegate = context.coordinator
        metalView.clearColor = MTLClearColor(red: 0, green: 0, blue: 0, alpha: 1)
        
        // Passing the device to the coordinator to avoid context construction
        context.coordinator.device = device
        
        // Setup Metal and Triangulation
        context.coordinator.setupMetal()
        context.coordinator.setupTriangulation()
        
        return metalView
    }

    func updateNSView(_ nsView: MTKView, context: Context) {}
    
    // Coordinator class to manage Metal setup and rendering
    class Coordinator: NSObject, MTKViewDelegate {
        var parent: MetalView
        var commandQueue: MTLCommandQueue?
        var pipelineState: MTLRenderPipelineState?
        var vertexBuffer: MTLBuffer?
        var indexBuffer: MTLBuffer?
        var colorBuffer: MTLBuffer? // Buffer to hold colors for each vertex
        var device: MTLDevice? // Store the device instance here
        
        init(_ parent: MetalView) {
            self.parent = parent
            super.init()
        }

        // MTKViewDelegate - Handles drawing content
        func draw(in view: MTKView) {
            guard let drawable = view.currentDrawable,
                  let commandQueue = commandQueue,
                  let renderPassDescriptor = view.currentRenderPassDescriptor,
                  let pipelineState = pipelineState,
                  let vertexBuffer = vertexBuffer,
                  let indexBuffer = indexBuffer,
                  let colorBuffer = colorBuffer else {
                return
            }

            let commandBuffer = commandQueue.makeCommandBuffer()
            let renderEncoder = commandBuffer?.makeRenderCommandEncoder(descriptor: renderPassDescriptor)
            renderEncoder?.setRenderPipelineState(pipelineState)
            renderEncoder?.setVertexBuffer(vertexBuffer, offset: 0, index: 0)
            renderEncoder?.setVertexBuffer(colorBuffer, offset: 0, index: 1)
            renderEncoder?.drawIndexedPrimitives(type: .triangle, indexCount: indexBuffer.length / MemoryLayout<UInt16>.stride, indexType: .uint16, indexBuffer: indexBuffer, indexBufferOffset: 0)
            renderEncoder?.endEncoding()

            commandBuffer?.present(drawable)
            commandBuffer?.commit()
        }

        func mtkView(_ view: MTKView, drawableSizeWillChange size: CGSize) {}

        // Configure Metal pipeline and setup the triangulation
        func setupMetal() {
            guard let device = device else { return }
            commandQueue = device.makeCommandQueue()

            let library = device.makeDefaultLibrary()
            let vertexFunction = library?.makeFunction(name: "vertex_main")
            let fragmentFunction = library?.makeFunction(name: "fragment_main")

            let pipelineDescriptor = MTLRenderPipelineDescriptor()
            pipelineDescriptor.vertexFunction = vertexFunction
            pipelineDescriptor.fragmentFunction = fragmentFunction
            pipelineDescriptor.colorAttachments[0].pixelFormat = .bgra8Unorm

            // Setup the vertex descriptor
            let vertexDescriptor = MTLVertexDescriptor()
            vertexDescriptor.attributes[0].format = .float2    // Position attribute
            vertexDescriptor.attributes[0].offset = 0
            vertexDescriptor.attributes[0].bufferIndex = 0

            vertexDescriptor.attributes[1].format = .float4    // Color attribute
            vertexDescriptor.attributes[1].offset = 0
            vertexDescriptor.attributes[1].bufferIndex = 1

            vertexDescriptor.layouts[0].stride = MemoryLayout<SIMD2<Float>>.stride
            vertexDescriptor.layouts[1].stride = MemoryLayout<SIMD4<Float>>.stride
            pipelineDescriptor.vertexDescriptor = vertexDescriptor

            do {
                pipelineState = try device.makeRenderPipelineState(descriptor: pipelineDescriptor)
            } catch {
                print("Failed to create pipeline state, error \(error)")
            }
        }

        // Setup triangulation and create Metal buffers
        func setupTriangulation() {
            // Generate random points in screen space
            var randomPoints = [DPoint]()
            
            for _ in 0..<10000 {
                let x = CGFloat.random(in: -1.0...1.0)
                let y = CGFloat.random(in: -1.0...1.0)
                randomPoints.append(DPoint(x: x, y: y))
            }
            
            print("Generated \(randomPoints.count) random points")
            
            let triangulationStartTime = DispatchTime.now()
            let triangles = triangulate(randomPoints)
            let triangulationEndTime = DispatchTime.now()
                    
            let triangulationTime = Double(triangulationEndTime.uptimeNanoseconds - triangulationStartTime.uptimeNanoseconds) / 1_000_000
            print("Triangulation Time: \(triangulationTime) ms")

            // Convert triangulation points to Metal-compatible vertices and build point map
            var vertices = [SIMD2<Float>]()
            var pointIndexMap = [DPoint: UInt16]()
            var currentIndex: UInt16 = 0
            
            // Create a set of all unique points from triangles
            var uniquePoints = Set<DPoint>()
            for triangle in triangles {
                uniquePoints.insert(triangle.point1)
                uniquePoints.insert(triangle.point2)
                uniquePoints.insert(triangle.point3)
            }
            
            // First, map all points from triangles
            for point in uniquePoints {
                let vertex = SIMD2<Float>(Float(point.x), Float(point.y))
                vertices.append(vertex)
                pointIndexMap[point] = currentIndex
                currentIndex += 1
            }
            
            print("Created map with \(pointIndexMap.count) unique points")
            
            // Create index buffer from Delaunay triangles with error checking
            var indices = [UInt16]()
            var missingPointCount = 0
            
            for (i, triangle) in triangles.enumerated() {
                let index1 = pointIndexMap[triangle.point1]
                let index2 = pointIndexMap[triangle.point2]
                let index3 = pointIndexMap[triangle.point3]
                
                if let idx1 = index1, let idx2 = index2, let idx3 = index3 {
                    indices.append(idx1)
                    indices.append(idx2)
                    indices.append(idx3)
                } else {
                    missingPointCount += 1
                    print("Missing point in triangle \(i):")
                    print("  Point1: \(triangle.point1) -> \(String(describing: index1))")
                    print("  Point2: \(triangle.point2) -> \(String(describing: index2))")
                    print("  Point3: \(triangle.point3) -> \(String(describing: index3))")
                }
            }
            
            print("Generated \(indices.count / 3) complete triangles")
            if missingPointCount > 0 {
                print("Warning: \(missingPointCount) triangles had missing points")
            }
            
            // Ensure we have enough vertices and indices
            guard indices.count > 0 else {
                print("Error: No indices generated!")
                return
            }
            
            // Create random colors for the actual number of triangles we have
            var colors = [SIMD4<Float>]()
            let numTriangles = indices.count / 3
            for _ in 0..<numTriangles {
                let r = Float.random(in: 0.0...1.0)
                let g = Float.random(in: 0.0...1.0)
                let b = Float.random(in: 0.0...1.0)
                let color = SIMD4<Float>(r, g, b, 1.0)
                colors.append(contentsOf: [color, color, color])
            }
            
            guard let device = device else { return }
            
            // Create Metal buffers with size checking
            print("Creating Metal buffers:")
            print("  Vertices: \(vertices.count)")
            print("  Indices: \(indices.count)")
            print("  Colors: \(colors.count)")
            
            vertexBuffer = device.makeBuffer(bytes: vertices,
                                           length: vertices.count * MemoryLayout<SIMD2<Float>>.stride,
                                           options: [])
            indexBuffer = device.makeBuffer(bytes: indices,
                                          length: indices.count * MemoryLayout<UInt16>.stride,
                                          options: [])
            colorBuffer = device.makeBuffer(bytes: colors,
                                          length: colors.count * MemoryLayout<SIMD4<Float>>.stride,
                                          options: [])
        }
    }
}
