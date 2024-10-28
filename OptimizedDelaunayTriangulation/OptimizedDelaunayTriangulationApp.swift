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
        metalView.preferredFramesPerSecond = 30 // Set the frame rate to 30 fps
        
        // Passing the device to the coordinator to avoid context construction
        context.coordinator.device = device
        
        // Setup Metal pipeline
        context.coordinator.setupMetal()
        
        return metalView
    }

    func updateNSView(_ nsView: MTKView, context: Context) {}
    
    // Coordinator class to manage Metal setup and rendering
    class Coordinator: NSObject, MTKViewDelegate {
        var parent: MetalView
        var commandQueue: MTLCommandQueue?
        var pipelineState: MTLRenderPipelineState?
        var device: MTLDevice? // Store the device instance here
        
        // Buffers to hold vertex and index data
        var vertexBuffer: MTLBuffer?
        var indexBuffer: MTLBuffer?
        var colorBuffer: MTLBuffer?
        
        // Preallocate arrays to avoid frequent allocations
        var vertices = [SIMD2<Float>]()
        var indices = [UInt16]()
        var colors = [SIMD4<Float>]()
        
        init(_ parent: MetalView) {
            self.parent = parent
            super.init()
        }

        // MTKViewDelegate - Handles drawing content
        func draw(in view: MTKView) {
            autoreleasepool {
                guard let drawable = view.currentDrawable,
                      let commandQueue = commandQueue,
                      let renderPassDescriptor = view.currentRenderPassDescriptor,
                      let pipelineState = pipelineState,
                      let device = device else {
                    return
                }

                let frameStartTime = DispatchTime.now()
                var processTimes = [String: Double]()
                
                // Step 1: Generate random points
                let startTime1 = DispatchTime.now()
                let pointCount = Int.random(in: 100...1000)
                var randomPoints = [DPoint]()
                randomPoints.reserveCapacity(pointCount)
                
                for _ in 0..<pointCount {
                    let x = Double.random(in: -1.0...1.0)
                    let y = Double.random(in: -1.0...1.0)
                    randomPoints.append(DPoint(x: x, y: y))
                }
                let endTime1 = DispatchTime.now()
                processTimes["Generating Random Points"] = Double(endTime1.uptimeNanoseconds - startTime1.uptimeNanoseconds) * 1e-6
                
                // Step 2: Perform triangulation
                let startTime2 = DispatchTime.now()
                let triangles = triangulate(randomPoints)
                let endTime2 = DispatchTime.now()
                processTimes["Triangulation"] = Double(endTime2.uptimeNanoseconds - startTime2.uptimeNanoseconds) * 1e-6
                
                // Step 3: Prepare vertex and index data
                let startTime3 = DispatchTime.now()
                self.vertices.removeAll(keepingCapacity: true)
                self.indices.removeAll(keepingCapacity: true)
                self.colors.removeAll(keepingCapacity: true)
                
                var pointIndexMap = [DPoint: UInt16]()
                var currentIndex: UInt16 = 0
                
                // Map points to vertices
                for triangle in triangles {
                    for point in [triangle.point1, triangle.point2, triangle.point3] {
                        if pointIndexMap[point] == nil {
                            let vertex = SIMD2<Float>(Float(point.x), Float(point.y))
                            self.vertices.append(vertex)
                            pointIndexMap[point] = currentIndex
                            currentIndex += 1
                        }
                    }
                }
                
                // Create index buffer
                for triangle in triangles {
                    if let idx1 = pointIndexMap[triangle.point1],
                       let idx2 = pointIndexMap[triangle.point2],
                       let idx3 = pointIndexMap[triangle.point3] {
                        self.indices.append(idx1)
                        self.indices.append(idx2)
                        self.indices.append(idx3)
                    }
                }
                let endTime3 = DispatchTime.now()
                processTimes["Preparing Vertex and Index Data"] = Double(endTime3.uptimeNanoseconds - startTime3.uptimeNanoseconds) * 1e-6
                
                // Step 4: Create random colors
                let startTime4 = DispatchTime.now()
                let numTriangles = self.indices.count / 3
                for _ in 0..<numTriangles {
                    let r = Float.random(in: 0.0...1.0)
                    let g = Float.random(in: 0.0...1.0)
                    let b = Float.random(in: 0.0...1.0)
                    let color = SIMD4<Float>(r, g, b, 1.0)
                    self.colors.append(contentsOf: [color, color, color])
                }
                let endTime4 = DispatchTime.now()
                processTimes["Creating Random Colors"] = Double(endTime4.uptimeNanoseconds - startTime4.uptimeNanoseconds) * 1e-6
                
                // Step 5: Create or update Metal buffers
                let startTime5 = DispatchTime.now()
                self.vertexBuffer = device.makeBuffer(bytes: self.vertices,
                                                      length: self.vertices.count * MemoryLayout<SIMD2<Float>>.stride,
                                                      options: .storageModeShared)
                self.indexBuffer = device.makeBuffer(bytes: self.indices,
                                                     length: self.indices.count * MemoryLayout<UInt16>.stride,
                                                     options: .storageModeShared)
                self.colorBuffer = device.makeBuffer(bytes: self.colors,
                                                     length: self.colors.count * MemoryLayout<SIMD4<Float>>.stride,
                                                     options: .storageModeShared)
                let endTime5 = DispatchTime.now()
                processTimes["Creating Metal Buffers"] = Double(endTime5.uptimeNanoseconds - startTime5.uptimeNanoseconds) * 1e-6
                
                // Step 6: Render
                let startTime6 = DispatchTime.now()
                let commandBuffer = commandQueue.makeCommandBuffer()
                let renderEncoder = commandBuffer?.makeRenderCommandEncoder(descriptor: renderPassDescriptor)
                renderEncoder?.setRenderPipelineState(pipelineState)
                renderEncoder?.setVertexBuffer(self.vertexBuffer, offset: 0, index: 0)
                renderEncoder?.setVertexBuffer(self.colorBuffer, offset: 0, index: 1)
                
                // Correct indexCount calculation
                let indexCount = self.indices.count
                
                if indexCount > 0 {
                    renderEncoder?.drawIndexedPrimitives(type: .triangle,
                                                         indexCount: indexCount,
                                                         indexType: .uint16,
                                                         indexBuffer: self.indexBuffer!,
                                                         indexBufferOffset: 0)
                }
                
                renderEncoder?.endEncoding()
                
                commandBuffer?.present(drawable)
                commandBuffer?.commit()
                let endTime6 = DispatchTime.now()
                processTimes["Rendering"] = Double(endTime6.uptimeNanoseconds - startTime6.uptimeNanoseconds) * 1e-6
                
                let frameEndTime = DispatchTime.now()
                processTimes["Total Frame Time"] = Double(frameEndTime.uptimeNanoseconds - frameStartTime.uptimeNanoseconds) * 1e-6
                
                // Optional: Print process times for debugging
                // print("Process Times (ms):", processTimes)
                
                // Optional: Log frame rate
                // print("Frame Time: \(processTimes["Total Frame Time"]!) ms")
            }
        }

        func mtkView(_ view: MTKView, drawableSizeWillChange size: CGSize) {}
        
        // Configure Metal pipeline
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
    }
}
