import SwiftUI
import MetalKit
import simd

@main
struct TriangulationApp: App {
    var body: some Scene {
        WindowGroup {
            ContentView()
        }
    }
}


struct ContentView: View {
    
    private static let MAX_POINTS = 1 << 11
    
    @State private var points: [Point] = []
    @State private var triangulation: Triangulation = Triangulation()
    @State private var triangulationTime: Double = 0.0  // Triangulation time in milliseconds
    @State private var numPoints: Int = 2
    
    @State private var delaunay = OptimizedDelaunay(maxPoints: MAX_POINTS) // Initialize once and reuse
    
    private let frameRate = 60.0
    
    var body: some View {
        VStack {
            MetalView(triangulation: $triangulation, points: $points)
                .frame(maxWidth: .infinity, maxHeight: .infinity)
                .background(Color.black)
            
            // Display number of points and triangulation time
            Text("Points: \(numPoints) - Triangulation time: \(String(format: "%.2f", triangulationTime)) ms")
                .font(.footnote)
                .padding()
        }
        .onAppear {
            generatePoints() // Initialize data before rendering starts
            startTimer()
        }
    }
    
    // Timer for updating points and triangulation every frame
    private func startTimer() {
        Timer.scheduledTimer(withTimeInterval: 1.0 / frameRate, repeats: true) { _ in
            generatePoints()
        }
    }
    
    func generatePoints() {
        numPoints += 1
        if numPoints >= ContentView.MAX_POINTS {
            numPoints = 3
        }
        // numPoints = Int.random(in: 1...ContentView.MAX_POINTS)
        // numPoints = ContentView.MAX_POINTS
        points = (0..<numPoints).map { _ in
            Point(x: Double.random(in: 0...1), y: Double.random(in: 0...1))
        }
        
        // Measure the time taken for triangulation
        let start = DispatchTime.now()
        
        // Reuse the delaunay instance
        delaunay.triangulate(points: points)
        
        let end = DispatchTime.now()
        let nanoTime = end.uptimeNanoseconds - start.uptimeNanoseconds
        triangulationTime = Double(nanoTime) / 1_000_000  // Convert to milliseconds
        
        // Update triangulation result
        triangulation = Triangulation(using: delaunay, with: points)
    }
}


struct MetalView: NSViewRepresentable {
    @Binding var triangulation: Triangulation
    @Binding var points: [Point]
    
    func makeCoordinator() -> Coordinator {
        Coordinator(triangulation: $triangulation, points: $points)
    }
    
    func makeNSView(context: Context) -> MTKView {
        let mtkView = MTKView()
        mtkView.device = MTLCreateSystemDefaultDevice()
        mtkView.delegate = context.coordinator
        mtkView.preferredFramesPerSecond = 30
        mtkView.enableSetNeedsDisplay = false
        mtkView.framebufferOnly = false
        mtkView.colorPixelFormat = .rgba8Unorm // Set pixel format to RGBA
        mtkView.clearColor = MTLClearColorMake(0, 0, 0, 1) // Black background
        return mtkView
    }
    
    func updateNSView(_ nsView: MTKView, context: Context) {
        context.coordinator.update(triangulation: triangulation, points: points)
    }
    
    class Coordinator: NSObject, MTKViewDelegate {
        var triangulation: Triangulation
        var points: [Point]
        
        var device: MTLDevice!
        var commandQueue: MTLCommandQueue!
        var pipelineState: MTLRenderPipelineState!
        
        var vertexBuffer: MTLBuffer?
        var numVertices: Int = 0
        
        init(triangulation: Binding<Triangulation>, points: Binding<[Point]>) {
            self.triangulation = triangulation.wrappedValue
            self.points = points.wrappedValue
            super.init()
            
            device = MTLCreateSystemDefaultDevice()
            commandQueue = device.makeCommandQueue()
            buildPipelineState()
            updateBuffers()
        }
        
        func update(triangulation: Triangulation, points: [Point]) {
            self.triangulation = triangulation
            self.points = points
            updateBuffers()
        }
        
        func buildPipelineState() {
            let library = device.makeDefaultLibrary()
            let pipelineDescriptor = MTLRenderPipelineDescriptor()
            pipelineDescriptor.vertexFunction = library?.makeFunction(name: "vertexShader")
            pipelineDescriptor.fragmentFunction = library?.makeFunction(name: "fragmentShader")
            pipelineDescriptor.colorAttachments[0].pixelFormat = .rgba8Unorm // Ensure RGBA format
            
            // Create and configure the vertex descriptor
            let vertexDescriptor = MTLVertexDescriptor()
            
            // Calculate the correct offset for the color attribute
            let colorOffset = MemoryLayout<Vertex>.offset(of: \Vertex.color)!
            
            // Position attribute (attribute 0)
            vertexDescriptor.attributes[0].format = .float2
            vertexDescriptor.attributes[0].offset = 0
            vertexDescriptor.attributes[0].bufferIndex = 0
            
            // Color attribute (attribute 1)
            vertexDescriptor.attributes[1].format = .float4
            vertexDescriptor.attributes[1].offset = colorOffset // Corrected offset
            vertexDescriptor.attributes[1].bufferIndex = 0
            
            // Layout for the vertex buffer
            vertexDescriptor.layouts[0].stride = MemoryLayout<Vertex>.size
            vertexDescriptor.layouts[0].stepRate = 1
            vertexDescriptor.layouts[0].stepFunction = .perVertex
            
            pipelineDescriptor.vertexDescriptor = vertexDescriptor
            
            do {
                pipelineState = try device.makeRenderPipelineState(descriptor: pipelineDescriptor)
            } catch let error {
                print("Failed to create pipeline state: \(error)")
            }
        }
        
        func updateBuffers() {
            guard !triangulation.triangles.isEmpty else { return }
            
            var vertexData: [Vertex] = []
            
            // Process each triangle
            for i in stride(from: 0, to: triangulation.triangles.count, by: 3) {
                let index0 = triangulation.triangles[i]
                let index1 = triangulation.triangles[i + 1]
                let index2 = triangulation.triangles[i + 2]
                
                let p0 = points[index0]
                let p1 = points[index1]
                let p2 = points[index2]
                
                // Transform points from [0,1] to [-1,1] coordinate space
                let pos0 = SIMD2<Float>(Float(p0.x * 2 - 1), Float(p0.y * 2 - 1))
                let pos1 = SIMD2<Float>(Float(p1.x * 2 - 1), Float(p1.y * 2 - 1))
                let pos2 = SIMD2<Float>(Float(p2.x * 2 - 1), Float(p2.y * 2 - 1))
                
                // Generate random colors for the vertices (RGBA format)
                let color0 = SIMD4<Float>(Float.random(in: 0...1), // Red
                                          Float.random(in: 0...1), // Green
                                          Float.random(in: 0...1), // Blue
                                          1.0) // Alpha
                let color1 = SIMD4<Float>(Float.random(in: 0...1),
                                          Float.random(in: 0...1),
                                          Float.random(in: 0...1),
                                          1.0)
                let color2 = SIMD4<Float>(Float.random(in: 0...1),
                                          Float.random(in: 0...1),
                                          Float.random(in: 0...1),
                                          1.0)
                
                // Create vertices
                let vertex0 = Vertex(position: pos0, color: color0)
                let vertex1 = Vertex(position: pos1, color: color1)
                let vertex2 = Vertex(position: pos2, color: color2)
                
                vertexData.append(vertex0)
                vertexData.append(vertex1)
                vertexData.append(vertex2)
            }
            
            numVertices = vertexData.count
            
            // Create vertex buffer
            vertexBuffer = device.makeBuffer(bytes: vertexData,
                                             length: vertexData.count * MemoryLayout<Vertex>.size,
                                             options: [])
        }
        
        func draw(in view: MTKView) {
            guard let drawable = view.currentDrawable,
                  let renderPassDescriptor = view.currentRenderPassDescriptor,
                  let commandBuffer = commandQueue.makeCommandBuffer(),
                  let renderEncoder = commandBuffer.makeRenderCommandEncoder(descriptor: renderPassDescriptor),
                  let vertexBuffer = vertexBuffer else {
                return
            }
            
            renderEncoder.setRenderPipelineState(pipelineState)
            renderEncoder.setVertexBuffer(vertexBuffer, offset: 0, index: 0)
            
            // Draw the triangle mesh
            renderEncoder.drawPrimitives(type: .triangle,
                                         vertexStart: 0,
                                         vertexCount: numVertices)
            
            renderEncoder.endEncoding()
            commandBuffer.present(drawable)
            commandBuffer.commit()
        }
        
        func mtkView(_ view: MTKView, drawableSizeWillChange size: CGSize) {
            // Handle view size changes if necessary
        }
    }
}


struct Vertex {
    var position: SIMD2<Float>
    var color: SIMD4<Float>
}
