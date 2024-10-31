import SwiftUI
import MetalKit

@main
struct OptimizedDelaunayTriangulationApp: App {
    var body: some Scene {
        WindowGroup {
            ContentView()
        }
    }
}

struct ContentView: View {
    var body: some View {
        MetalView()
            .frame(maxWidth: .infinity, maxHeight: .infinity)
    }
}

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
        metalView.preferredFramesPerSecond = 30
        
        context.coordinator.device = device
        context.coordinator.setupMetal()
        
        Task { @MainActor in
            await context.coordinator.startRenderLoop()
        }
        
        return metalView
    }
    
    func updateNSView(_ nsView: MTKView, context: Context) {}
}

class Coordinator: NSObject, MTKViewDelegate {
    var parent: MetalView
    var device: MTLDevice?
    var commandQueue: MTLCommandQueue?
    var pipelineState: MTLRenderPipelineState?
    
    @MainActor private var currentVertexBuffer: MTLBuffer?
    @MainActor private var currentIndexBuffer: MTLBuffer?
    @MainActor private var currentColorBuffer: MTLBuffer?
    @MainActor private var currentIndexCount: Int = 0
    
    private var isProcessing = false
    private var needsUpdate = true
    
    init(_ parent: MetalView) {
        self.parent = parent
        super.init()
    }
    
    @MainActor
    func startRenderLoop() async {
        while true {
            if needsUpdate && !isProcessing {
                await generateNewTriangulation()
            }
            try? await Task.sleep(nanoseconds: UInt64(1.0 / 30.0 * Double(NSEC_PER_SEC)))
        }
    }
    
    private func generateNewTriangulation() async {
        isProcessing = true
        needsUpdate = false
        
        // let pointCount = Int.random(in: 100...1000)
        let pointCount = 1000
        var randomPoints = [DPoint]()
        randomPoints.reserveCapacity(pointCount)
        
        for _ in 0..<pointCount {
            let x = Double.random(in: -1.0...1.0)
            let y = Double.random(in: -1.0...1.0)
            randomPoints.append(DPoint(x: x, y: y))
        }
        
        let startTime = CFAbsoluteTimeGetCurrent()
        let triangles = await parallelTriangulate(randomPoints)
        let endTime = CFAbsoluteTimeGetCurrent()
        
        print("Triangulation completed in \((endTime - startTime) * 1000)ms for \(pointCount) points")
        
        await updateBuffers(with: triangles)
        isProcessing = false
    }
    
    private func updateBuffers(with triangles: [DTriangle]) async {
        guard let device = device else { return }
        
        var newVertices = [SIMD2<Float>]()
        var newIndices = [UInt16]()
        var newColors = [SIMD4<Float>]()
        
        var pointIndexMap = [DPoint: UInt16]()
        var currentIndex: UInt16 = 0
        
        for triangle in triangles {
            for point in [triangle.point1, triangle.point2, triangle.point3] {
                if pointIndexMap[point] == nil {
                    let vertex = SIMD2<Float>(Float(point.x), Float(point.y))
                    newVertices.append(vertex)
                    pointIndexMap[point] = currentIndex
                    currentIndex += 1
                }
            }
        }
        
        for triangle in triangles {
            if let idx1 = pointIndexMap[triangle.point1],
               let idx2 = pointIndexMap[triangle.point2],
               let idx3 = pointIndexMap[triangle.point3] {
                newIndices.append(idx1)
                newIndices.append(idx2)
                newIndices.append(idx3)
            }
        }
        
        let numTriangles = newIndices.count / 3
        for _ in 0..<numTriangles {
            let r = Float.random(in: 0.0...1.0)
            let g = Float.random(in: 0.0...1.0)
            let b = Float.random(in: 0.0...1.0)
            let color = SIMD4<Float>(r, g, b, 1.0)
            newColors.append(contentsOf: [color, color, color])
        }
        
        let vertexBufferSize = newVertices.count * MemoryLayout<SIMD2<Float>>.stride
        let indexBufferSize = newIndices.count * MemoryLayout<UInt16>.stride
        let colorBufferSize = newColors.count * MemoryLayout<SIMD4<Float>>.stride
        
        let newVertexBuffer = device.makeBuffer(bytes: newVertices, length: vertexBufferSize, options: .storageModeShared)
        let newIndexBuffer = device.makeBuffer(bytes: newIndices, length: indexBufferSize, options: .storageModeShared)
        let newColorBuffer = device.makeBuffer(bytes: newColors, length: colorBufferSize, options: .storageModeShared)
        
        await MainActor.run {
            self.currentVertexBuffer = newVertexBuffer
            self.currentIndexBuffer = newIndexBuffer
            self.currentColorBuffer = newColorBuffer
            self.currentIndexCount = newIndices.count
            
            print("Generated \(numTriangles) triangles from \(pointIndexMap.count) points")
        }
    }
    
    func draw(in view: MTKView) {
        guard let drawable = view.currentDrawable,
              let commandQueue = commandQueue,
              let renderPassDescriptor = view.currentRenderPassDescriptor,
              let pipelineState = pipelineState else {
            return
        }
        
        let vertexBuffer = currentVertexBuffer
        let indexBuffer = currentIndexBuffer
        let colorBuffer = currentColorBuffer
        let indexCount = currentIndexCount
        
        let commandBuffer = commandQueue.makeCommandBuffer()
        let renderEncoder = commandBuffer?.makeRenderCommandEncoder(descriptor: renderPassDescriptor)
        renderEncoder?.setRenderPipelineState(pipelineState)
        
        if let vertexBuffer = vertexBuffer,
           let indexBuffer = indexBuffer,
           let colorBuffer = colorBuffer,
           indexCount > 0 {
            
            renderEncoder?.setVertexBuffer(vertexBuffer, offset: 0, index: 0)
            renderEncoder?.setVertexBuffer(colorBuffer, offset: 0, index: 1)
            
            renderEncoder?.drawIndexedPrimitives(type: .triangle,
                                               indexCount: indexCount,
                                               indexType: .uint16,
                                               indexBuffer: indexBuffer,
                                               indexBufferOffset: 0)
        }
        
        renderEncoder?.endEncoding()
        commandBuffer?.present(drawable)
        commandBuffer?.commit()
        
        if !isProcessing {
            needsUpdate = true
        }
    }
    
    func mtkView(_ view: MTKView, drawableSizeWillChange size: CGSize) {}
    
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
        
        let vertexDescriptor = MTLVertexDescriptor()
        vertexDescriptor.attributes[0].format = .float2
        vertexDescriptor.attributes[0].offset = 0
        vertexDescriptor.attributes[0].bufferIndex = 0
        
        vertexDescriptor.attributes[1].format = .float4
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
