import SwiftUI

@main
struct JCVDelaunayApp: App {
    var body: some Scene {
        WindowGroup {
            ContentView()
        }
    }
}

struct ContentView: View {
    @State private var points: [Point] = []
    @State private var triangulation: Triangulation?
    @State private var executionTime: Double = 0
    @State private var isRunning = false
    @State private var timer: Timer?
    @State private var frameCount = 0
    
    var body: some View {
        VStack {
            ZStack {
                if let triangulation = triangulation {
                    TriangulationView(triangulation: triangulation)
                        .frame(maxWidth: .infinity, maxHeight: .infinity)
                        .border(Color.gray)
                } else {
                    Text("No triangulation yet")
                        .frame(maxWidth: .infinity, maxHeight: .infinity)
                        .border(Color.gray)
                }
                
                VStack {
                    Text("FPS: \(Int(1000/max(1, executionTime)))")
                        .font(.headline)
                    Text("Points: \(points.count)")
                        .font(.headline)
                    Text("Execution Time: \(String(format: "%.2f", executionTime)) ms")
                        .font(.headline)
                }
                .padding()
                .background(Color.white.opacity(0.8))
                .cornerRadius(8)
            }
        }
        .onAppear {
            startAnimation()
        }
    }
    
    private func startAnimation() {
        isRunning = true
        timer = Timer.scheduledTimer(withTimeInterval: 1.0/30.0, repeats: true) { _ in
            generateAndTriangulate()
        }
    }
    
    private func generateAndTriangulate() {
        // Generate random points
        var newPoints: [Point] = []
        for _ in 0..<1000 {
            let x = Double.random(in: 0...1000)
            let y = Double.random(in: 0...1000)
            newPoints.append(Point(x: x, y: y))
        }
        
        // Measure triangulation time
        let start = DispatchTime.now()
        let delaunay = JCVDelaunay(from: newPoints)
        let end = DispatchTime.now()
        
        let nanoTime = end.uptimeNanoseconds - start.uptimeNanoseconds
        let timeInMS = Double(nanoTime) / 1_000_000 // Convert to milliseconds
        
        // Create triangulation result
        let result = Triangulation(using: delaunay, with: newPoints)
        
        self.points = newPoints
        self.triangulation = result
        self.executionTime = timeInMS
    }
}

struct TriangulationView: View {
    let triangulation: Triangulation
    
    // Color configuration for visualization
    private let triangleStrokeColor = Color.blue.opacity(0.6)
    private let triangleFillColor = Color.blue.opacity(0.1)
    private let pointColor = Color.red
    private let hullColor = Color.green.opacity(0.8)
    
    // Drawing constants
    private let pointSize: CGFloat = 3
    private let lineWidth: CGFloat = 0.5
    private let hullWidth: CGFloat = 1.5
    
    var body: some View {
        GeometryReader { geometry in
            Canvas { context, size in
                // Calculate scaling factors to fit the view
                let scaleX = size.width / 1000
                let scaleY = size.height / 1000
                let scale = min(scaleX, scaleY)
                
                // Center the drawing
                let offsetX = (size.width - 1000 * scale) / 2
                let offsetY = (size.height - 1000 * scale) / 2
                
                // Transform point coordinates
                func transformPoint(_ point: Point) -> CGPoint {
                    CGPoint(
                        x: point.x * scale + offsetX,
                        y: point.y * scale + offsetY
                    )
                }
                
                // Draw triangles
                for i in stride(from: 0, to: triangulation.triangles.count, by: 3) {
                    let p1 = transformPoint(triangulation.points[triangulation.triangles[i]])
                    let p2 = transformPoint(triangulation.points[triangulation.triangles[i + 1]])
                    let p3 = transformPoint(triangulation.points[triangulation.triangles[i + 2]])
                    
                    var path = Path()
                    path.move(to: p1)
                    path.addLine(to: p2)
                    path.addLine(to: p3)
                    path.closeSubpath()
                    
                    // Fill triangle with a light color
                    context.fill(path, with: .color(triangleFillColor))
                    // Stroke triangle edges
                    context.stroke(path, with: .color(triangleStrokeColor), lineWidth: lineWidth)
                }
                
                // Draw convex hull
                if !triangulation.hull.isEmpty {
                    var hullPath = Path()
                    let firstPoint = transformPoint(triangulation.points[triangulation.hull[0]])
                    hullPath.move(to: firstPoint)
                    
                    for i in 1..<triangulation.hull.count {
                        let point = transformPoint(triangulation.points[triangulation.hull[i]])
                        hullPath.addLine(to: point)
                    }
                    hullPath.addLine(to: firstPoint)
                    
                    context.stroke(hullPath, with: .color(hullColor), lineWidth: hullWidth)
                }
                
                // Draw points
                for point in triangulation.points {
                    let transformedPoint = transformPoint(point)
                    let rect = CGRect(
                        x: transformedPoint.x - pointSize/2,
                        y: transformedPoint.y - pointSize/2,
                        width: pointSize,
                        height: pointSize
                    )
                    context.fill(Path(ellipseIn: rect), with: .color(pointColor))
                }
            }
        }
    }
}
