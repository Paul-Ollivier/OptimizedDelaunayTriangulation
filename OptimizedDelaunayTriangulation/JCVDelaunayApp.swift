import SwiftUI
import simd


// MARK: - SwiftUI View

struct ContentView: View {
    @State private var points: [Point] = []
    @State private var triangulation: Triangulation = Triangulation()
    @State private var triangulationTime: Double = 0.0  // Triangulation time in milliseconds
    @State private var numPoints: Int = 0               // Number of points generated each frame
    
    @State private var delaunay = JCVDelaunay(maxPoints: 1000) // Initialize once and reuse

    private let frameRate = 30.0 // Target 30 FPS

    var body: some View {
        VStack {
            Canvas { context, size in
                // Draw the triangles
                if !triangulation.triangles.isEmpty {
                    for i in stride(from: 0, to: triangulation.triangles.count, by: 3) {
                        let index0 = triangulation.triangles[i]
                        let index1 = triangulation.triangles[i + 1]
                        let index2 = triangulation.triangles[i + 2]
                        
                        let p0 = CGPoint(
                            x: points[index0].x * size.width,
                            y: points[index0].y * size.height
                        )
                        let p1 = CGPoint(
                            x: points[index1].x * size.width,
                            y: points[index1].y * size.height
                        )
                        let p2 = CGPoint(
                            x: points[index2].x * size.width,
                            y: points[index2].y * size.height
                        )
                        
                        var path = Path()
                        path.move(to: p0)
                        path.addLine(to: p1)
                        path.addLine(to: p2)
                        path.closeSubpath()
                        
                        context.stroke(path, with: .color(.blue), lineWidth: 0.5)
                    }
                }
                
                // Draw the points
                for point in points {
                    let p = CGPoint(
                        x: point.x * size.width,
                        y: point.y * size.height
                    )
                    let rect = CGRect(x: p.x - 2, y: p.y - 2, width: 4, height: 4)
                    context.fill(Path(ellipseIn: rect), with: .color(.red))
                }
            }
            .frame(maxWidth: .infinity, maxHeight: .infinity)
            .background(Color.white)
            
            // Display number of points and triangulation time
            Text("n points: \(numPoints) - triangulation time: \(String(format: "%.2f", triangulationTime)) ms")
                .font(.footnote)
                .padding()
        }
        .onAppear {
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
        // Generate between 1 and 1000 random points
        numPoints = Int.random(in: 1...1000)
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

// MARK: - Main App

@main
struct TriangulationApp: App {
    var body: some Scene {
        WindowGroup {
            ContentView()
        }
    }
}
