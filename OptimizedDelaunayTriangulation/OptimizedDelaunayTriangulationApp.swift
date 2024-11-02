import SwiftUI


@main
struct DelaunatorApp: App {
    var body: some Scene {
        WindowGroup {
            DelaunatorView()
        }
    }
}


struct DelaunatorView: View {
    @State private var points: [Point] = []
    @State private var triangles: [UInt32] = []
    @State private var lastExecutionTime: Double = 0
    private let size: CGFloat = 600
    private let pointCount = 10000
    
    var body: some View {
        VStack {
            Canvas { context, size in
                // Draw triangles
                for i in stride(from: 0, to: triangles.count, by: 3) {
                    if i + 2 < triangles.count {
                        let p1 = points[Int(triangles[i])]
                        let p2 = points[Int(triangles[i + 1])]
                        let p3 = points[Int(triangles[i + 2])]
                        
                        let path = Path { path in
                            path.move(to: CGPoint(x: p1.x, y: p1.y))
                            path.addLine(to: CGPoint(x: p2.x, y: p2.y))
                            path.addLine(to: CGPoint(x: p3.x, y: p3.y))
                            path.closeSubpath()
                        }
                        
                        context.stroke(path, with: .color(.blue.opacity(0.5)), lineWidth: 0.5)
                    }
                }
                
                // Draw points
                for point in points {
                    let rect = CGRect(x: point.x - 1, y: point.y - 1, width: 2, height: 2)
                    context.fill(Path(ellipseIn: rect), with: .color(.red.opacity(0.5)))
                }
            }
            .frame(width: size, height: size)
            .border(.gray, width: 1)
            
            Text("Points: \(pointCount)")
                .padding(.top)
            
            Text("Triangulation time: \(String(format: "%.2f", lastExecutionTime)) ms")
                .padding(.bottom)
            
            Button(action: generatePoints) {
                Text("Generate New Points")
                    .padding()
                    .background(Color.blue)
                    .foregroundColor(.white)
                    .cornerRadius(8)
            }
            .padding()
        }
        .onAppear {
            generatePoints()
        }
    }
    
    private func generatePoints() {
        // Generate random points
        var newPoints: [Point] = []
        for _ in 0..<pointCount {
            let x = Double.random(in: 50..<(size-50))
            let y = Double.random(in: 50..<(size-50))
            newPoints.append(Point(x: x, y: y))
        }
        
        // Measure triangulation time
        let startTime = CFAbsoluteTimeGetCurrent()
        
        let delaunator = Delaunator.from(points: newPoints)
        
        let endTime = CFAbsoluteTimeGetCurrent()
        lastExecutionTime = (endTime - startTime) * 1000 // Convert to milliseconds
        
        points = newPoints
        triangles = delaunator.triangles
        
        print("Generated \(pointCount) points")
        print("Number of triangles: \(triangles.count / 3)")
        print("Triangulation time: \(String(format: "%.2f", lastExecutionTime)) ms")
    }
}
