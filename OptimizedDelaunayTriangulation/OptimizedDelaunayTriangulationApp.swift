import SwiftUI
import simd

#if os(iOS)
import UIKit
typealias PlatformColor = UIColor
#else
import AppKit
typealias PlatformColor = NSColor
#endif

class TriangulationModel: ObservableObject {
    @Published var points: [Point] = []
    @Published var triangles: [UInt] = []
    @Published var triangulationTime: TimeInterval = 0
    @Published var pointCount: Int = 0
    
    let delaunator: CachedDelaunator
    private let maxPoints = 1000
    private var timer: Timer?
    private let bounds: CGRect
    
    init(bounds: CGRect) {
        self.bounds = bounds
        self.delaunator = CachedDelaunator(maxPoints: maxPoints)
        startAnimation()
    }
    
    func startAnimation() {
        timer = Timer.scheduledTimer(withTimeInterval: 1.0/30.0, repeats: true) { [weak self] _ in
            self?.updatePoints()
        }
    }
    
    private func updatePoints() {
        // Randomly choose number of points between 0 and maxPoints
        pointCount = Int.random(in: 0...maxPoints)
        
        // Generate new random points
        points = (0..<pointCount).map { _ in
            Point(DPoint(
                Double.random(in: bounds.minX...bounds.maxX),
                Double.random(in: bounds.minY...bounds.maxY)
            ))
        }
        
        updateTriangulation()
    }
    
    private func updateTriangulation() {
        let startTime = CACurrentMediaTime()
        delaunator.update(from: points)
        let endTime = CACurrentMediaTime()
        
        triangulationTime = endTime - startTime
        triangles = delaunator.triangles
    }
    
    deinit {
        timer?.invalidate()
    }
}

struct TriangulationView: View {
    @StateObject private var model: TriangulationModel
    
    init(bounds: CGRect) {
        _model = StateObject(wrappedValue: TriangulationModel(bounds: bounds))
    }
    
    var body: some View {
        ZStack {
            Canvas { context, size in
                context.withCGContext { ctx in
                    ctx.setStrokeColor(PlatformColor.white.cgColor)
                    ctx.setLineWidth(0.5)
                    
                    for i in stride(from: 0, to: model.triangles.count, by: 3) {
                        guard i + 2 < model.triangles.count else { break }
                        
                        let p1 = model.points[Int(model.triangles[i])].coords
                        let p2 = model.points[Int(model.triangles[i + 1])].coords
                        let p3 = model.points[Int(model.triangles[i + 2])].coords
                        
                        ctx.move(to: CGPoint(x: p1.x, y: p1.y))
                        ctx.addLine(to: CGPoint(x: p2.x, y: p2.y))
                        ctx.addLine(to: CGPoint(x: p3.x, y: p3.y))
                        ctx.closePath()
                    }
                    ctx.strokePath()
                    
                    ctx.setFillColor(PlatformColor.red.cgColor)
                    for point in model.points {
                        let rect = CGRect(
                            x: point.coords.x - 1,
                            y: point.coords.y - 1,
                            width: 2,
                            height: 2
                        )
                        ctx.fillEllipse(in: rect)
                    }
                }
            }
            .background(Color.black)
            .ignoresSafeArea()
            
            VStack {
                Spacer()
                VStack(spacing: 8) {
                    Text("Points: \(model.pointCount)")
                        .foregroundColor(.white)
                        .font(.system(.title2, design: .monospaced))
                    Text("Triangulation Time: \(String(format: "%.2f", model.triangulationTime * 1000)) ms")
                        .foregroundColor(.white)
                        .font(.system(.title2, design: .monospaced))
                }
                .padding()
                .background(Color.black.opacity(0.7))
                .cornerRadius(10)
                .padding(.bottom)
            }
        }
    }
}

struct ContentView: View {
    var body: some View {
        GeometryReader { geometry in
            TriangulationView(bounds: CGRect(
                x: 20,
                y: 20,
                width: geometry.size.width - 40,
                height: geometry.size.height - 40
            ))
        }
    }
}

@main
struct DelaunayApp: App {
    var body: some Scene {
        WindowGroup {
            ContentView()
        }
    }
}
