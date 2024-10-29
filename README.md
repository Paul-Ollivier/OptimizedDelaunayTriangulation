# Optimized Delaunay Triangulation Visualization

A high-performance, real-time Delaunay triangulation visualization tool built with Swift & Metal.

<img width="1792" alt="Screenshot 2024-10-29 at 23 33 48" src="https://github.com/user-attachments/assets/0da58c2d-9e96-442d-82b3-fb07a1b0ae7c">


## Features

- Real-time Delaunay triangulation computation ( CPU )
- Hardware-accelerated rendering using Metal

## Performance Optimizations

- Preallocation of buffers and collections
- SIMD operations for geometric calculations
- Binary search for ordered insertions
- Efficient memory reuse
- Autoreleasepool usage for consistent performance
- Performance monitoring and timing for each processing stage

## Requirements

- macOS 11.0 or later
- Xcode 13.0 or later
- Metal-capable Mac

## Building and Running

1. Clone the repository:
```bash
git clone https://github.com/yourusername/optimized-delaunay-triangulation.git
```

2. Open the project in Xcode:
```bash
cd optimized-delaunay-triangulation
open OptimizedDelaunayTriangulation.xcodeproj
```

3. Build and run the project (⌘R)

## Project Structure

```
OptimizedDelaunayTriangulation/
├── OptimizedDelaunayTriangulationApp.swift  # Main app entry point and SwiftUI setup
├── OptimizedDelaunayTriangulation.swift     # Core triangulation algorithm
└── Shaders.metal                            # Metal shader code
```

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request. For major changes, please open an issue first to discuss what you would like to change.

## License

[MIT](https://choosealicense.com/licenses/mit/)

## Acknowledgments

- The triangulation algorithm is based on the Bowyer-Watson algorithm
- Special thanks to the Metal and Swift communities for their resources and documentation

## Author

Paul Ollivier

## Contact

For questions and feedback:
- LinkedIn: https://www.linkedin.com/in/paulollivier/
