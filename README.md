CC-ICP: Cylinder-Constrained ICP for Forest Point Cloud Registration
This is the official implementation of the Cylinder-Constrained Iterative Closest Point (CC-ICP) algorithm for fine registration of forest point clouds.

Forest point clouds obtained from multi-scan terrestrial laser scanning (TLS) suffer from repetitive geometries, severe occlusions, noise, and limited overlap, making conventional ICP-based fine registration methods inadequate. CC-ICP addresses these challenges by jointly leveraging ground plane constraints (point-to-plane residuals for vertical alignment) and tree trunk cylinder constraints (point-to-point residuals from cylinder axes for horizontal alignment) within a unified weighted least-squares optimization framework.

The method first extracts ground points via grid-minimum filtering and detects tree trunks as cylinders using DBSCAN clustering followed by RANSAC-based cylinder fitting. Cylinder axis correspondences between source and target scans are established through a projection-based optimization that refines circle centers in the projected plane. Both residual types are then integrated into a joint objective function and solved via the Gauss–Newton method.

Reference: [Author names], "CC-ICP: Cylinder-Constrained Iterative Closest Point for Fine Registration of Forest Point Clouds," [Journal], Year.

Build
Dependencies
The current implementation depends on the following libraries:

Point Cloud Library (PCL) (>= 1.8)

Eigen (>= 3.3)

A C++ compiler supporting C++14 or later

CMake (>= 3.12)

Building with CMake
On Linux or macOS:

bash
cd path-to-CC-ICP
mkdir build && cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make
On Windows with Microsoft Visual Studio, use the x64 Native Tools Command Prompt:

bash
cd path-to-CC-ICP
mkdir build && cd build
cmake -G "NMake Makefiles" -DCMAKE_BUILD_TYPE=Release ..
nmake
Alternatively, you can use any IDE that supports CMake (e.g., CLion, Visual Studio) to open the CMakeLists.txt file directly.

Usage
A basic usage example is provided in main.cpp. The program reads two PCD point clouds (source and target), performs ground extraction, cylinder detection, and runs the CC-ICP fine registration.

Quick Start
cpp
#include "cc_icp.h"

// Load point clouds
PointCloudT::Ptr source(new PointCloudT);
PointCloudT::Ptr target(new PointCloudT);
readPointCloud("source.pcd", source);
readPointCloud("target.pcd", target);

// Configure CC-ICP parameters
CCICPParams params;
params.slice_height = 0.1f;        // height interval for slicing (m)
params.max_iterations = 100;       // maximum ICP iterations
params.tolerance = 1e-5;           // convergence threshold
params.corr_dist = 0.5;            // correspondence search radius (m)
params.resolution = 0.1;           // grid size for ground extraction (m)
params.num_cylinder = 10;          // minimum points per cylinder cluster

// Run registration
CCICP cc_icp(params);
Eigen::Matrix4f transform = cc_icp.align(source, target);
Key Parameters
Parameter	Description	Default
slice_height	Height interval for point cloud slicing (m)	0.1
max_iterations	Maximum number of ICP iterations	100
tolerance	Convergence threshold for termination	1e-5
corr_dist	Maximum distance for establishing correspondences (m)	0.5
resolution	Grid cell size for ground extraction (m)	0.1
num_cylinder	Minimum number of points to validate a cylinder cluster	10
threshold_radius	Maximum radius difference for cylinder correspondence (m)	0.05
threshold_axisdirection	Maximum axis angle difference for cylinder correspondence (deg)	10.0
Input / Output
Input: Two PCD point cloud files (source and target scans)

Output: Transformed source point cloud aligned to the target, with the final 4×4 transformation matrix printed to console

Example
bash
./cc_icp
The default implementation reads 300-1-0614.pcd (source) and 300-2-0614.pcd (target). Modify the file paths in main.cpp according to your data.

Experimental Results
Extensive experiments were conducted on 24 plots across two mountainous forests (Guanshan and Tianmushan). CC-ICP consistently outperforms five widely-used algorithms—standard ICP, Point-to-Plane ICP, Generalized ICP (GICP), Symmetric ICP, and Normal Distributions Transform (NDT)—in terms of robustness, accuracy, and computational efficiency.

The method demonstrates strong performance under diverse forest conditions, including varying terrain slopes, tree species, and point densities.

Citation
If you find this work useful for your research, please cite:

bibtex
@article{...,
  title   = {CC-ICP: Cylinder-Constrained Iterative Closest Point for Fine Registration of Forest Point Clouds},
  author  = {...},
  journal = {...},
  year    = {...},
  volume  = {...},
  pages   = {...}
}
License
This project is licensed under the MIT License - see the LICENSE file for details.

