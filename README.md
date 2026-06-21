# CC-ICP: Cylinder-Constrained Iterative Closest Point for Forest Point Cloud Fine Registration

This repository contains the official C++ implementation of **CC-ICP (Cylinder-Constrained Iterative Closest Point)**, a dedicated fine registration algorithm for forest point clouds, as described in the paper. 

## Overview

### Background

Laser scanning technology enables accurate, non-destructive 3D forest structure measurement for forest resource inventory and monitoring. Precise multi-view and multi-scan point cloud registration is a fundamental prerequisite for unified forest scene reconstruction and reliable forest parameter interpretation. The complete point cloud registration pipeline consists of two core stages:

1. **Coarse registration**: Roughly calibrates the relative rotation and translation offsets between different point cloud scans to eliminate large pose deviations;

2. **Fine registration**: Further optimizes and refines alignment accuracy to achieve high-precision fusion of overlapping point cloud regions.

Traditional fine registration methods (ICP, GICP, NDT, etc.) are originally designed for structured artificial scenes. When applied to complex forest environments, they suffer from poor robustness and low precision due to repetitive tree trunk geometric features, severe occlusion, sensor noise, and limited overlapping areas between point clouds.

### Core Idea of CC-ICP

Targeting the inherent geometric characteristics of forest scenes, CC-ICP constructs dual constraint mechanisms to optimize the traditional ICP objective function, effectively solving the poor adaptation problem of classic algorithms in forest point cloud registration:

1. **Ground plane constraint (vertical stabilization)**: Extract high-precision ground points, construct point-to-plane residual terms, and provide stable vertical direction constraint for registration;

2. **Tree cylinder constraint (horizontal stabilization)**: Perform adaptive height slicing on point clouds based on ground interpolation, extract tree trunk primitives through DBSCAN clustering combined with RANSAC cylinder fitting, and build point-to-point residual terms relying on trunk cylinder center features to constrain horizontal offset errors.

The algorithm adaptively fuses two types of residual information through a weighted least squares framework, iteratively solves the rigid transformation matrix via the Gauss-Newton method, and dynamically updates ground features and trunk cylinder primitives in each iteration to realize robust and high-precision fine alignment of forest point clouds.

### Key Pipeline

1. Ground point extraction from raw forest point clouds;

2. Point cloud slicing;

3. Tree trunk detection: DBSCAN spatial clustering + RANSAC cylinder fitting;

4. Multi-feature matching: fusion of ground plane correspondences and trunk cylinder correspondences;

5. Weighted joint least-squares optimization to solve rigid transformation parameters;

6. Iterative transformation update until convergence to obtain the optimal registration result.


## Dependencies

This implementation depends on classic point cloud processing and linear algebra libraries, with good compatibility and portability:

- PCL >= 1.8 (Tested on PCL 1.11.1)

- Eigen >= 3.3

- C++14 and above compiler

## Usage

### Input Requirement

- Input: Two PCD format forest point cloud files

- `source.pcd`: Source point cloud to be transformed and aligned

- `target.pcd`: Reference target point cloud

- Output: Standard 4×4 rigid transformation matrix 

### Adjustable Core Parameters

All parameters are centrally defined in the `main()` function, which can be flexibly tuned according to the expected accuracy of coarse registration and the quality of the point cloud data,:

|Parameter|Description|
|---|---|
|slice_height|Vertical slice thickness for trunk segmentation (m)|
|max_iterations|Maximum iterations of ICP optimization|
|tolerance|Convergence threshold for transformation update|
|corr_dist|Maximum spatial distance threshold for trunk feature matching|
|threshold_radius|Allowable radius difference for matched trunk cylinders|
|threshold_axisdirection|Maximum axis angle difference (°) for trunk matching|
|num_cylinder|Minimum valid point count for a single trunk cylinder|
|resolution|Grid size for ground point extraction|
|num_groundPnt|Number of nearest ground points for IDW height interpolation|

## Citation

If you use this code, algorithm ideas or experimental results in your academic research, please cite our work.

## Contact

For code bugs, parameter tuning guidance, academic communication and cooperation, please contact:

Xusong Dai <daixusong@whu.edu.cn>
