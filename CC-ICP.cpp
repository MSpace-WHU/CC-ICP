/**
 * CC-ICP: Cylinder-Constrained ICP for Forest Point Cloud Registration
 *
 * This implementation combines ground-plane constraints and cylinder-axis constraints
 * for fine registration of forest point clouds. It extracts ground points via grid-minimum
 * filtering, detects tree trunks as cylinders using DBSCAN + RANSAC, and iteratively
 * optimizes the transformation using a weighted least-squares formulation.
 *
 * Dependencies: PCL 1.8+, Eigen 3.3+
 *
 * Usage: This is a demonstration version with hard-coded file paths.
 *        Modify file names and parameters as needed.
 */

#include <iostream>
#include <vector>
#include <unordered_map>
#include <cmath>

#include <pcl/io/pcd_io.h>
#include <pcl/point_types.h>
#include <pcl/filters/voxel_grid.h>
#include <pcl/filters/grid_minimum.h>
#include <pcl/filters/extract_indices.h>
#include <pcl/features/normal_3d.h>
#include <pcl/segmentation/sac_segmentation.h>
#include <pcl/sample_consensus/method_types.h>
#include <pcl/sample_consensus/model_types.h>
#include <pcl/ModelCoefficients.h>
#include <pcl/console/time.h>
#include <pcl/common/transforms.h>      
#include <pcl/search/kdtree.h>          

#include <Eigen/Dense>

using namespace pcl;
using namespace std;

typedef PointXYZ PointT;
typedef PointCloud<PointT> PointCloudT;

// ----------------------------------------------------------------------------
// Data structures
// ----------------------------------------------------------------------------
struct ClusterResult {
    PointCloudT::Ptr cluster_points;   // points belonging to the cylinder
    float radius;                      // cylinder radius
    Eigen::Vector3f axis;              // cylinder axis direction (unit vector)
    Eigen::Vector3f center;            // centroid of the cluster
};

// ----------------------------------------------------------------------------
// Helper functions
// ----------------------------------------------------------------------------

/**
 * Read a PCD point cloud from file.
 */
void readPointCloud(const string& filename, PointCloudT::Ptr& cloud) {
    if (io::loadPCDFile<PointT>(filename, *cloud) == -1) {
        PCL_ERROR("Couldn't read file %s \n", filename.c_str());
        return;
    }
}

/**
 * Compute height above ground for a given point using IDW interpolation from k nearest ground points.
 */
float computeHeightAboveGroundKDTree(const PointXYZ& point,
    KdTreeFLANN<PointXYZ>& kdtree,
    const PointCloudT::Ptr& ground_points,
    int k) {
    vector<int> nearest_indices(k);
    vector<float> nearest_distances(k);

    if (kdtree.nearestKSearch(point, k, nearest_indices, nearest_distances) > 0) {
        float total_weight = 0.0f;
        float weighted_height = 0.0f;

        for (int i = 0; i < k; ++i) {
            float horizontal_dist = sqrt(
                pow(point.x - ground_points->points[nearest_indices[i]].x, 2) +
                pow(point.y - ground_points->points[nearest_indices[i]].y, 2));

            float weight = 1.0f / (horizontal_dist + 1e-6f);
            weighted_height += weight * ground_points->points[nearest_indices[i]].z;
            total_weight += weight;
            if (horizontal_dist == 0.0) {
                weighted_height = 0.0;
                break;
            }
        }
        return (point.z - (weighted_height / total_weight));
    }
    else {
        return numeric_limits<float>::max();
    }
}

/**
 * Slice the point cloud based on height above ground using KD-tree interpolation.
 */
void slicePointCloudBasedOnGroundKDTree(const PointCloudT::Ptr& cloud,
    const PointCloudT::Ptr& ground_points,
    float slice_height,
    vector<PointCloudT::Ptr>& slices,
    int k) {
    KdTreeFLANN<PointXYZ> kdtree;
    kdtree.setInputCloud(ground_points);

    unordered_map<int, PointCloudT::Ptr> slice_map;

    for (const auto& point : cloud->points) {
        float height = computeHeightAboveGroundKDTree(point, kdtree, ground_points, k);
        if (height < 0) continue;

        int slice_idx = static_cast<int>(height / slice_height);
        if (slice_map.find(slice_idx) == slice_map.end()) {
            slice_map[slice_idx] = PointCloudT::Ptr(new PointCloudT);
        }
        slice_map[slice_idx]->points.push_back(point);
    }

    for (auto& pair : slice_map) {
        slices.push_back(pair.second);
    }
}

/**
 * Density-Based Spatial Clustering of Applications with Noise (DBSCAN).
 * Returns true if at least one cluster is found.
 */
template<typename PointT>
bool DBSCAN(const PointCloud<PointT>& cloud_in,
    vector<Indices>& cluster_idx,
    const double& epsilon,
    const int& minpts) {
    if (cloud_in.empty()) return false;

    typename search::KdTree<PointT>::Ptr tree(new search::KdTree<PointT>());
    typename PointCloud<PointT>::ConstPtr cloud_ptr = cloud_in.makeShared();
    tree->setInputCloud(cloud_ptr);

    vector<bool> cloud_processed(cloud_in.size(), false);

    for (size_t i = 0; i < cloud_in.size(); ++i) {
        if (cloud_processed[i]) continue;

        Indices seed_queue;
        Indices k_indices;
        vector<float> k_distances;

        if (tree->radiusSearch(cloud_in.points[i], epsilon, k_indices, k_distances) >= minpts) {
            seed_queue.push_back(i);
            cloud_processed[i] = true;
        }
        else {
            continue;
        }

        int seed_index = 0;
        while (seed_index < seed_queue.size()) {
            Indices indices;
            vector<float> dists;
            if (tree->radiusSearch(cloud_in.points[seed_queue[seed_index]], epsilon, indices, dists) < minpts) {
                ++seed_index;
                continue;
            }
            for (size_t j = 0; j < indices.size(); ++j) {
                if (cloud_processed[indices[j]]) continue;
                seed_queue.push_back(indices[j]);
                cloud_processed[indices[j]] = true;
            }
            ++seed_index;
        }
        cluster_idx.push_back(seed_queue);
    }

    return !cluster_idx.empty();
}

/**
 * Fit a cylinder to a point cluster using RANSAC from normals.
 * Returns true if a valid cylinder is found and the cluster size is sufficient.
 */
bool checkCylinder(PointCloud<PointXYZ>::Ptr cloud,
    Indices idx,
    PointCloud<PointXYZ>::Ptr cluster,
    Eigen::Vector3f& axis,
    Eigen::Vector3f& center,
    float& radius,
    int num_cylinder) {
    if (idx.empty()) return false;

    PointCloud<PointXYZ>::Ptr selected_cloud(new PointCloud<PointXYZ>);
    PointCloud<Normal>::Ptr selected_normals(new PointCloud<Normal>);
    copyPointCloud(*cloud, idx, *selected_cloud);

    NormalEstimation<PointT, Normal> ne;
    ne.setInputCloud(selected_cloud);
    ne.setKSearch(8);
    ne.compute(*selected_normals);

    SACSegmentationFromNormals<PointXYZ, Normal> seg;
    ModelCoefficients::Ptr coefficients(new ModelCoefficients);
    PointIndices::Ptr inliers(new PointIndices);

    seg.setOptimizeCoefficients(true);
    seg.setModelType(SACMODEL_CYLINDER);
    seg.setMethodType(SAC_RANSAC);
    seg.setNormalDistanceWeight(0.1);
    seg.setMaxIterations(500);
    seg.setDistanceThreshold(0.05);
    seg.setRadiusLimits(0.02, 0.5);
    seg.setInputCloud(selected_cloud);
    seg.setInputNormals(selected_normals);
    seg.segment(*inliers, *coefficients);

    if (inliers->indices.empty()) return false;

    ExtractIndices<PointXYZ> extract;
    extract.setInputCloud(selected_cloud);
    extract.setIndices(inliers);
    extract.setNegative(false);
    extract.filter(*cluster);

    if (cluster->points.size() < num_cylinder) return false;

    // Cylinder axis
    axis = Eigen::Vector3f(coefficients->values[3],
        coefficients->values[4],
        coefficients->values[5]);
    if (coefficients->values[5] < 0) axis *= -1;

    radius = coefficients->values[6];

    // Compute centroid of the cylinder cluster as its center
    center = Eigen::Vector3f::Zero();
    for (const auto& pt : cluster->points) {
        center += Eigen::Vector3f(pt.x, pt.y, pt.z);
    }
    center /= static_cast<float>(cluster->points.size());

    return true;
}

/**
 * Find correspondences between two sets of cylinder clusters based on spatial proximity,
 * radius similarity and axis direction similarity.
 */
void findCorrespondences(const vector<ClusterResult>& target_clusters,
    const vector<ClusterResult>& source_clusters,
    vector<pair<int, int>>& correspondences,
    float corr_dist,
    float corr_radius,
    float corr_angle) {
    search::KdTree<PointT>::Ptr kdtree(new search::KdTree<PointT>);

    PointCloudT::Ptr target_cloud(new PointCloudT);
    for (const auto& result : target_clusters) {
        PointXYZ p;
        p.x = result.center[0];
        p.y = result.center[1];
        p.z = result.center[2];
        target_cloud->points.push_back(p);
    }
    kdtree->setInputCloud(target_cloud);

    for (size_t i = 0; i < source_clusters.size(); ++i) {
        PointT source_center;
        source_center.x = source_clusters[i].center[0];
        source_center.y = source_clusters[i].center[1];
        source_center.z = source_clusters[i].center[2];

        vector<int> nearest_indices;
        vector<float> distances;
        if (kdtree->nearestKSearch(source_center, 1, nearest_indices, distances) > 0) {
            if (!nearest_indices.empty() && distances[0] <= corr_dist * corr_dist) {
                int tidx = nearest_indices[0];
                float dr = fabs(source_clusters[i].radius - target_clusters[tidx].radius);
                float angle = acos(source_clusters[i].axis.dot(target_clusters[tidx].axis)) * 180 / M_PI;
                if (dr <= corr_radius && angle <= corr_angle) {
                    correspondences.emplace_back(tidx, static_cast<int>(i));
                }
            }
        }
    }
}

/**
 * Compute rotation matrix that aligns a given normal to the Z-axis (0,0,1).
 */
Eigen::Affine3f computeRotationMatrix(Eigen::Vector3f normal) {
    normal.normalize();
    Eigen::Vector3f target(0, 0, 1);
    Eigen::Vector3f axis = normal.cross(target);
    float angle = acos(normal.dot(target));

    if (axis.norm() < 1e-6) {
        return Eigen::Affine3f::Identity();
    }

    axis.normalize();
    Eigen::Matrix3f K;
    K << 0, -axis.z(), axis.y(),
        axis.z(), 0, -axis.x(),
        -axis.y(), axis.x(), 0;

    Eigen::Matrix3f R = Eigen::Matrix3f::Identity() + sin(angle) * K + (1 - cos(angle)) * K * K;
    Eigen::Affine3f affine = Eigen::Affine3f::Identity();
    affine.linear() = R;
    return affine;
}

/**
 * Least-squares circle fitting in 2D (X-Y plane). Returns center (x,y,z) and radius.
 */
bool fitCircleLeastSquares(const PointCloud<PointXYZ>::Ptr& cloud,
    Eigen::Vector3f& center,
    float& radius) {
    if (cloud->points.size() < 3) {
        cerr << "Error: At least 3 points required to fit a circle." << endl;
        return false;
    }

    Eigen::Vector3f centroid(0, 0, 0);
    for (const auto& pt : cloud->points) {
        centroid += Eigen::Vector3f(pt.x, pt.y, pt.z);
    }
    centroid /= cloud->points.size();

    Eigen::MatrixXf A(cloud->points.size(), 3);
    Eigen::VectorXf b(cloud->points.size());

    for (size_t i = 0; i < cloud->points.size(); ++i) {
        float x = cloud->points[i].x - centroid.x();
        float y = cloud->points[i].y - centroid.y();
        A(i, 0) = x;
        A(i, 1) = y;
        A(i, 2) = 1.0;
        b(i) = x * x + y * y;
    }

    Eigen::Vector3f result = (A.transpose() * A).ldlt().solve(A.transpose() * b);

    float cx = result(0) / 2.0 + centroid.x();
    float cy = result(1) / 2.0 + centroid.y();
    float r = sqrt(result(2) + pow(result(0) / 2.0, 2) + pow(result(1) / 2.0, 2));

    center[0] = cx;
    center[1] = cy;
    center[2] = centroid.z();
    radius = r;
    return true;
}

/**
 * Least-squares circle fitting with a fixed radius. Returns only the center.
 */
bool fitCircleWithFixedRadius(const PointCloud<PointXYZ>::Ptr& cloud,
    float fixed_radius,
    Eigen::Vector3f& center) {
    if (cloud->points.size() < 3) {
        cerr << "Error: At least 3 points required to fit a circle." << endl;
        return false;
    }

    Eigen::Vector3f centroid(0, 0, 0);
    for (const auto& pt : cloud->points) {
        centroid += Eigen::Vector3f(pt.x, pt.y, pt.z);
    }
    centroid /= cloud->points.size();

    Eigen::MatrixXf A(cloud->points.size(), 2);
    Eigen::VectorXf b(cloud->points.size());

    for (size_t i = 0; i < cloud->points.size(); ++i) {
        float x = cloud->points[i].x - centroid.x();
        float y = cloud->points[i].y - centroid.y();
        A(i, 0) = -2.0 * x;
        A(i, 1) = -2.0 * y;
        b(i) = fixed_radius * fixed_radius - x * x - y * y;
    }

    Eigen::Vector2f result = (A.transpose() * A).ldlt().solve(A.transpose() * b);

    center[0] = result(0) + centroid.x();
    center[1] = result(1) + centroid.y();
    center[2] = centroid.z();
    return true;
}



/**
 * Construct the weighted linear system A * x = b for joint optimization.
 * x = [alpha, beta, gamma, tx, ty, tz]^T (rotation angles and translation).
 */
void constructWeightedErrorSystem(const vector<Eigen::Vector3f>& src_points,
    const vector<Eigen::Vector3f>& tgt_points,
    const vector<string>& types,
    const vector<Eigen::Vector3f>& normals,
    Eigen::MatrixXf& A,
    Eigen::VectorXf& b) {
    vector<Eigen::VectorXf> A_rows;
    vector<float> b_values;

    size_t count_point = 0, count_plane = 0;
    for (const auto& t : types) {
        if (t == "cylinderpoint") ++count_point;
        else if (t == "plane") ++count_plane;
    }

    float w_point = count_point > 0 ? 1.0f / count_point : 1.0f;
    float w_plane = count_plane > 0 ? 1.0f / count_plane : 1.0f;

    for (size_t i = 0; i < src_points.size(); ++i) {
        const Eigen::Vector3f& p = src_points[i];
        const Eigen::Vector3f& q = tgt_points[i];
        const string& type = types[i];

        if (type == "plane") {
            Eigen::Vector3f n = normals[i];
            Eigen::VectorXf row(6);
            row << n[2] * p[1] - n[1] * p[2],
                n[0] * p[2] - n[2] * p[0],
                n[1] * p[0] - n[0] * p[1],
                n[0], n[1], n[2];
            float rhs = n.dot(q - p);

            row *= w_plane;
            rhs *= w_plane;

            A_rows.push_back(row);
            b_values.push_back(rhs);
        }
        else if (type == "cylinderpoint") {
            Eigen::Vector3f d = q - p;
            for (int j = 0; j < 3; ++j) {
                Eigen::VectorXf row(6);
                row.setZero();

                if (j == 0) row.head<3>() << 0, p[2], -p[1];
                if (j == 1) row.head<3>() << -p[2], 0, p[0];
                if (j == 2) row.head<3>() << p[1], -p[0], 0;

                row(3 + j) = 1;
                row *= w_point;
                float rhs = d(j) * w_point;

                A_rows.push_back(row);
                b_values.push_back(rhs);
            }
        }
    }

    A.resize(A_rows.size(), 6);
    b.resize(A_rows.size());
    for (size_t i = 0; i < A_rows.size(); ++i) {
        A.row(i) = A_rows[i];
        b(i) = b_values[i];
    }
}

// ----------------------------------------------------------------------------
// Main registration pipeline
// ----------------------------------------------------------------------------
int main() {
    // Load point clouds (hard-coded paths for demonstration)
    PointCloudT::Ptr targetcloud(new PointCloudT);
    PointCloudT::Ptr sourcecloud(new PointCloudT);
    readPointCloud("17-2.pcd", targetcloud);
    readPointCloud("17-1.pcd", sourcecloud);

    // Algorithm parameters
    float slice_height = 0.05f;
    int max_iterations = 100;
    float tolerance = 1e-5;
    float corr_dist = 0.5;          // max distance for cylinder correspondence
    float threshold_radius = 0.02;  // max radius difference for cylinder correspondence
    float threshold_axisdirection = 10.0; // max axis angle difference (degrees)
    int num_cylinder = 10;          // minimum number of points in a cylinder cluster
    float resolution = 0.05;         // grid size for ground extraction
    int num_groundPnt = 3;          // number of nearest ground points for height interpolation

    console::TicToc time;
    time.tic();

    // ------------------------------------------------------------------------
    // 1. Ground extraction using grid minimum
    // ------------------------------------------------------------------------
    PointCloud<PointXYZ>::Ptr targetground(new PointCloud<PointXYZ>);
    PointCloud<PointXYZ>::Ptr sourceground(new PointCloud<PointXYZ>);

    GridMinimum<PointXYZ> gm(resolution);
    gm.setInputCloud(sourcecloud);
    gm.filter(*sourceground);
    gm.setInputCloud(targetcloud);
    gm.filter(*targetground);


    // Compute normals for target ground points (for point-to-plane metric)
    NormalEstimation<PointXYZ, Normal> ne;
    ne.setInputCloud(targetground);
    search::KdTree<PointXYZ>::Ptr tree(new search::KdTree<PointXYZ>);
    ne.setSearchMethod(tree);
    ne.setKSearch(10);
    PointCloud<Normal>::Ptr tgtcloud_normals(new PointCloud<Normal>);
    ne.compute(*tgtcloud_normals);

    vector<Eigen::Vector3f> tgtnormals;
    tgtnormals.reserve(tgtcloud_normals->size());
    for (const auto& n : tgtcloud_normals->points) {
        if (isFinite(n)) {
            tgtnormals.emplace_back(n.normal_x, n.normal_y, n.normal_z);
        }
        else {
            tgtnormals.emplace_back(0.f, 0.f, 0.f);
        }
    }

    // ------------------------------------------------------------------------
    // 2. Cylinder detection from target cloud
    // ------------------------------------------------------------------------
    vector<PointCloudT::Ptr> target_slices;
    slicePointCloudBasedOnGroundKDTree(targetcloud, targetground, slice_height,
        target_slices, num_groundPnt);

    vector<ClusterResult> target_clusters;
    for (auto& slice : target_slices) {
        if (!slice || slice->points.empty()) continue;

        vector<Indices> cluster_indices;
        if (DBSCAN(*slice, cluster_indices, 0.1, 20)) {
            for (auto& idx : cluster_indices) {
                PointCloudT::Ptr cluster(new PointCloudT);
                Eigen::Vector3f axis, center;
                float radius;
                if (checkCylinder(slice, idx, cluster, axis, center, radius, num_cylinder)) {
                    target_clusters.push_back({ cluster, radius, axis, center });
                }
            }
        }
    }
    cout << "Time after target cylinder detection: " << time.toc() << " ms" << endl;

    // ------------------------------------------------------------------------
    // 3. Cylinder detection from source cloud
    // ------------------------------------------------------------------------
    vector<PointCloudT::Ptr> source_slices;
    slicePointCloudBasedOnGroundKDTree(sourcecloud, sourceground, slice_height,
        source_slices, num_groundPnt);

    vector<ClusterResult> source_clusters;
    for (auto& slice : source_slices) {
        if (!slice || slice->points.empty()) continue;

        vector<Indices> cluster_indices;
        if (DBSCAN(*slice, cluster_indices, 0.1, 20)) {
            for (auto& idx : cluster_indices) {
                PointCloudT::Ptr cluster(new PointCloudT);
                Eigen::Vector3f axis, center;
                float radius;
                if (checkCylinder(slice, idx, cluster, axis, center, radius, num_cylinder)) {
                    source_clusters.push_back({ cluster, radius, axis, center });
                }
            }
        }
    }
    cout << "Time after source cylinder detection: " << time.toc() << " ms" << endl;
    cout << "Target cylinders: " << target_clusters.size()
        << ", Source cylinders: " << source_clusters.size() << endl;

    // ------------------------------------------------------------------------
    // 4. Weighted ICP optimization
    // ------------------------------------------------------------------------
    Eigen::Matrix3f final_rotation = Eigen::Matrix3f::Identity();
    Eigen::Vector3f final_translation = Eigen::Vector3f::Zero();

    search::KdTree<PointXYZ>::Ptr kdtreePlane(new search::KdTree<PointXYZ>);
    kdtreePlane->setInputCloud(targetground);

    for (int iter = 0; iter < max_iterations; ++iter) {
        vector<Eigen::Vector3f> source_corr, target_corr, target_normal;
        vector<string> types;

        // 4a. Ground point correspondences (point-to-plane)
        for (size_t i = 0; i < sourceground->points.size(); ++i) {
            vector<int> idx;
            vector<float> dist;
            if (kdtreePlane->nearestKSearch(sourceground->points[i], 1, idx, dist) > 0) {
                if (dist[0] <= corr_dist * corr_dist) {
                    target_corr.emplace_back(targetground->points[idx[0]].x,
                        targetground->points[idx[0]].y,
                        targetground->points[idx[0]].z);
                    target_normal.push_back(tgtnormals[idx[0]]);
                    source_corr.emplace_back(sourceground->points[i].x,
                        sourceground->points[i].y,
                        sourceground->points[i].z);
                    types.push_back("plane");
                }
            }
        }

        // 4b. Cylinder axis correspondences (point-to-point)
        vector<pair<int, int>> correspondences;
        findCorrespondences(target_clusters, source_clusters, correspondences,
            corr_dist, threshold_radius, threshold_axisdirection);

        cout << "Cylinder correspondences: " << correspondences.size() << endl;

        for (auto& corr : correspondences) {
            int tidx = corr.first;
            int sidx = corr.second;

            Eigen::Vector3f target_axis = target_clusters[tidx].axis;
            Eigen::Vector3f source_axis = source_clusters[sidx].axis;

            Eigen::Affine3f source_transform = computeRotationMatrix(source_axis);
            Eigen::Affine3f target_transform = computeRotationMatrix(target_axis);

            PointCloudT::Ptr source_Z_project(new PointCloudT);
            PointCloudT::Ptr target_Z_project(new PointCloudT);
            transformPointCloud(*source_clusters[sidx].cluster_points, *source_Z_project, source_transform);
            transformPointCloud(*target_clusters[tidx].cluster_points, *target_Z_project, target_transform);

            // Fit circles in projected plane
            Eigen::Vector3f source_center, target_center;
            float source_radius, target_radius;
            fitCircleLeastSquares(source_Z_project, source_center, source_radius);
            fitCircleLeastSquares(target_Z_project, target_center, target_radius);

            if (fabs(target_radius - source_radius) > threshold_radius) continue;

            // Merge projected points to compute common radius
            PointCloudT::Ptr merged_project(new PointCloudT);
            PointXYZ pt;
            for (const auto& p : target_Z_project->points) {
                pt.x = p.x - target_center[0];
                pt.y = p.y - target_center[1];
                pt.z = p.z - target_center[2];
                merged_project->points.push_back(pt);
            }
            for (const auto& p : source_Z_project->points) {
                pt.x = p.x - source_center[0];
                pt.y = p.y - source_center[1];
                pt.z = p.z - source_center[2];
                merged_project->points.push_back(pt);
            }

            Eigen::Vector3f merged_center;
            float merged_radius;
            fitCircleLeastSquares(merged_project, merged_center, merged_radius);

            // Refit circles with fixed radius to get precise centers
            Eigen::Vector3f target_xoy_center, source_xoy_center;
            fitCircleWithFixedRadius(target_Z_project, merged_radius, target_xoy_center);
            fitCircleWithFixedRadius(source_Z_project, merged_radius, source_xoy_center);

            // Back-project to original 3D coordinates
            Eigen::Vector3f final_target_center, final_source_center;
            transformPoint(target_xoy_center, final_target_center, target_transform.inverse());
            transformPoint(source_xoy_center, final_source_center, source_transform.inverse());

            target_corr.push_back(final_target_center);
            source_corr.push_back(final_source_center);
            target_normal.emplace_back(0.f, 0.f, 0.f);
            types.push_back("cylinderpoint");
        }

        if (source_corr.size() < 20) {
            cout << "Too few correspondences, skipping iteration." << endl;
            continue;
        }

        // 4c. Solve weighted least squares
        Eigen::MatrixXf A;
        Eigen::VectorXf b;
        constructWeightedErrorSystem(source_corr, target_corr, types, target_normal, A, b);

        Eigen::VectorXf x = A.colPivHouseholderQr().solve(b);

        float alpha = x(0), beta = x(1), gamma = x(2);
        Eigen::Vector3f translation = x.segment<3>(3);

        Eigen::Matrix3f Rx = Eigen::AngleAxisf(alpha, Eigen::Vector3f::UnitX()).toRotationMatrix();
        Eigen::Matrix3f Ry = Eigen::AngleAxisf(beta, Eigen::Vector3f::UnitY()).toRotationMatrix();
        Eigen::Matrix3f Rz = Eigen::AngleAxisf(gamma, Eigen::Vector3f::UnitZ()).toRotationMatrix();
        Eigen::Matrix3f rotation = Rz * Ry * Rx;

        // 4d. Update source point clouds and cylinder parameters
        Eigen::Affine3f iter_transform = Eigen::Affine3f::Identity();
        iter_transform.linear() = rotation;
        iter_transform.translation() = translation;

        transformPointCloud(*sourceground, *sourceground, iter_transform);
        for (auto& cluster : source_clusters) {
            cluster.center = rotation * cluster.center + translation;
            cluster.axis = rotation * cluster.axis;
            transformPointCloud(*cluster.cluster_points, *cluster.cluster_points, iter_transform);
        }

        // 4e. Accumulate global transformation
        final_rotation = rotation * final_rotation;
        final_translation = rotation * final_translation + translation;

        cout << "Iteration " << iter + 1 << endl;
        cout << "  Rotation:\n" << final_rotation << endl;
        cout << "  Translation: " << final_translation.transpose() << endl;

        // Check convergence
        float change = translation.norm() + (rotation - Eigen::Matrix3f::Identity()).norm();
        if (change < tolerance) {
            cout << "Converged after " << iter + 1 << " iterations." << endl;
            break;
        }
    }

    // ------------------------------------------------------------------------
    // 5. Output final transformation
    // ------------------------------------------------------------------------
    cout << "Total time: " << time.toc() << " ms" << endl;
    Eigen::Affine3f final_transform = Eigen::Affine3f::Identity();
    final_transform.linear() = final_rotation;
    final_transform.translation() = final_translation;

    cout << "Final Rotation Matrix:\n" << final_rotation << endl;
    cout << "Final Translation Vector:\n" << final_translation.transpose() << endl;
    cout << "Final Transformation Matrix:\n" << final_transform.matrix() << endl;

    return 0;
}