#pragma once
#include "igl/opengl/glfw/Viewer.h"
#include "igl/aabb.h"

class SandBox : public igl::opengl::glfw::Viewer
{
public:
	bool draw_inter = false;
	Eigen::AlignedBox<double, 3> box_1_to_remove;
	Eigen::AlignedBox<double, 3> box_2_to_remove;

	SandBox();
	~SandBox();
	void Init(const std::string& config);
	double doubleVariable;
	bool check_collision_cond(Eigen::AlignedBox<double, 3> box_1, Eigen::AlignedBox<double, 3> box_2, Eigen::Matrix3d &R_mat, Eigen::Matrix4d& A_trans, Eigen::Matrix4d& B_trans, Eigen::Matrix3d& A_rot, Eigen::Matrix3d& B_rot);
	bool check_collision_rec(igl::AABB<Eigen::MatrixXd, 3>* node1, igl::AABB<Eigen::MatrixXd, 3>* node2, Eigen::Matrix3d &R_mat, Eigen::Matrix4d& A_trans, Eigen::Matrix4d& B_trans, Eigen::Matrix3d& A_rot, Eigen::Matrix3d& B_rot);
	void check_collision();
		

private:
	// Prepare array-based edge data structures and priority queue
	
	
	void Animate();
};

