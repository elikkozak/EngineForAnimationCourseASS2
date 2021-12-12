#include "tutorial/sandBox/sandBox.h"
#include "igl/edge_flaps.h"
#include "igl/collapse_edge.h"
#include "Eigen/dense"
#include <functional>



SandBox::SandBox()
{
	

}

bool SandBox::check_collision_cond(Eigen::AlignedBox<double, 3> box_1, Eigen::AlignedBox<double, 3> box_2, Eigen::Matrix3d &R_mat, Eigen::Matrix4d& A_trans, Eigen::Matrix4d& B_trans, Eigen::Matrix3d& A_rot, Eigen::Matrix3d& B_rot)
{
	double R, R_0, R_1;

	Eigen::Vector4d A_vec_center = Eigen::Vector4d(box_1.center()[0], box_1.center()[1], box_1.center()[2], 1);
	Eigen::Vector4d A_vec_center_trans = A_trans * A_vec_center;
	Eigen::RowVector3d P_A = Eigen::RowVector3d(A_vec_center_trans.x(), A_vec_center_trans.y(), A_vec_center_trans.z());

	Eigen::RowVector3d A_x = A_rot * Eigen::Vector3d(1, 0, 0);
	Eigen::RowVector3d A_y = A_rot * Eigen::Vector3d(0, 1, 0);
	Eigen::RowVector3d A_z = A_rot * Eigen::Vector3d(0, 0, 1);
	double W_A = box_1.sizes()[0] / 2;
	double H_A = box_1.sizes()[1] / 2;
	double D_A = box_1.sizes()[2] / 2;

	

	Eigen::Vector4d B_vec_center = Eigen::Vector4d(box_2.center()[0], box_2.center()[1], box_2.center()[2], 1);
	Eigen::Vector4d B_vec_center_trans = B_trans * B_vec_center;
	Eigen::RowVector3d P_B = Eigen::RowVector3d(B_vec_center_trans.x(), B_vec_center_trans.y(), B_vec_center_trans.z());

	Eigen::RowVector3d B_x = B_rot * Eigen::Vector3d(1, 0, 0);
	Eigen::RowVector3d B_y = B_rot * Eigen::Vector3d(0, 1, 0);
	Eigen::RowVector3d B_z = B_rot * Eigen::Vector3d(0, 0, 1);
	double W_B = box_2.sizes()[0] / 2;
	double H_B = box_2.sizes()[1] / 2;
	double D_B = box_2.sizes()[2] / 2;



	Eigen::RowVector3d T = P_B - P_A;
	//Eigen::Matrix3d R_mat = A * B;

	//std::cout << T << std::endl;

	//CASE 1:
	R_0 = W_A;
	R_1 = W_B * abs(R_mat(0, 0)) + H_B * abs(R_mat(0, 1)) + D_B * abs(R_mat(0, 2));
	R = abs(T.dot(A_x));
	if (R > R_0 + R_1)
		return false;

	//CASE 2:
	R_0 = H_A;
	R_1 = W_B * abs(R_mat(1, 0)) + H_B * abs(R_mat(1, 1)) + D_B * abs(R_mat(1, 2));
	R = abs(T.dot(A_y));
	if (R > R_0 + R_1)
		return false;

	//CASE 3:
	R_0 = D_A;
	R_1 = W_B * abs(R_mat(2, 0)) + H_B * abs(R_mat(2, 1)) + D_B * abs(R_mat(2, 2));
	R = abs(T.dot(A_z));
	if (R > R_0 + R_1)
		return false;

	//CASE 4:
	R_0 = W_A * abs(R_mat(0, 0)) + H_A * abs(R_mat(1, 0)) + D_A * abs(R_mat(2, 0));
	R_1 = W_B;
	R = abs(T.dot(B_x));
	if (R > R_0 + R_1)
		return false;

	//CASE 5:
	R_0 = W_A * abs(R_mat(0, 1)) + H_A * abs(R_mat(1, 1)) + D_A * abs(R_mat(2, 1));
	R_1 = H_B;
	R = abs(T.dot(B_y));
	if (R > R_0 + R_1)
		return false;

	//CASE 6:
	R_0 = W_A * abs(R_mat(0, 2)) + H_A * abs(R_mat(1, 2)) + D_A * abs(R_mat(2, 2));
	R_1 = D_B;
	R = abs(T.dot(B_z));
	if (R > R_0 + R_1)
		return false;

	//CASE 7:
	R_0 = H_A * abs(R_mat(2, 0)) + D_A * abs(R_mat(1, 0));
	R_1 = H_B * abs(R_mat(0, 2)) + D_B * abs(R_mat(0, 1));
	R = abs(R_mat(1, 0) * T.dot(A_z) - R_mat(2, 0) * T.dot(A_y));
	if (R > R_0 + R_1)
		return false;

	//CASE 8:
	R_0 = H_A * abs(R_mat(2, 1)) + D_A * abs(R_mat(1, 1));
	R_1 = W_B * abs(R_mat(0, 2)) + D_B * abs(R_mat(0, 0));
	R = abs(R_mat(1, 1) * T.dot(A_z) - R_mat(2, 1) * T.dot(A_y));
	if (R > R_0 + R_1)
		return false;

	//CASE 9:
	R_0 = H_A * abs(R_mat(2, 2)) + D_A * abs(R_mat(1, 2));
	R_1 = W_B * abs(R_mat(0, 1)) + H_B * abs(R_mat(0, 0));
	R = abs(R_mat(1, 2) * T.dot(A_z) - R_mat(2, 2) * T.dot(A_y));
	if (R > R_0 + R_1)
		return false;

	//CASE 10:
	R_0 = W_A * abs(R_mat(2, 0)) + D_A * abs(R_mat(0, 0));
	R_1 = H_B * abs(R_mat(1, 2)) + D_B * abs(R_mat(1, 1));
	R = abs(R_mat(2, 0) * T.dot(A_x) - R_mat(0, 0) * T.dot(A_z));
	if (R > R_0 + R_1)
		return false;

	//CASE 11:
	R_0 = W_A * abs(R_mat(2, 1)) + D_A * abs(R_mat(0, 1));
	R_1 = W_B * abs(R_mat(1, 2)) + D_B * abs(R_mat(1, 0));
	R = abs(R_mat(2, 1) * T.dot(A_x) - R_mat(0, 1) * T.dot(A_z));
	if (R > R_0 + R_1)
		return false;

	//CASE 12:
	R_0 = W_A * abs(R_mat(2, 2)) + D_A * abs(R_mat(0, 2));
	R_1 = W_B * abs(R_mat(1, 1)) + H_B * abs(R_mat(1, 0));
	R = abs(R_mat(2, 2) * T.dot(A_x) - R_mat(0, 2) * T.dot(A_z));
	if (R > R_0 + R_1)
		return false;

	//CASE 13:
	R_0 = W_A * abs(R_mat(1, 0)) + H_A * abs(R_mat(0, 0));
	R_1 = H_B * abs(R_mat(2, 2)) + D_B * abs(R_mat(2, 1));
	R = abs(R_mat(0, 0) * T.dot(A_y) - R_mat(1, 0) * T.dot(A_x));
	if (R > R_0 + R_1)
		return false;

	//CASE 14:
	R_0 = W_A * abs(R_mat(1, 1)) + H_A * abs(R_mat(0, 1));
	R_1 = W_B * abs(R_mat(2, 2)) + D_B * abs(R_mat(2, 0));
	R = abs(R_mat(0, 1) * T.dot(A_y) - R_mat(1, 1) * T.dot(A_x));
	if (R > R_0 + R_1)
		return false;

	//CASE 15:
	R_0 = W_A * abs(R_mat(1, 2)) + H_A * abs(R_mat(0, 2));
	R_1 = W_B * abs(R_mat(2, 1)) + H_B * abs(R_mat(2, 0));
	R = abs(R_mat(0, 2) * T.dot(A_y) - R_mat(1, 2) * T.dot(A_x));
	if (R > R_0 + R_1)
		return false;

	return true;
}

bool SandBox::check_collision_rec(igl::AABB<Eigen::MatrixXd, 3>* node1, igl::AABB<Eigen::MatrixXd, 3>* node2, Eigen::Matrix3d &R_mat, Eigen::Matrix4d& A_trans, Eigen::Matrix4d& B_trans, Eigen::Matrix3d& A_rot, Eigen::Matrix3d& B_rot)
{
	if(check_collision_cond(node1->m_box,node2->m_box, R_mat,  A_trans,  B_trans,  A_rot,  B_rot))
	{
		if (node1->is_leaf() && node2->is_leaf())
		{
			draw_inter = true;
			data_list[0].drawBox(node1->m_box);
			data_list[1].drawBox(node2->m_box);

			box_1_to_remove = node1->m_box;
			box_2_to_remove = node2->m_box;
			
			dir = -1 * dir;
			return true;
		}
		
		igl::AABB<Eigen::MatrixXd, 3>* left1 = node1->is_leaf() ? node1 : node1->m_left;
		igl::AABB<Eigen::MatrixXd, 3>* right1 = node1->is_leaf() ? node1 : node1->m_right;
		igl::AABB<Eigen::MatrixXd, 3>* left2 = node2->is_leaf() ? node2 : node2->m_left;;
		igl::AABB<Eigen::MatrixXd, 3>* right2 = node2->is_leaf() ? node2 : node2->m_right;

		if (check_collision_rec(left1, left2, R_mat, A_trans, B_trans, A_rot, B_rot))
			return true;
		else if (check_collision_rec(left1, right2, R_mat, A_trans, B_trans, A_rot, B_rot))
			return true;
		else if (check_collision_rec(right1, left2, R_mat, A_trans, B_trans, A_rot, B_rot))
			return true;
		else if (check_collision_rec(right1, right2, R_mat, A_trans, B_trans, A_rot, B_rot))
			return true;
		else return false;
		
	}
	return false;
}

void SandBox::check_collision()
{
	igl::AABB<Eigen::MatrixXd, 3> *node1 = &data_list[0].tree;
	igl::AABB<Eigen::MatrixXd, 3> *node2 = &data_list[1].tree;

	Eigen::Matrix4d A_trans = data_list[0].MakeTransd();
	Eigen::Matrix3d A_rot = data_list[0].GetRotation();

	Eigen::Matrix4d B_trans = data_list[1].MakeTransd();
	Eigen::Matrix3d B_rot = data_list[1].GetRotation();

	Eigen::RowVector3d A_x = A_rot * Eigen::Vector3d(1, 0, 0);
	Eigen::RowVector3d A_y = A_rot * Eigen::Vector3d(0, 1, 0);
	Eigen::RowVector3d A_z = A_rot * Eigen::Vector3d(0, 0, 1);


	Eigen::Matrix3d A;
	A << A_x[0], A_x[1], A_x[2],
		A_y[0], A_y[1], A_y[2],
		A_z[0], A_z[1], A_z[2];


	Eigen::RowVector3d B_x = B_rot * Eigen::Vector3d(1, 0, 0);
	Eigen::RowVector3d B_y = B_rot * Eigen::Vector3d(0, 1, 0);
	Eigen::RowVector3d B_z = B_rot * Eigen::Vector3d(0, 0, 1);

	Eigen::Matrix3d B;
	B << B_x[0], B_y[0], B_z[0],
		B_x[1], B_y[1], B_z[1],
		B_x[2], B_y[2], B_z[2];

	Eigen::Matrix3d R_mat = A * B;

	



	if (check_collision_rec(node1, node2, R_mat, A_trans, B_trans, A_rot, B_rot)) {
		std::cout << "Collision detected! " << std::endl;
		isMoving = false;
	}
}

void SandBox::Init(const std::string &config)
{
	std::string item_name;
	std::ifstream nameFileout;
	double i = -0.5;
	doubleVariable = 0;
	nameFileout.open(config);
	if (!nameFileout.is_open())
	{
		std::cout << "Can't open file "<<config << std::endl;
	}
	else
	{
		
		while (nameFileout >> item_name)
		{
			std::cout << "openning " << item_name << std::endl;
			load_mesh_from_file(item_name);
			
			parents.push_back(-1);
			data().add_points(Eigen::RowVector3d(0, 0, 0), Eigen::RowVector3d(0, 0, 1));
			//data().show_overlay = true;
			//data().show_lines = true;
			data().show_overlay_depth = true;
			data().point_size = 10;
			data().line_width = 2;
			data().set_visible(false, 1);

			//SET LOCATION NO TOUCHY
			data().MyTranslate( Eigen::Vector3d(i*2, 0, 0),true);
			i++;

			//DRAW A BOX
			data().tree.init(data().V, data().F);
			igl::AABB<Eigen::MatrixXd, 3> tree = data().tree;
			Eigen::AlignedBox<double, 3> box = tree.m_box;
			data().drawBox(box);
		}
		nameFileout.close();
	}
	MyTranslate(Eigen::Vector3d(0, 0, -1), true);
	
	data().set_colors(Eigen::RowVector3d(0.9, 0.1, 0.1));

}

SandBox::~SandBox()
{

}

void SandBox::Animate()
{
	if (isActive)
	{

	}

	if(draw_inter && isMoving)
	{
		data_list[0].removeDrawBox(box_1_to_remove);
		data_list[1].removeDrawBox(box_2_to_remove);
		data_list[0].drawBox(data_list[0].tree.m_box);
		data_list[1].drawBox(data_list[1].tree.m_box);
		//data_list[1].removeDrawBox(box_2_to_remove);
		draw_inter = false;
	}
	
	if (isMoving) {
		data_list[0].TranslateInVelocity(GetRotation(),dir);
		check_collision();
		
	}
	
		
		
}


