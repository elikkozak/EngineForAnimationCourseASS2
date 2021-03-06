#include "igl/opengl/glfw/renderer.h"

#include <GLFW/glfw3.h>
#include <igl/unproject_onto_mesh.h>
#include "igl/look_at.h"
//#include <Eigen/Dense>

Renderer::Renderer() : selected_core_index(0),
next_core_id(2)
{
	core_list.emplace_back(igl::opengl::ViewerCore());
	core_list.front().id = 1;
	// C-style callbacks
	callback_init = nullptr;
	callback_pre_draw = nullptr;
	callback_post_draw = nullptr;
	callback_mouse_down = nullptr;
	callback_mouse_up = nullptr;
	callback_mouse_move = nullptr;
	callback_mouse_scroll = nullptr;
	callback_key_down = nullptr;
	callback_key_up = nullptr;

	callback_init_data = nullptr;
	callback_pre_draw_data = nullptr;
	callback_post_draw_data = nullptr;
	callback_mouse_down_data = nullptr;
	callback_mouse_up_data = nullptr;
	callback_mouse_move_data = nullptr;
	callback_mouse_scroll_data = nullptr;
	callback_key_down_data = nullptr;
	callback_key_up_data = nullptr;
	//highdpi = 1;

	xold = 0;
	yold = 0;
}

IGL_INLINE void Renderer::draw( GLFWwindow* window)
{
	using namespace std;
	using namespace Eigen;

	int width, height;
	glfwGetFramebufferSize(window, &width, &height);

	int width_window, height_window;
	glfwGetWindowSize(window, &width_window, &height_window);
	
	auto highdpi_tmp = (width_window == 0 || width == 0) ? highdpi : (width / width_window);

	if (fabs(highdpi_tmp - highdpi) > 1e-8)
	{
		post_resize(window,width, height);
		highdpi = highdpi_tmp;
	}

	for (auto& core : core_list)
	{
		core.clear_framebuffers();
	}
	int coreIndx = 1;
	if (menu)
	{
		menu->pre_draw();
		menu->callback_draw_viewer_menu();
	}
	for (auto& core : core_list)
	{
		int indx = 0;
		for (auto& mesh : scn->data_list)
		{
			
			if (mesh.is_visible & core.id)
			{// for kinematic chain change scn->MakeTrans to parent matrix
				
				core.draw(scn->MakeTransScale()*scn->CalcParentsTrans(indx).cast<float>(),mesh);
			}
			indx++;
		}

		
	}
	if (menu)
	{
		menu->post_draw();
		
	}

}

void Renderer::SetScene(igl::opengl::glfw::Viewer* viewer)
{
	scn = viewer;
}

IGL_INLINE void Renderer::init(igl::opengl::glfw::Viewer* viewer,int coresNum, igl::opengl::glfw::imgui::ImGuiMenu* _menu)
{
	scn = viewer;
	
	doubleVariable = 0;
	core().init(); 
	menu = _menu;
	core().align_camera_center(scn->data().V, scn->data().F);

	if (coresNum > 1)
	{	
		int width = core().viewport[2];
		int height = core().viewport[3];
		
		core().viewport = Eigen::Vector4f(0, 0, width/4, height);
		left_view = core_list[0].id;
		right_view = append_core(Eigen::Vector4f(width/4, 0, width*3/4, height));
		core_index(right_view - 1);
		for (size_t i = 0; i < scn->data_list.size(); i++)
		{
			core().toggle(scn->data_list[i].show_faces);
			core().toggle(scn->data_list[i].show_lines);
			core().toggle(scn->data_list[i].show_texture );
		}
		//Eigen::Vector3d v = -scn->GetCameraPosition();
		//TranslateCamera(v.cast<float>());

		core_index(left_view - 1);
	}

	if (menu)
	{
		menu->callback_draw_viewer_menu = [&]()
		{
			// Draw parent menu content
			menu->draw_viewer_menu(scn,core_list);


		};
	}
}

void Renderer::UpdatePosition(double xpos, double ypos)
{
	xrel = xold - xpos;
	yrel = yold - ypos;
	xold = xpos;
	yold = ypos;
}

void Renderer::MouseProcessing(int button)
{
	
	if (scn->isPicked )
	{
		if (button == 1)
		{
			float near = core().camera_dnear, far = core().camera_dfar, angle = core().camera_view_angle;
			//float z = far + depth * (near - far);
			
			Eigen::Matrix4f tmpM = core().proj;
			double xToMove = -(double)xrel / core().viewport[3] * (z+2*near) * (far) / (far + 2*near) * 2.0 * tanf(angle / 360 * M_PI) / (core().camera_zoom * core().camera_base_zoom);
			double yToMove = (double)yrel / core().viewport[3] *(z+2*near) * (far ) / (far+ 2*near) * 2.0 * tanf(angle / 360 * M_PI) / (core().camera_zoom * core().camera_base_zoom);
		
			scn->data().TranslateInSystem(scn->GetRotation(), Eigen::Vector3d(xToMove, 0, 0));
			scn->data().TranslateInSystem(scn->GetRotation(), Eigen::Vector3d(0, yToMove, 0));
			scn->WhenTranslate();
		}
		else
		{
			scn->data().RotateInSystem(Eigen::Vector3d(1, 0, 0), yrel / 100.0);

			scn->data().RotateInSystem(Eigen::Vector3d(0, 1, 0), xrel / 100.0);
		}
	}
	else
	{
		if (button == 1)
		{
			float near = core().camera_dnear, far = core().camera_dfar, angle = core().camera_view_angle;
			float z = far + 0.5f * (near - far);
			

			double xToMove = -(double)xrel / core().viewport[3] * far / z * near * 2.0f * tanf(angle / 360 * M_PI) / (core().camera_zoom * core().camera_base_zoom);
			double yToMove = (double)yrel / core().viewport[3] * far / z * near * 2.0f * tanf(angle / 360 * M_PI) / (core().camera_zoom * core().camera_base_zoom);
			scn->MyTranslate(Eigen::Vector3d(xToMove, 0, 0), true);
			scn->MyTranslate(Eigen::Vector3d(0, yToMove, 0), true);

		}
		else
		{
			scn->RotateInSystem(Eigen::Vector3d(1, 0, 0), yrel / 100.0);
			scn->RotateInSystem(Eigen::Vector3d(0, 1, 0), xrel / 100.0);

		}
	}
}

void Renderer::TranslateCamera(Eigen::Vector3f amt)
{
	core().camera_translation += amt;
}



void Renderer::RotateCamera(float amtX, float amtY)
{
	core().camera_eye = core().camera_eye + Eigen::Vector3f(0,amtY,0);
	Eigen::Matrix3f Mat;
		Mat << cos(amtY),0,sin(amtY),  0, 1, 0 ,  -sin(amtY), 0, cos(amtY) ;
	core().camera_eye = Mat* core().camera_eye;
	
}

Renderer::~Renderer()
{
	//if (scn)
	//	delete scn;
}

double Renderer::Picking(double newx, double newy)
{
		int fid;
		//Eigen::MatrixXd C = Eigen::MatrixXd::Constant(scn->data().F.rows(), 3, 1);
		Eigen::Vector3f bc;
		double x = newx;
		double y = core().viewport(3) - newy;
		
		Eigen::Matrix4f view = Eigen::Matrix4f::Identity();
		
		igl::look_at(core().camera_eye, core().camera_center, core().camera_up, view);
		//std::cout << "view matrix\n" << view << std::endl;
		view = view * (core().trackball_angle * Eigen::Scaling(core().camera_zoom * core().camera_base_zoom)
				* Eigen::Translation3f(core().camera_translation + core().camera_base_translation)).matrix() * scn->MakeTransScale()* scn->CalcParentsTrans(scn->selected_data_index).cast<float>() * scn->data().MakeTransScale();
		bool picked = igl::unproject_onto_mesh(Eigen::Vector2f(x, y), view,
			core().proj, core().viewport, scn->data().V, scn->data().F, fid, bc);
		scn->isPicked = scn->isPicked | picked;
		if (picked)
		{
			Eigen::Vector3i face = scn->data().F.row(fid);
			Eigen::Matrix3d vertices ;
			Eigen::Vector4f p,pp ;

			vertices.col(0) = scn->data().V.row(face(0));
			vertices.col(1) =  scn->data().V.row(face(1));
			vertices.col(2) = scn->data().V.row(face(2));
		
			p <<  vertices.cast<float>() * bc ,1;
			p = view * p;
			//std::cout << scn->data().V.row(face(0)) << std::endl;
			pp = core().proj * p;
			//glReadPixels(x,  y, 1, 1, GL_DEPTH_COMPONENT, GL_FLOAT, &depth);
			z = pp(2);
			return p(2);
			
		}
		return 0;
}

IGL_INLINE void Renderer::resize(GLFWwindow* window,int w, int h)
	{
		if (window) {
			glfwSetWindowSize(window, w / highdpi, h / highdpi);
		}
		post_resize(window,w, h);
	}

	IGL_INLINE void Renderer::post_resize(GLFWwindow* window, int w, int h)
	{
		if (core_list.size() == 1)
		{
	
			core().viewport = Eigen::Vector4f(0, 0, w, h);
		}
		else
		{
			// It is up to the user to define the behavior of the post_resize() function
			// when there are multiple viewports (through the `callback_post_resize` callback)
			core(left_view).viewport = Eigen::Vector4f(0, 0, w / 4, h);
			core(right_view).viewport = Eigen::Vector4f(w / 4, 0, w - (w / 4), h);

		}
		//for (unsigned int i = 0; i < plugins.size(); ++i)
		//{
		//	plugins[i]->post_resize(w, h);
		//}
		if (callback_post_resize)
		{
			callback_post_resize(window, w, h);
		}
	}

	IGL_INLINE igl::opengl::ViewerCore& Renderer::core(unsigned core_id /*= 0*/)
	{
		assert(!core_list.empty() && "core_list should never be empty");
		int core_index;
		if (core_id == 0)
			core_index = selected_core_index;
		else
			core_index = this->core_index(core_id);
		assert((core_index >= 0 && core_index < core_list.size()) && "selected_core_index should be in bounds");
		return core_list[core_index];
	}

	IGL_INLINE const igl::opengl::ViewerCore& Renderer::core(unsigned core_id /*= 0*/) const
	{
		assert(!core_list.empty() && "core_list should never be empty");
		int core_index;
		if (core_id == 0)
			core_index = selected_core_index;
		else
			core_index = this->core_index(core_id);
		assert((core_index >= 0 && core_index < core_list.size()) && "selected_core_index should be in bounds");
		return core_list[core_index];
	}

	IGL_INLINE bool Renderer::erase_core(const size_t index)
	{
		assert((index >= 0 && index < core_list.size()) && "index should be in bounds");
		//assert(data_list.size() >= 1);
		if (core_list.size() == 1)
		{
			// Cannot remove last viewport
			return false;
		}
		core_list[index].shut(); // does nothing
		core_list.erase(core_list.begin() + index);
		if (selected_core_index >= index && selected_core_index > 0)
		{
			selected_core_index--;
		}
		return true;
	}

	IGL_INLINE size_t Renderer::core_index(const int id) const {
		for (size_t i = 0; i < core_list.size(); ++i)
		{
			if (core_list[i].id == id)
				return i;
		}
		return 0;
	}

	IGL_INLINE int Renderer::append_core(Eigen::Vector4f viewport, bool append_empty /*= false*/)
	{
		core_list.push_back(core()); // copies the previous active core and only changes the viewport
		core_list.back().viewport = viewport;
		core_list.back().id = next_core_id;
		next_core_id <<= 1;
		if (!append_empty)
		{
			for (auto& data : scn->data_list)
			{
				data.set_visible(true, core_list.back().id);
				//data.copy_options(core(), core_list.back());
			}
		}
		selected_core_index = core_list.size() - 1;
		return core_list.back().id;
	}

	//IGL_INLINE void Viewer::select_hovered_core()
	//{
	//	int width_window, height_window = 800;
	//   glfwGetFramebufferSize(window, &width_window, &height_window);
	//	for (int i = 0; i < core_list.size(); i++)
	//	{
	//		Eigen::Vector4f viewport = core_list[i].viewport;

	//		if ((current_mouse_x > viewport[0]) &&
	//			(current_mouse_x < viewport[0] + viewport[2]) &&
	//			((height_window - current_mouse_y) > viewport[1]) &&
	//			((height_window - current_mouse_y) < viewport[1] + viewport[3]))
	//		{
	//			selected_core_index = i;
	//			break;
	//		}
	//	}
	//}

	bool Renderer::check_collision_cond(Eigen::AlignedBox<double, 3> box_1, Eigen::AlignedBox<double, 3> box_2) {
		double R, R_0, R_1;

		Eigen::Matrix4d A_trans = scn->data_list[0].MakeTransd();
		Eigen::Vector4d A_vec_center = Eigen::Vector4d(box_1.center()[0], box_1.center()[1], box_1.center()[2], 1);
		Eigen::Vector4d A_vec_center_trans = A_trans * A_vec_center;
		Eigen::RowVector3d P_A = Eigen::RowVector3d(A_vec_center_trans.x(), A_vec_center_trans.y(), A_vec_center_trans.z());

		Eigen::RowVector3d A_x = scn->data_list[0].GetRotation() * Eigen::Vector3d(1, 0, 0);
		Eigen::RowVector3d A_y = scn->data_list[0].GetRotation() * Eigen::Vector3d(0, 1, 0);
		Eigen::RowVector3d A_z = scn->data_list[0].GetRotation() * Eigen::Vector3d(0, 0, 1);
		double W_A = box_1.sizes()[0] / 2;
		double H_A = box_1.sizes()[1] / 2;
		double D_A = box_1.sizes()[2] / 2;
		
		Eigen::Matrix3d A;
		A << A_x[0], A_x[1], A_x[2],
			 A_y[0], A_y[1], A_y[2],
			 A_z[0], A_z[1], A_z[2];

		Eigen::Matrix4d B_trans = scn->data_list[1].MakeTransd();
		Eigen::Vector4d B_vec_center = Eigen::Vector4d(box_2.center()[0], box_2.center()[1], box_2.center()[2], 1);
		Eigen::Vector4d B_vec_center_trans = A_trans * A_vec_center;
		Eigen::RowVector3d P_B = Eigen::RowVector3d(B_vec_center_trans.x(), B_vec_center_trans.y(), B_vec_center_trans.z());

		Eigen::RowVector3d B_x = scn->data_list[1].GetRotation() * Eigen::Vector3d(1, 0, 0);
		Eigen::RowVector3d B_y = scn->data_list[1].GetRotation() * Eigen::Vector3d(0, 1, 0);
		Eigen::RowVector3d B_z = scn->data_list[1].GetRotation() * Eigen::Vector3d(0, 0, 1);
		double W_B = box_2.sizes()[0] / 2;
		double H_B = box_2.sizes()[1] / 2;
		double D_B = box_2.sizes()[2] / 2;
		Eigen::Matrix3d B;
		B << B_x[0], B_y[0], B_z[0],
			 B_x[1], B_y[1], B_z[1],
			 B_x[2], B_y[2], B_z[2];
		
		 
		Eigen::RowVector3d T = P_B - P_A;
		Eigen::Matrix3d R_mat = A * B;

		//CASE 1:
		R_0 = W_A;
		R_1 = W_B * abs(R_mat(0,0)) + H_B * abs(R_mat(0, 1)) + D_B * abs(R_mat(0, 2));
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
		R_0 = W_A * abs(R_mat(2, 1)) + D_A * abs(R_mat(0,1 ));
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

	void Renderer::check_collition()
	{
		igl::AABB<Eigen::MatrixXd, 3> node1 = scn->data_list[0].tree;
		igl::AABB<Eigen::MatrixXd, 3> node2 = scn->data_list[1].tree;

		if (check_collision_cond(node1.m_box, node2.m_box)) {
			std::cout << "Collision detected! " << std::endl;
			scn->isMoving = false;
		}
	}
