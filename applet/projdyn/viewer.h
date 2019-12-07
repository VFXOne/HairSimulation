//=============================================================================
//
//   Code framework for the lecture
//
//   "Digital 3D Geometry Processing"
//
//   Gaspard Zoss
//
//   Copyright (C) 2016 by Computer Graphics and Geometry Laboratory,
//         EPF Lausanne
//
//   Edited 2017
//-----------------------------------------------------------------------------
#pragma once

#define _USE_MATH_DEFINES
#include <math.h>

#include <nanogui/opengl.h>
#include <nanogui/glutil.h>
#include <nanogui/screen.h>
#include <nanogui/window.h>
#include <nanogui/layout.h>
#include <nanogui/popupbutton.h>
#include <nanogui/label.h>
#include <nanogui/button.h>
#include <nanogui/textbox.h>
#include <nanogui/tabwidget.h>
#include <surface_mesh/Surface_mesh.h>

#if defined(__GNUC__)
#  pragma GCC diagnostic ignored "-Wmissing-field-initializers"
#endif
#if defined(_WIN32)
#  pragma warning(push)
#  pragma warning(disable: 4457 4456 4005 4312)
#endif

#if defined(_WIN32)
#  pragma warning(pop)
#endif
#if defined(_WIN32)
#  if defined(APIENTRY)
#    undef APIENTRY
#  endif
#  include <windows.h>
#endif

using std::cout;
using std::cerr;
using std::endl;
using std::string;
using std::vector;
using std::pair;
using std::to_string;
using std::min;
using std::max;
using namespace surface_mesh;
using namespace nanogui;

#ifndef WIN32
    typedef unsigned long size_t;
#endif

constexpr float INITIAL_GRAB_RADIUS = 20;

class Viewer : public nanogui::Screen {
public:
    Viewer(std::string title, bool (*pre_draw_callback)(Viewer*) = nullptr, bool (*mesh_load_callback)(Viewer*) = nullptr)
        :
        nanogui::Screen(Eigen::Vector2i(1024, 768), title) {

        initGUI();
        initShaders();

		m_pre_draw_callback = pre_draw_callback;
		m_mesh_load_callback = mesh_load_callback;
    }

    void loadMesh(string filename) {
        if (!mesh.read(filename)) {
            std::cerr << "Mesh not found, exiting." << std::endl;
            exit(-1);
        }

        cout << "Mesh "<< filename << " loaded." << endl;
        n_vertices = mesh.n_vertices();
        cout << "# of vertices : " << n_vertices << endl;
        n_faces = mesh.n_faces();
        cout << "# of faces : " << n_faces << endl;
        n_edges = mesh.n_edges();
        cout << "# of edges : " << n_edges << endl;

        mesh_center = computeCenter(&mesh);
        float dist_max = 0.0f;
        for (auto v: mesh.vertices()) {
            if (distance(mesh_center, mesh.position(v))> dist_max) {
                dist_max = distance(mesh_center, mesh.position(v));
            }
        }

        mCamera.arcball = Arcball(2.);
        mCamera.arcball.setSize(mSize);
        mCamera.modelZoom = 2/dist_max;
        mCamera.modelTranslation = -Vector3f(mesh_center.x, mesh_center.y, mesh_center.z);

		meshProcess();

		if (m_mesh_load_callback) {
			if (!m_mesh_load_callback(this)) {
				std::cout << "Error on callback after loading mesh!" << std::endl;
			}
		}
	}

    void default_rods(const size_t res = 5, const size_t num_rods = 3) {
        using_rods = true;
        r_res = res;
        r_num_rods = num_rods;
        size_t nPoints = res*num_rods;
        m_updated_rods_pos.resize(3, nPoints);
        m_updated_rods_tangents.resize(3, nPoints);
        m_updated_rods_normals.resize(3, nPoints);
        m_rod_indices.clear();
        const float distance_between_rods = 0.3;
        for (size_t j = 0; j < num_rods; j++) {
            for (int i = 0; i < res; i++) {
                m_updated_rods_pos.col(i+j*res) << distance_between_rods*j, -i,0;
                m_updated_rods_tangents.col(i+j*res) << 0,-1,0;
                m_updated_rods_normals.col(i+j*res) << -1,0,0;
            }
            m_rod_indices.push_back(j*res);
        }
    }

	void rodMesh() {
        cout << "Creating rod mesh" << endl;
        r_vertices = r_num_rods * (ndivs*(r_res - 1) + 1);
        cout << "# of vertices : " << r_vertices << endl;
        r_faces = ndivs == 2 ? r_num_rods * (2*(r_res - 2) + 1) : r_num_rods * (2*ndivs*(r_res-2) + ndivs);
        cout << "# of faces : " << r_faces << endl;
        r_edges = r_num_rods * (4*(r_res-2) + 3 + 4*ndivs);
        cout << "# of edges : " << r_edges << endl;

        rod_mesh.clear();

        MatrixXf vertices;
        vertices = createRodMesh();

        //First add every vertex to the mesh and store them in an array.
        //Then use this array to add the faces so we don't need to add them in order
        std::vector<Surface_mesh::Vertex> v_m;
        for (size_t i = 0; i < r_vertices; i++) {
            Surface_mesh::Vertex v = rod_mesh.add_vertex(Point(vertices.coeff(0,i), vertices.coeff(1, i), vertices.coeff(2, i)));
            v_m.push_back(v);
        }

        Surface_mesh::Vertex u, v, w, x;
        for (size_t ind = 0; ind < m_rod_face_indices.size()-1; ind++) {
            size_t rod_index = m_rod_face_indices.at(ind);
            size_t next_index = m_rod_face_indices.at(ind+1);

            size_t nfaces = next_index - rod_index;
            for (size_t i = 0; i < (r_faces/r_num_rods-ndivs)/2; i++) {
                u = v_m.at(rod_index+i);
                v = i == 0 or (i+1) % ndivs != 0 ? v_m.at(rod_index+i+1) : v_m.at(rod_index+i-(ndivs-1));
                w = v_m.at(rod_index+i+ndivs);
                x = i == 0 or (i+1) % ndivs != 0 ? v_m.at(rod_index+i+ndivs+1) : v_m.at(rod_index+i+1);
                rod_mesh.add_triangle(u, v, w);
                rod_mesh.add_triangle(v, x, w);
            }
            for (size_t i = 0; i < ndivs; i++) {
                u = v_m.at(next_index-1);
                v = v_m.at(next_index-ndivs+i-1);
                w = i == ndivs-1 ? v_m.at(next_index-ndivs-1) : v_m.at(next_index-ndivs+i);
                rod_mesh.add_triangle(v, w, u);
            }
        }

        mesh_center = computeCenter(&rod_mesh);
        float dist_max = 0.0f;
        for (auto v: rod_mesh.vertices()) {
            if ( distance(mesh_center, rod_mesh.position(v))> dist_max) {
                dist_max = distance(mesh_center, rod_mesh.position(v));
            }
        }

        mCamera.arcball = Arcball(2.);
        mCamera.arcball.setSize(mSize);
        mCamera.modelZoom = 2/dist_max;
        mCamera.modelTranslation = -Vector3f(mesh_center.x, mesh_center.y, mesh_center.z);

        rodProcess();

        if (m_mesh_load_callback) {
            if (!m_mesh_load_callback(this)) {
                std::cout << "Error on callback after loading mesh!" << std::endl;
            }
        }
    }

    void add_rods_on_ball(float radius, size_t res = 5, size_t num_rods = 1, float seg_length = 0.5) {
        using_rods = true;
        r_res = res;
        r_num_rods = num_rods;
        size_t nPoints = res*num_rods;
        m_updated_rods_pos.resize(3, nPoints);
        m_updated_rods_tangents.resize(3, nPoints);
        m_updated_rods_normals.resize(3, nPoints);
        m_rod_indices.clear();
        Vector3f default_tangent = Vector3f::UnitX();
        Vector3f default_normal = Vector3f::UnitY();
        auto randVal = [](float max){ return float(rand()) / float(RAND_MAX) * max; };
        for (size_t j = 0; j < num_rods; j++) {
            Quaternionf rotQuat;
            rotQuat.FromTwoVectors(default_tangent, Vector3f(randVal(1), randVal(1), randVal(1)));
            rotQuat.normalize();
            Vector3f tangent = rotQuat.toRotationMatrix() * default_tangent;
            tangent = default_tangent;
            for (size_t i = 0; i < res; i++) {
                Vector3f p = tangent * radius + tangent * i * seg_length;
                m_updated_rods_pos.col(i+j*res) << p;
                m_updated_rods_tangents.col(i+j*res) << p.normalized();
                m_updated_rods_normals.col(i+j*res) << rotQuat.toRotationMatrix() * default_normal;
            }
            m_rod_indices.push_back(j*res);
        }

        rodMesh();
        loadMesh("../data/small_sphere.obj");
        Point center = computeCenter(&mesh);
        center += Point(1+radius, -2, 0);
        for (auto v : mesh.vertices()) {
            mesh.position(v) = (mesh.position(v) - center) * radius + center;
        }
    }

	void updateShaderVertices(const MatrixXf& vPos) {
        m_updated_shader_verts = vPos;
        m_reupload_vertices = true;
    }

	void updateShaderRods(const MatrixXf& rPos, const MatrixXf& rTan, const MatrixXf& rNorm, const std::vector<size_t> rod_indices) {
        m_updated_rods_pos = rPos;
        m_updated_rods_tangents = rTan;
        m_updated_rods_normals = rNorm;
        m_updated_rods_verts = createRodMesh();
        m_rod_indices = rod_indices;
        m_reupload_vertices = true;
    }

    void meshProcess() {
        Point default_normal(0.0, 1.0, 0.0);
        Surface_mesh::Vertex_property<Point> vertex_normal =
                mesh.vertex_property<Point>("v:normal");
        mesh.update_face_normals();
        mesh.update_vertex_normals();

		int j = 0;
        MatrixXf mesh_points(3, n_vertices);
        MatrixXu indices(3, n_faces);

        for(auto f: mesh.faces()) {
            vector<float> vv(3.0f);
            int k = 0;
            for (auto v: mesh.vertices(f)) {
                vv[k] = v.idx();
                ++k;
            }
            indices.col(j) << vv[0], vv[1], vv[2];
            ++j;
        }

        // Create big matrices to send the data to the GPU with the required
        // format
        MatrixXf normals_attrib(3, n_vertices);

        j = 0;
        for (auto v: mesh.vertices()) {
            mesh_points.col(j) << mesh.position(v).x,
                                  mesh.position(v).y,
                                  mesh.position(v).z;

            normals_attrib.col(j) << vertex_normal[v].x,
                                     vertex_normal[v].y,
                                     vertex_normal[v].z;
            ++j;
        }

        mShader.bind();
        mShader.uploadIndices(indices);
        mShader.uploadAttrib("position", mesh_points);
        mShader.uploadAttrib("normal", normals_attrib);
        mShader.setUniform("intensity", Vector3f(0.98, 0.59, 0.04));

        mShaderNormals.bind();
        mShaderNormals.uploadIndices(indices);
        mShaderNormals.uploadAttrib("position", mesh_points);
        mShaderNormals.uploadAttrib("normal", normals_attrib);
    }

    void rodProcess() {
        Point default_normal(0.0, 1.0, 0.0);
        Surface_mesh::Vertex_property<Point> vertex_normal =
                rod_mesh.vertex_property<Point>("v:normal");
        rod_mesh.update_face_normals();
        rod_mesh.update_vertex_normals();

        int j = 0;
        MatrixXf mesh_points(3, r_vertices);
        MatrixXu indices(3, r_faces);

        for(auto f: rod_mesh.faces()) {
            vector<float> vv(3.0f);
            int k = 0;
            for (auto v: rod_mesh.vertices(f)) {
                vv[k] = v.idx();
                ++k;
            }
            indices.col(j) << vv[0], vv[1], vv[2];
            ++j;
        }

        // Create big matrices to send the data to the GPU with the required
        // format
        MatrixXf normals_attrib(3, r_vertices);

        j = 0;
        for (auto v: rod_mesh.vertices()) {
            mesh_points.col(j) << rod_mesh.position(v).x,
                    rod_mesh.position(v).y,
                    rod_mesh.position(v).z;

            normals_attrib.col(j) << vertex_normal[v].x,
                    vertex_normal[v].y,
                    vertex_normal[v].z;
            ++j;
        }

        rShader.bind();
        rShader.uploadIndices(indices);
        rShader.uploadAttrib("position", mesh_points);
        rShader.uploadAttrib("normal", normals_attrib);
        rShader.setUniform("intensity", Vector3f(0.0, 0.0, 0.88));

        rShaderNormals.bind();
        rShaderNormals.uploadIndices(indices);
        rShaderNormals.uploadAttrib("position", mesh_points);
        rShaderNormals.uploadAttrib("normal", normals_attrib);
    }

    void initGUI() {
        window = new Window(this, "Controls");
        window->setPosition(Vector2i(15, 15));
        window->setLayout(new GroupLayout());

        PopupButton *popupBtn = new PopupButton(window, "Open a mesh", ENTYPO_ICON_EXPORT);
        Popup *popup = popupBtn->popup();
		popup->setLayout(new GroupLayout());

        Button* b = new Button(popup, "Eight");
        b->setCallback([this,popupBtn]() {
            loadMesh("../data/eight.off");
            meshProcess();
            popupBtn->setPushed(false);
        });

        b = new Button(popup, "Small Sphere");
        b->setCallback([this,popupBtn]() {
            loadMesh("../data/small_sphere.obj");
            meshProcess();
            popupBtn->setPushed(false);
        });

        b = new Button(popup, "Default rods");
        b->setCallback([this,popupBtn]() {
            default_rods();
            rodMesh();
            rodProcess();
            popupBtn->setPushed(false);
        });

        new Label(window, "Display Control", "sans-bold");

        b = new Button(window, "Wireframe");
        b->setFlags(Button::ToggleButton);
        b->setChangeCallback([this](bool wireframe) {
            this->wireframe =! this->wireframe;
        });
        b = new Button(window, "Normals");
        b->setFlags(Button::ToggleButton);
        b->setChangeCallback([this](bool normals) {
            this->normals =! this->normals;
        });

        performLayout();
    }

    void initShaders() {
        // Shaders
        mShader.init(
            "a_simple_shader",

            /* Vertex shader */
            "#version 330\n"
            "uniform mat4 MV;\n"
            "uniform mat4 P;\n"
            "uniform vec3 intensity;\n"

            "in vec3 position;\n"
			"in vec3 normal;\n"

            "out vec3 fcolor;\n"
            "out vec3 fnormal;\n"
            "out vec3 view_dir;\n"
            "out vec3 light_dir;\n"

            "void main() {\n"
            "    vec4 vpoint_mv = MV * vec4(position, 1.0);\n"
            "    gl_Position = P * vpoint_mv;\n"
            "    fcolor = intensity;\n"
			"    fnormal = mat3(transpose(inverse(MV))) * normal;\n"
            "    light_dir = vec3(0.0, 3.0, 3.0) - vpoint_mv.xyz;\n"
            "    view_dir = -vpoint_mv.xyz;\n"
            "}",

            /* Fragment shader */
            "#version 330\n"
            "uniform vec3 intensity;\n"

            "in vec3 fcolor;\n"
            "in vec3 fnormal;\n"
            "in vec3 view_dir;\n"
            "in vec3 light_dir;\n"

            "out vec4 color;\n"

            "void main() {\n"
            "    vec3 c = vec3(0.0);\n"
            "    c += vec3(1.0)*vec3(0.18, 0.1, 0.1);\n"
            "    vec3 n = normalize(fnormal);\n"
            "    vec3 v = normalize(view_dir);\n"
            "    vec3 l = normalize(light_dir);\n"
            "    float lambert = dot(n,l);\n"
            "    if(lambert > 0.0) {\n"
            "        c += vec3(1.0)*vec3(0.9, 0.5, 0.5)*lambert;\n"
            "        vec3 v = normalize(view_dir);\n"
            "        vec3 r = reflect(-l,n);\n"
            "        c += vec3(1.0)*vec3(0.8, 0.8, 0.8)*pow(max(dot(r,v), 0.0), 90.0);\n"
            "    }\n"
            "    c *= fcolor;\n"
            "    if (intensity == vec3(0.0)) {\n"
            "        c = intensity;\n"
            "    }\n"
            "    color = vec4(c, 1.0);\n"
            "}"
        );

        rShader.init(
            "a_simple_rod_shader",

            /* Vertex shader */
            "#version 330\n"
            "uniform mat4 MV;\n"
            "uniform mat4 P;\n"
            "uniform vec3 intensity;\n"

            "in vec3 position;\n"
			"in vec3 normal;\n"

            "out vec3 fcolor;\n"
            "out vec3 fnormal;\n"
            "out vec3 view_dir;\n"
            "out vec3 light_dir;\n"

            "void main() {\n"
            "    vec4 vpoint_mv = MV * vec4(position, 1.0);\n"
            "    gl_Position = P * vpoint_mv;\n"
            "    fcolor = intensity;\n"
			"    fnormal = mat3(transpose(inverse(MV))) * normal;\n"
            "    light_dir = vec3(0.0, 3.0, 3.0) - vpoint_mv.xyz;\n"
            "    view_dir = -vpoint_mv.xyz;\n"
            "}",

            /* Fragment shader */
            "#version 330\n"
            "uniform vec3 intensity;\n"

            "in vec3 fcolor;\n"
            "in vec3 fnormal;\n"
            "in vec3 view_dir;\n"
            "in vec3 light_dir;\n"

            "out vec4 color;\n"

            "void main() {\n"
            "    vec3 c = vec3(0.0);\n"
            "    c += vec3(1.0)*vec3(0.18, 0.1, 0.1);\n"
            "    vec3 n = normalize(fnormal);\n"
            "    vec3 v = normalize(view_dir);\n"
            "    vec3 l = normalize(light_dir);\n"
            "    float lambert = dot(n,l);\n"
            "    if(lambert > 0.0) {\n"
            "        c += vec3(1.0)*vec3(0.9, 0.5, 0.5)*lambert;\n"
            "        vec3 v = normalize(view_dir);\n"
            "        vec3 r = reflect(-l,n);\n"
            "        c += vec3(1.0)*vec3(0.8, 0.8, 0.8)*pow(max(dot(r,v), 0.0), 90.0);\n"
            "    }\n"
            "    c *= fcolor;\n"
            "    if (intensity == vec3(0.0)) {\n"
            "        c = intensity;\n"
            "    }\n"
            "    color = vec4(c, 1.0);\n"
            "}"
        );

        mShaderNormals.init(
            "normal_shader",
            /* Vertex shader */
            "#version 330\n\n"
            "in vec3 position;\n"
			"in vec3 normal;\n"
            "uniform mat4 MV;\n"
            "uniform mat4 P;\n"

			"out VS_OUT {\n"
            "    mat3 normal_mat;\n"
            "    vec3 normal;\n"
            "} vs_out;\n"
            "void main() {\n"
            "    gl_Position = vec4(position, 1.0);\n"
			"    vs_out.normal = normal;\n"
            "    vs_out.normal_mat = mat3(transpose(inverse(MV)));\n"
            "}",
            /* Fragment shader */
            "#version 330\n\n"
            "out vec4 frag_color;\n"
            "void main() {\n"
            "   frag_color = vec4(0.0, 1.0, 0.0, 1.0);\n"
            "}",
            /* Geometry shader */
            "#version 330\n\n"
            "layout (triangles) in;\n"
            "layout (line_strip, max_vertices = 6) out;\n"
            "uniform mat4 MV;\n"
            "uniform mat4 P;\n"
            "in VS_OUT {\n"
            "    mat3 normal_mat;\n"
            "    vec3 normal;\n"
            "} gs_in[];\n"
            "void createline(int index) {\n"
            "   gl_Position = P * MV * gl_in[index].gl_Position;\n"
            "   EmitVertex();\n"
            "   vec4 normal_mv = vec4(normalize(gs_in[index].normal_mat *\n"
            "                                   gs_in[index].normal), 1.0f);\n"
            "   gl_Position = P * (MV * gl_in[index].gl_Position\n"
            "                      + normal_mv * 0.035f);\n"
            "   EmitVertex();\n"
            "   EndPrimitive();\n"
            "}\n"
            "void main() {\n"
            "   createline(0);\n"
            "   createline(1);\n"
            "   createline(2);\n"
            "}"
        );

        rShaderNormals.init(
            "normal_rod_shader",
            /* Vertex shader */
            "#version 330\n\n"
            "in vec3 position;\n"
			"in vec3 normal;\n"
            "uniform mat4 MV;\n"
            "uniform mat4 P;\n"

			"out VS_OUT {\n"
            "    mat3 normal_mat;\n"
            "    vec3 normal;\n"
            "} vs_out;\n"
            "void main() {\n"
            "    gl_Position = vec4(position, 1.0);\n"
			"    vs_out.normal = normal;\n"
            "    vs_out.normal_mat = mat3(transpose(inverse(MV)));\n"
            "}",
            /* Fragment shader */
            "#version 330\n\n"
            "out vec4 frag_color;\n"
            "void main() {\n"
            "   frag_color = vec4(0.0, 1.0, 0.0, 1.0);\n"
            "}",
            /* Geometry shader */
            "#version 330\n\n"
            "layout (triangles) in;\n"
            "layout (line_strip, max_vertices = 6) out;\n"
            "uniform mat4 MV;\n"
            "uniform mat4 P;\n"
            "in VS_OUT {\n"
            "    mat3 normal_mat;\n"
            "    vec3 normal;\n"
            "} gs_in[];\n"
            "void createline(int index) {\n"
            "   gl_Position = P * MV * gl_in[index].gl_Position;\n"
            "   EmitVertex();\n"
            "   vec4 normal_mv = vec4(normalize(gs_in[index].normal_mat *\n"
            "                                   gs_in[index].normal), 1.0f);\n"
            "   gl_Position = P * (MV * gl_in[index].gl_Position\n"
            "                      + normal_mv * 0.05f);\n"
            "   EmitVertex();\n"
            "   EndPrimitive();\n"
            "}\n"
            "void main() {\n"
            "   createline(0);\n"
            "   createline(1);\n"
            "   createline(2);\n"
            "}"
        );
    }

    ~Viewer() {
        mShader.free();
        mShaderNormals.free();
        rShader.free();
        rShaderNormals.free();
    }

    Point computeCenter(Surface_mesh *mesh) {
        Point center = Point(0.0f);

        for (auto v: mesh->vertices()) {
            center += mesh->position(v);
        }

        return center/mesh->n_vertices();
    }

    virtual bool keyboardEvent(int key, int scancode, int action, int modifiers) {
        if (Screen::keyboardEvent(key, scancode, action, modifiers)) {
            return true;
        }
        if (key == GLFW_KEY_ESCAPE && action == GLFW_PRESS) {
            setVisible(false);
            return true;
        }
        return false;
    }

    virtual void draw(NVGcontext *ctx) {
        /* Draw the user interface */
        Screen::draw(ctx);
    }

    Vector2f getScreenCoord() {
        Vector2i pos = mousePos();
        return Vector2f(2.0f * (float)pos.x() / width() - 1.0f,
                        1.0f - 2.0f * (float)pos.y() / height());
    }

    void repaint() {
        //glfwPostEmptyEvent();
    }

    virtual void drawContents() {
        using namespace nanogui;

		/* Pre-draw callback */
		if (m_pre_draw_callback) {
			if (!m_pre_draw_callback(this)) {
				return;
			}
		}

		/* Draw the window contents using OpenGL */
		mShader.bind();

		if (m_reupload_vertices) {
			mShader.uploadAttrib("position", m_updated_shader_verts);
		}

        Eigen::Matrix4f model, view, proj;
        computeCameraMatrices(model, view, proj);

        Matrix4f mv = view*model;
        Matrix4f p = proj;

        /* MVP uniforms */
        mShader.setUniform("MV", mv);
        mShader.setUniform("P", p);

        /* Setup OpenGL (making sure the GUI doesn't disable these */
        glEnable(GL_DEPTH_TEST);
        glDisable(GL_CULL_FACE);

        /* Render everything */
        if (wireframe) {
            glEnable(GL_POLYGON_OFFSET_FILL);
            glPolygonOffset(1.0, 1.0);
            glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
        }

        Vector3f colors(0.98, 0.59, 0.04);
        mShader.setUniform("intensity", colors);
        mShader.drawIndexed(GL_TRIANGLES, 0, n_faces);

        if (wireframe) {
            glDisable(GL_POLYGON_OFFSET_FILL);
            glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
            colors << 0.0, 0.0, 0.0;
            mShader.setUniform("intensity", colors);
            mShader.drawIndexed(GL_TRIANGLES, 0, n_faces);
            glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
        }

        if (normals) {
            mShaderNormals.bind();
            mShaderNormals.setUniform("MV", mv);
            mShaderNormals.setUniform("P", p);
            mShaderNormals.drawIndexed(GL_TRIANGLES, 0, n_faces);
        }


		/* Draw the rods */
		rShader.bind();

        MatrixXf normals_attrib(3, r_vertices);

		if (m_reupload_vertices) {
            Surface_mesh::Vertex_property<Point> vertex_normal =
                    rod_mesh.vertex_property<Point>("v:normal");
            rod_mesh.update_face_normals();
            rod_mesh.update_vertex_normals();

            size_t j = 0;
            for (auto v: rod_mesh.vertices()) {
                normals_attrib.col(j) << vertex_normal[v].x,
                        vertex_normal[v].y,
                        vertex_normal[v].z;
                ++j;
            }

			rShader.uploadAttrib("position", m_updated_rods_verts);
            rShader.uploadAttrib("normal", normals_attrib);

			m_reupload_vertices = false;
		}

        /* MVP uniforms */
        rShader.setUniform("MV", mv);
        rShader.setUniform("P", p);

        /* Setup OpenGL (making sure the GUI doesn't disable these */
        glEnable(GL_DEPTH_TEST);
        glDisable(GL_CULL_FACE);

        /* Render everything */
        if (wireframe) {
            glEnable(GL_POLYGON_OFFSET_FILL);
            glPolygonOffset(1.0, 1.0);
            glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
        }

        colors << 0.0, 0.0, 0.88;
        rShader.setUniform("intensity", colors);
        rShader.drawIndexed(GL_TRIANGLES, 0, r_faces);

        if (wireframe) {
            glDisable(GL_POLYGON_OFFSET_FILL);
            glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
            colors << 0.0, 0.0, 0.0;
            rShader.setUniform("intensity", colors);
            rShader.drawIndexed(GL_TRIANGLES, 0, r_faces);
            glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
        }

        if (normals) {
            rShaderNormals.bind();
            rShaderNormals.uploadAttrib("position", m_updated_rods_verts);
            rShaderNormals.uploadAttrib("normal", normals_attrib);
            rShaderNormals.setUniform("MV", mv);
            rShaderNormals.setUniform("P", p);
            rShaderNormals.drawIndexed(GL_TRIANGLES, 0, r_faces);
        }
    }

    bool scrollEvent(const Vector2i &p, const Vector2f &rel) {
        if (!Screen::scrollEvent(p, rel)) {
            mCamera.zoom = max(0.1, mCamera.zoom * (rel.y() > 0 ? 1.1 : 0.9));
            repaint();
        }
        return true;
    }

    void reset() {
        // reset all the components
        // recenter the mesh (maybe keep the original mesh somewhere so that if
        // we modify - smooth or else - it we can restore the original one.)
    }

    bool mouseMotionEvent(const Vector2i &p, const Vector2i &rel,
                          int button, int modifiers) {
        if (!Screen::mouseMotionEvent(p, rel, button, modifiers)) {
            if (m_isGrabbing) {
                grabMove(p);
            } else if (mCamera.arcball.motion(p)) {
                repaint();
            } else if (mTranslate) {
                Eigen::Matrix4f model, view, proj;
                computeCameraMatrices(model, view, proj);
                float zval = nanogui::project(Vector3f(mesh_center.x,
                                                       mesh_center.y,
                                                       mesh_center.z),
                                              view * model, proj, mSize).z();
                Eigen::Vector3f pos1 = nanogui::unproject(
                        Eigen::Vector3f(p.x(), mSize.y() - p.y(), zval), view * model, proj, mSize);
                Eigen::Vector3f pos0 = nanogui::unproject(
                        Eigen::Vector3f(mTranslateStart.x(), mSize.y() -
                           mTranslateStart.y(), zval), view * model, proj, mSize);
                mCamera.modelTranslation = mCamera.modelTranslation_start + (pos1-pos0);
                repaint();
            }
        }
        return true;
    }

    bool mouseButtonEvent(const Vector2i &p, int button, bool down, int modifiers) {
        if (!Screen::mouseButtonEvent(p, button, down, modifiers)) {
            // In grab mode the left mouse click selects/drags vertices of the mesh,
            if (button == GLFW_MOUSE_BUTTON_1) {
                if (modifiers == GLFW_MOD_ALT) {
                    grabButton(p, down);
                } else if ( modifiers == 0) {
                    mCamera.arcball.button(p, down);
                }
            }
            if (button == GLFW_MOUSE_BUTTON_2 ||
                  (button == GLFW_MOUSE_BUTTON_1 && modifiers == GLFW_MOD_SHIFT)) {
                mCamera.modelTranslation_start = mCamera.modelTranslation;
                mTranslate = true;
                mTranslateStart = p;
            }
        }

        if (button == GLFW_MOUSE_BUTTON_1 && !down) {
            grabButton(p, down);
            mCamera.arcball.button(p, down);
        }

        if (!down) {
            mDrag = false;
            mTranslate = false;
        }
        return true;
    }

    void setGrabCallbacks(void (*grabCallback)(const std::vector<size_t>&,const std::vector<Vector3f>&), void (*grabReleaseCallback)()) {
        m_grab_callback = grabCallback;
        m_grab_release_callback = grabReleaseCallback;
    }

	Surface_mesh* getMesh() {
		return &mesh;
	}

	bool is_using_rods() {
        return using_rods;
    }

    MatrixXf* getRodsPos() {
        return &m_updated_rods_pos;
    }

    std::vector<size_t> getRodIndices() {
        return m_rod_indices;
    }

private:

    MatrixXf createRodMesh() {
        MatrixXf vertices;
        vertices.resize(3, r_vertices);
        size_t j = 0;
        m_rod_face_indices.clear();
        m_rod_face_indices.push_back(0);
        for (size_t ind = 0; ind < m_rod_indices.size(); ind++) {
            size_t rod_index = m_rod_indices.at(ind);
            size_t next_index = ind == m_rod_indices.size()-1 ? m_updated_rods_pos.cols() : m_rod_indices.at(ind+1);

            size_t nPoints = next_index - rod_index;
            for (size_t i = rod_index; i < next_index - 1; i++) {
                const float width = m_thickness * (nPoints - (i - rod_index)) / nPoints;
                Vector3f tangent = m_updated_rods_tangents.col(i).normalized();
                Vector3f normal = m_updated_rods_normals.col(i).normalized();

                for (size_t k = 0; k < ndivs; k++) {
                    float theta = k / (float)ndivs * 2 * M_PI;

                    Eigen::AngleAxisf rotMat(theta, tangent);

                    vertices.col(j++) = m_updated_rods_pos.col(i) + width * rotMat._transformVector(normal);
                }
            }
            vertices.col(j++) = m_updated_rods_pos.col(next_index - 1);
            m_rod_face_indices.push_back(j);
        }
        return vertices;
    }

    void grabButton(const Vector2i& mousePos, bool down) {
        if (m_isGrabbing && !down) {
            grabFree();
        }

        if (!m_isGrabbing && down) {
            grabSelect(mousePos);
        }
    }

    void grabFree() {
        m_isGrabbing = false;
        m_grabSelection.clear();
        if (m_grab_release_callback) m_grab_release_callback();
    }

    void grabSelect(const Vector2i& mousePos)
    {
    	if ( !m_currentMVP.hasNaN() && m_updated_shader_verts.rows() > 0 && m_grab_callback) {
    		int winWidth = mSize.x();
    		int winHeight = mSize.y();

    		m_grabSelection.clear();
    		m_grabVertPos.clear();
    		m_grabOrigPos.setZero();
    		m_grabOrigProj.setZero();
            // Check which projections of used vertices are close to mouse
    		for (size_t i = 0; i < m_updated_shader_verts.cols(); i++) {
    			Vector3f vert = m_updated_shader_verts.col(i);
    			Vector4f vertProjection(vert(0), vert(1), vert(2), 1.f);
    			vertProjection = m_currentMVP * vertProjection;
    			Vector2i screenPos;
    			screenPos(0) = (((vertProjection(0) / vertProjection(3)) + 1.f) / 2.f) * (float)winWidth;
    			screenPos(1) = (((vertProjection(1) / -vertProjection(3)) + 1.f) / 2.f) * (float)winHeight;
    			float distToMouse = (screenPos - mousePos).norm();
    			if (distToMouse < m_grabRadius) {
    				m_grabSelection.push_back(i);
    				m_grabVertPos.push_back(vert);
    				m_grabOrigProj += vertProjection;
    			}
    		}
    		if (m_grabSelection.size() > 0) {
    			m_grabOrigProj *= 1.f / (float)m_grabSelection.size();
                m_grabOrigProj(0) = ((2.f * (float)mousePos.x() / (float)winWidth) - 1.f) * m_grabOrigProj(3);
                m_grabOrigProj(1) = ((2.f * (float)mousePos.y() / (float)winHeight) - 1.f) * (-m_grabOrigProj(3));
                m_grabOrigPos = (m_currentMVP.inverse() * m_grabOrigProj).template segment<3>(0);

                m_isGrabbing = true;
    			m_grab_callback(m_grabSelection, m_grabVertPos);
    		}
    		else {
    			m_isGrabbing = false;
    		}
    	}

    	m_lastMousePosition = mousePos;
    }


    void grabMove(const Vector2i& mousePos)
    {
        if ( m_isGrabbing && m_grabSelection.size() > 0 && !m_currentMVP.hasNaN() && m_updated_shader_verts.rows() > 0 && m_grab_callback) {
            int winWidth = mSize.x();
    		int winHeight = mSize.y();

    		Vector4f newProj = m_grabOrigProj;
    		newProj(0) = ((2.f * (float)mousePos.x() / (float)winWidth) - 1.f) * newProj(3);
    		newProj(1) = ((2.f * (float)mousePos.y() / (float)winHeight) - 1.f) * (-newProj(3));
    		Vector3f newPos = (m_currentMVP.inverse() * newProj).template segment<3>(0);

    		Vector3f trans = newPos - m_grabOrigPos;

    		std::vector<Vector3f> newVertPos;
    		for (size_t v = 0; v < m_grabVertPos.size(); v++) {
    			newVertPos.push_back(m_grabVertPos[v] + trans);
    		}
    		m_grab_callback(m_grabSelection, newVertPos);
    	}

    	m_lastMousePosition = mousePos;
    }

    struct CameraParameters {
        nanogui::Arcball arcball;
        float zoom = 1.0f, viewAngle = 45.0f;
        float dnear = 0.05f, dfar = 100.0f;
        Eigen::Vector3f eye = Eigen::Vector3f(0.0f, 0.0f, 5.0f);
        Eigen::Vector3f center = Eigen::Vector3f(0.0f, 0.0f, 0.0f);
        Eigen::Vector3f up = Eigen::Vector3f(0.0f, 1.0f, 0.0f);
        Eigen::Vector3f modelTranslation = Eigen::Vector3f::Zero();
        Eigen::Vector3f modelTranslation_start = Eigen::Vector3f::Zero();
        float modelZoom = 1.0f;
    };

    CameraParameters mCamera;
    bool mTranslate = false;
    bool mDrag = false;
    Vector2i mTranslateStart = Vector2i(0,0);

    void computeCameraMatrices(Eigen::Matrix4f &model,
                               Eigen::Matrix4f &view,
                               Eigen::Matrix4f &proj) {

        view = nanogui::lookAt(mCamera.eye, mCamera.center, mCamera.up);

        float fH = std::tan(mCamera.viewAngle / 360.0f * M_PI) * mCamera.dnear;
        float fW = fH * (float) mSize.x() / (float) mSize.y();

        proj = nanogui::frustum(-fW, fW, -fH, fH, mCamera.dnear, mCamera.dfar);
        model = mCamera.arcball.matrix();
        model *= nanogui::scale(Eigen::Vector3f::Constant(mCamera.zoom * mCamera.modelZoom));
        model *= nanogui::translate(mCamera.modelTranslation);

        m_currentMVP = proj * view * model;
    }

    // Variables for the viewer
    nanogui::GLShader mShader;
    nanogui::GLShader rShader;
    nanogui::GLShader mShaderNormals;
    nanogui::GLShader rShaderNormals;
    nanogui::Window *window;
    nanogui::Arcball arcball;

    Point mesh_center = Point(0.0f, 0.0f, 0.0f);
    Surface_mesh mesh, rod_mesh;

    enum COLOR_MODE : int { NORMAL = 0, VALENCE = 1, CURVATURE = 2 };
    enum CURVATURE_TYPE : int { UNIMEAN = 2, LAPLACEBELTRAMI = 3, GAUSS = 4 };

    // Boolean for the viewer
    bool wireframe = false;
    bool normals = false;

    // Grab variables
    bool m_isGrabbing = false;
    std::vector<size_t> m_grabSelection;
    Eigen::Matrix4f m_currentMVP; // Current Model-View-Projection
    float m_grabRadius = INITIAL_GRAB_RADIUS;
    std::vector<Vector3f> m_grabVertPos;
    Vector3f m_grabOrigPos;
    Vector4f m_grabOrigProj;
    Vector2i m_lastMousePosition;
    void (*m_grab_callback)(const std::vector<size_t>&,const std::vector<Vector3f>&) = nullptr;
    void (*m_grab_release_callback)() = nullptr;

    PopupButton *popupCurvature;
    FloatBox<float>* coefTextBox;
    IntBox<int>* iterationTextBox;

    // Mesh informations
    int n_vertices = 0;
    int n_faces = 0;
    int n_edges = 0;

    // Rod informations
    size_t r_res = 0;
    size_t r_num_rods = 0;
    const size_t ndivs = 8;
    const float m_thickness = 0.1f;
    size_t r_vertices = 0;
    size_t r_faces = 0;
    size_t r_edges = 0;

	// Temporary storage for updated vertex positions that have to be uploaded to the shader
	// in the next drawContents() call.
	// Can be set using updateShaderVertices()
	MatrixXf m_updated_shader_verts;

	//Same for the rods variables
	bool using_rods = false;
	MatrixXf m_updated_rods_pos;
	MatrixXf m_updated_rods_tangents;
	MatrixXf m_updated_rods_normals;
	MatrixXf m_updated_rods_verts;
	std::vector<size_t> m_rod_indices;
	std::vector<size_t> m_rod_face_indices;

	// Flag that will be set to true when new vertex positions have been povided via updateShaderVertices
	// which need to be re-uploaded at the beginning of the next drawContents() call
	bool m_reupload_vertices = false;

	// Callback function to be called before each draw
	bool (*m_pre_draw_callback)(Viewer*) = nullptr;

	// Callback function to be called after loading a mesh
	bool (*m_mesh_load_callback)(Viewer*) = nullptr;


};
