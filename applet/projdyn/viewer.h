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

typedef unsigned long size_t;

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
            if ( distance(mesh_center, mesh.position(v))> dist_max) {
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

	void updateShaderVertices(const MatrixXf& vPos) {
        m_updated_shader_verts = vPos;
        m_reupload_vertices = true;
    }

	void updateShaderRods(const MatrixXf& vTan, const MatrixXf& rPos) {
        m_updated_rods_pos = rPos;
        m_updated_rods_tangents = vTan;
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
            "a_simple_shader_for_rods",

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

        rShaderLines.init(
            "normal_shader",
            /* Vertex shader */
            "#version 330\n\n"
            "in vec3 position;\n"
            "uniform mat4 MV;\n"
            "uniform mat4 P;\n"

			"out VS_OUT {\n"
            "    mat3 normal_mat;\n"
            "} vs_out;\n"
            "void main() {\n"
            "    gl_Position = vec4(position, 1.0);\n"
            "    vs_out.normal_mat = mat3(transpose(inverse(MV)));\n"
            "}",
            /* Fragment shader */
            "#version 330\n\n"
            "out vec4 frag_color;\n"
            "void main() {\n"
            "   frag_color = vec4(0.0, 0.0, 1.0, 1.0);\n"
            "}",
            /* Geometry shader */
            "#version 330\n\n"
            "layout (triangles) in;\n"
            "layout (line_strip, max_vertices = 6) out;\n"
            "uniform mat4 MV;\n"
            "uniform mat4 P;\n"
            "in VS_OUT {\n"
            "    mat3 normal_mat;\n"
            "} gs_in[];\n"
            "void createline(int index) {\n"
            "   gl_Position = P * MV * gl_in[index].gl_Position;\n"
            "   EmitVertex();\n"
            "   gl_Position = P * MV * gl_in[index+1].gl_Position;\n"
            "   EmitVertex();\n"
            "   EndPrimitive();\n"
            "}\n"
            "void main() {\n"
            "   createline(0);\n"
            "}"
        );
    }

    ~Viewer() {
        mShader.free();
        mShaderNormals.free();
        rShader.free();
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
			//m_reupload_vertices = false;
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

        //rShader.bind();
        unsigned int numFaces;

        MatrixXf lIndices;
        lIndices.resize(3, m_updated_rods_pos.cols());
        for (size_t i = 0; i < m_updated_rods_pos.cols(); i++) {
            lIndices.coeffRef(0, i) = i;
        }

        rShaderLines.bind();
        if (m_reupload_vertices) {
            rShaderLines.uploadIndices(lIndices);
            rShaderLines.uploadAttrib("position", m_updated_rods_pos);
        }
        rShaderLines.setUniform("MV", mv);
        rShaderLines.setUniform("P", p);
        //rShaderLines.uploadAttrib("normal", normals_attrib);
        rShaderLines.drawIndexed(GL_LINES, 0, m_updated_rods_pos.cols() - 1);

        //Draw the rods if necessary
        if (false && m_reupload_vertices) {

            size_t ndivs = 3;
            size_t ncurves = m_updated_rods_pos.cols() - 1;
            size_t curveNumPts = m_updated_rods_pos.size();

            //std::vector<Vector3f []> P(new Vector3f[(ndivs + 1) * ndivs * ncurves + 1]);
            MatrixXf P_, N_, st_;
            P_.resize(3, ndivs * ndivs * ncurves + ndivs);
            //std::vector<Vector3f []> N(new Vector3f[(ndivs + 1) * ndivs * ncurves + 1]);
            N_.resize(3, ndivs * ndivs * ncurves + ndivs);
            //std::vector<Vector2f []> st(new Vector2f[(ndivs + 1) * ndivs * ncurves + 1]);
            st_.resize(2, ndivs * ndivs * ncurves + ndivs);

            size_t count = 0;
            for (size_t i = 0; i < ncurves; ++i) {
                for (size_t j = 0; j < ndivs; ++j) {
                    Vector3f p0, p1;
                    p0 = m_updated_rods_pos.col(i);
                    p1 = m_updated_rods_pos.col(i + 1);
                    float s = j / (float) ndivs;
                    Vector3f pt = evalCurve(&p0, &p1, s);
                    Vector3f tangent = m_updated_rods_tangents.col(i).normalized();

                    uint8_t maxAxis;
                    if (std::abs(tangent.x()) > std::abs(tangent.y()))
                        if (std::abs(tangent.x()) > std::abs(tangent.z()))
                            maxAxis = 0;
                        else
                            maxAxis = 2;
                    else if (std::abs(tangent.y()) > std::abs(tangent.z()))
                        maxAxis = 1;
                    else
                        maxAxis = 2;

                    Vector3f up, forward, right;

                    switch (maxAxis) {
                        case 0:
                        case 1:
                            up = tangent;
                            forward = Vector3f(0, 0, 1);
                            right = up.cross(forward);
                            forward = right.cross(up);
                            break;
                        case 2:
                            up = tangent;
                            right = Vector3f(0, 0, 1);
                            forward = right.cross(up);
                            right = up.cross(forward);
                            break;
                        default:
                            break;
                    };

                    float sNormalized = (i * ndivs + j) / float(ndivs * ncurves);
                    float rad = 2 * (1 - sNormalized);
                    for (size_t k = 0; k <= ndivs; ++k) {
                        count++;
                        float t = k / (float) ndivs;
                        float theta = t * 2 * M_PI;
                        Vector3f pc(cos(theta) * rad, 0, sin(theta) * rad);
                        float x = pc.x() * right.x() + pc.y() * up.x() + pc.z() * forward.x();
                        float y = pc.x() * right.y() + pc.y() * up.y() + pc.z() * forward.y();
                        float z = pc.x() * right.z() + pc.y() * up.z() + pc.z() * forward.z();
                        //P[i * (ndivs + 1) * ndivs + j * (ndivs + 1) + k] = Vector3f(pt.x() + x, pt.y() + y, pt.z() + z);
                        size_t index = i * ndivs * ndivs + j * ndivs + k;
                        P_.col(index) = Vector3f(pt.x() + x, pt.y() + y,pt.z() + z);
                        P_.col(index) = Vector3f(10, 10, 10);
                        //N[i * (ndivs + 1) * ndivs + j * (ndivs + 1) + k] = Vector3f(x, y, z).normalize();
                        N_.col(i * ndivs * ndivs + j * ndivs + k) = Vector3f(x, y, z).normalized();
                        //st[i * (ndivs + 1) * ndivs + j * (ndivs + 1) + k] = Vector2f(sNormalized , t);
                        st_.col(i * ndivs * ndivs + j * ndivs + k) = Vector2f(sNormalized, t);
                    }
                }
            }
            std::cout << "count: " << count << std::endl;
            P_.col(ndivs * ndivs * ncurves + ndivs - 1) = m_updated_rods_pos.col(ncurves);
            N_.col(ndivs * ndivs * ncurves + ndivs - 1) = (m_updated_rods_pos.col(ncurves - 1) -
                                                     m_updated_rods_pos.col(ncurves)).normalized();
            st_.col(ndivs * ndivs * ncurves + ndivs - 1) = Vector2f(1, 0.5);
            numFaces = ndivs * ndivs * ncurves;
            MatrixXf verts;
            verts.resize(3, numFaces);
            for (size_t i = 0; i < numFaces; ++i)
                verts.coeffRef(i) = (i < (numFaces - ndivs)) ? 4 : 3;
            MatrixXf vertIndices;
            vertIndices.resize(3, 2 * ndivs * ndivs * ncurves);
            unsigned int nf = 0, ix = 0;
            for (size_t k = 0; k < ncurves; ++k) {
                for (size_t j = 0; j < ndivs; ++j) {
                    //if (k == (ncurves - 1) && j == (ndivs - 1)) { break; }
                    for (size_t i = 0; i < ndivs; ++i) {
                        //vertIndices[ix] = nf;
                        //vertIndices[ix + 1] = nf + (ndivs + 1);
                        //vertIndices[ix + 2] = nf + (ndivs + 1) + 1;
                        //vertIndices[ix + 3] = nf + 1;
                        vertIndices.col(ix) << nf, nf + ndivs, nf + ndivs + 1;
                        ix += 1;
                        vertIndices.col(ix) << nf, nf + (ndivs + 1) + 1, nf + 1;
                        ix += 1;
                        ++nf;
                    }
                    nf++;
                }
            }

//            MatrixXu indices(3, 2); /* Draw 2 triangles */
//            indices.col(0) << 0, 1, 2;
//            indices.col(1) << 2, 3, 0;
//
//            MatrixXf positions(3, 4);
//            positions.col(0) << -1, -1, 0;
//            positions.col(1) << 1, -1, 0;
//            positions.col(2) << 1, 1, 0;
//            positions.col(3) << -1, 1, 0;

            rShader.uploadIndices(vertIndices);
            rShader.uploadAttrib("position", P_);

            m_reupload_vertices = false;
            std::cout << "Face index: \n" << vertIndices << std::endl;
            std::cout << "Positions: \n" << P_ << std::endl;
        }


//        rShader.setUniform("MV", mv);
//        rShader.setUniform("P", p);
//        Vector3f color(0.5, 0.1, 0.88);
//        rShader.setUniform("intensity", color);
//        rShader.drawIndexed(GL_TRIANGLES, 0, numFaces);
    }

    Vector3f evalCurve(const Vector3f *p0, const Vector3f *p1, const float &t)
    {
        //Linear interpolation between the two points
        Vector3f v;
        v << p0->x()*(1-t) + t*p1->x(), p0->y()*(1-t) + t*p1->y(), p0->z()*(1-t) + t*p1->z();
        return v;
    }

    void createCurveGeometry() {
        std::cout << "yay" << std::endl;
        uint32_t ndivs = 3;
        uint32_t ncurves = m_updated_rods_pos.cols() - 1;
        size_t curveNumPts = m_updated_rods_pos.size();

        //std::vector<Vector3f []> P(new Vector3f[(ndivs + 1) * ndivs * ncurves + 1]);
        MatrixXf P_, N_, st_;
        P_.resize(3, (ndivs + 1) * ndivs * ncurves + 1);
        //std::vector<Vector3f []> N(new Vector3f[(ndivs + 1) * ndivs * ncurves + 1]);
        N_.resize(3, (ndivs + 1) * ndivs * ncurves + 1);
        //std::vector<Vector2f []> st(new Vector2f[(ndivs + 1) * ndivs * ncurves + 1]);
        st_.resize(2, (ndivs + 1) * ndivs * ncurves + 1);

        for (uint32_t i = 0; i < ncurves; ++i) {
            for (uint32_t j = 0; j < ndivs; ++j) {
                Vector3f p0, p1;
                p0 = m_updated_rods_pos.col(i);
                p1 = m_updated_rods_pos.col(i+1);
                float s = j / (float)ndivs;
                Vector3f pt = evalCurve(&p0, &p1, s);
                Vector3f tangent = m_updated_rods_tangents.col(i).normalized();

                uint8_t maxAxis;
                if (std::abs(tangent.x()) > std::abs(tangent.y()))
                    if (std::abs(tangent.x()) > std::abs(tangent.z()))
                        maxAxis = 0;
                    else
                        maxAxis = 2;
                else if (std::abs(tangent.y()) > std::abs(tangent.z()))
                    maxAxis = 1;
                else
                    maxAxis = 2;

                Vector3f up, forward, right;

                switch (maxAxis) {
                    case 0:
                    case 1:
                        up = tangent;
                        forward = Vector3f(0, 0, 1);
                        right = up.cross(forward);
                        forward = right.cross(up);
                        break;
                    case 2:
                        up = tangent;
                        right = Vector3f(0, 0, 1);
                        forward = right.cross(up);
                        right = up.cross(forward);
                        break;
                    default:
                        break;
                };

                float sNormalized = (i * ndivs + j) / float(ndivs * ncurves);
                float rad = 2 * (1 - sNormalized);
                for (uint32_t k = 0; k <= ndivs; ++k) {
                    float t = k / (float)ndivs;
                    float theta = t * 2 * M_PI;
                    Vector3f pc(cos(theta) * rad, 0, sin(theta) * rad);
                    float x = pc.x() * right.x() + pc.y() * up.x() + pc.z() * forward.x();
                    float y = pc.x() * right.y() + pc.y() * up.y() + pc.z() * forward.y();
                    float z = pc.x() * right.z() + pc.y() * up.z() + pc.z() * forward.z();
                    //P[i * (ndivs + 1) * ndivs + j * (ndivs + 1) + k] = Vector3f(pt.x() + x, pt.y() + y, pt.z() + z);
                    P_.col(i * (ndivs + 1) * ndivs + j * (ndivs + 1) + k) = Vector3f(pt.x() + x, pt.y() + y, pt.z() + z);
                    //N[i * (ndivs + 1) * ndivs + j * (ndivs + 1) + k] = Vector3f(x, y, z).normalize();
                    N_.col(i * (ndivs + 1) * ndivs + j * (ndivs + 1) + k) = Vector3f(x, y, z).normalized();
                    //st[i * (ndivs + 1) * ndivs + j * (ndivs + 1) + k] = Vector2f(sNormalized , t);
                    st_.col(i * (ndivs + 1) * ndivs + j * (ndivs + 1) + k) = Vector2f(sNormalized , t);
                }
            }
        }
        P_.col((ndivs + 1) * ndivs * ncurves) = m_updated_rods_pos.col(ncurves);
        N_.col((ndivs + 1) * ndivs * ncurves) = (m_updated_rods_pos.col(ncurves - 1) - m_updated_rods_pos.col(ncurves)).normalized();
        st_.col((ndivs + 1) * ndivs * ncurves) = Vector2f(1, 0.5);
        uint32_t numFaces = ndivs * ndivs * ncurves;
        MatrixXf verts;
        verts.resize(3, numFaces);
        for (uint32_t i = 0; i < numFaces; ++i)
            verts.coeffRef(i) = (i < (numFaces - ndivs)) ? 4 : 3;
        MatrixXf vertIndices;
        vertIndices.resize(4, ndivs * ndivs * ncurves * 4 + ndivs * 3);
        uint32_t nf = 0, ix = 0;
        //TODO: Create triangles instead of quads (shouldn't be too difficult, just need to understand this piece of code
        for (uint32_t k = 0; k < ncurves; ++k) {
            for (uint32_t j = 0; j < ndivs; ++j) {
                if (k == (ncurves - 1) && j == (ndivs - 1)) { break; }
                for (uint32_t i = 0; i < ndivs; ++i) {
                    //vertIndices[ix] = nf;
                    //vertIndices[ix + 1] = nf + (ndivs + 1);
                    //vertIndices[ix + 2] = nf + (ndivs + 1) + 1;
                    //vertIndices[ix + 3] = nf + 1;
                    vertIndices.col(ix) << nf, nf + (ndivs + 1), nf + (ndivs + 1) + 1, nf + 1;
                    ix += 4;
                    ++nf;
                }
                nf++;
            }
        }

        for (uint32_t i = 0; i < ndivs; ++i) {
            vertIndices.coeffRef(ix) = nf;
            vertIndices.coeffRef(ix + 1) = (ndivs + 1) * ndivs * ncurves;
            vertIndices.coeffRef(ix + 2) = nf + 1;
            ix += 3;
            nf++;
        }

//        //Now we bind the shader !
//        mShader.bind();
//        mShader.uploadAttrib("position", P_);
//
//        Eigen::Matrix4f model, view, proj;
//        computeCameraMatrices(model, view, proj);
//
//        Matrix4f mv = view*model;
//        Matrix4f p = proj;
//
//        /* MVP uniforms */
//        mShader.setUniform("MV", mv);
//        mShader.setUniform("P", p);
//
//        /* Setup OpenGL (making sure the GUI doesn't disable these */
//        glEnable(GL_DEPTH_TEST);
//        glDisable(GL_CULL_FACE);
//
//        mShader.uploadIndices(vertIndices);
//        mShader.uploadAttrib("normal", N_);
//        Vector3f colors(0.0, 0.0, 0.88);
//        mShader.setUniform("intensity", colors);
//        mShader.drawIndexed(GL_QUADS, 0, numFaces);

        Eigen::Matrix4f model, view, proj;

        Matrix4f mv = view*model;
        Matrix4f p = proj;

        /* Setup OpenGL (making sure the GUI doesn't disable these */

        MatrixXu indices(3, 2); /* Draw 2 triangles */
        indices.col(0) << 0, 1, 2;
        indices.col(1) << 2, 3, 0;

        MatrixXf positions(3, 4);
        positions.col(0) << -1, -1, 0;
        positions.col(1) <<  1, -1, 0;
        positions.col(2) <<  1,  1, 0;
        positions.col(3) << -1,  1, 0;

        rShader.bind();
        rShader.uploadIndices(indices);
        rShader.uploadAttrib("position", positions);
        rShader.setUniform("MV", mv);
        rShader.setUniform("P", p);
        glEnable(GL_DEPTH_TEST);
        glDisable(GL_CULL_FACE);
        Vector3f colors(0.5, 0.1, 0.88);
        rShader.setUniform("intensity", colors);
        rShader.drawIndexed(GL_TRIANGLES, 0, 2);

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

private:

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
    nanogui::GLShader rShaderLines;
    nanogui::Window *window;
    nanogui::Arcball arcball;

    Point mesh_center = Point(0.0f, 0.0f, 0.0f);
    Surface_mesh mesh;

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

	// Temporary storage for updated vertex positions that have to be uploaded to the shader
	// in the next drawContents() call.
	// Can be set using updateShaderVertices()
	MatrixXf m_updated_shader_verts;

	//Same for the rods variables
	MatrixXf m_updated_rods_pos;
	MatrixXf m_updated_rods_tangents;

	// Flag that will be set to true when new vertex positions have been povided via updateShaderVertices
	// which need to be re-uploaded at the beginning of the next drawContents() call
	bool m_reupload_vertices = false;

	// Callback function to be called before each draw
	bool (*m_pre_draw_callback)(Viewer*) = nullptr;

	// Callback function to be called after loading a mesh
	bool (*m_mesh_load_callback)(Viewer*) = nullptr;


};
