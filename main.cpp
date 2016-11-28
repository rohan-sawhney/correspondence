#include "Mesh.h"
#include "RenderData.h"
#include "Camera.h"

#define ESCAPE 27
#define DIGIT_OFFSET 48

int gridX = 800;
int gridY = 600;

GLuint transformUbo;
GLuint lightUbo;

std::vector<std::string> paths;
std::string shaderPath;
Shader meshShader;
Shader normalShader;
Shader wireframeShader;
Shader pointShader;

Camera camera;
int t = 0;
float lastTime = 0.0;
float dt = 0.0;
float lastX = 0.0, lastY = 0.0;
bool keys[256];
bool firstMouse = true;

std::vector<Mesh> meshes = {Mesh(), Mesh()};
std::vector<GLMesh> glMeshes = {GLMesh(&meshes[0]), GLMesh(&meshes[1])};

std::vector<std::unordered_map<int, bool>> featureMaps(2);
std::vector<std::vector<GLPoint>> glPoints(2);

const Eigen::Vector3f defaultColor(1.0, 0.5, 0.2);
std::vector<std::vector<Eigen::Vector3f>> colors(meshes.size());

const Eigen::Vector3f lightPosition(0.0, 3.0, -3.0);
const Eigen::Vector3f lightColor(1.0, 1.0, 1.0);

bool success = true;
bool showNormals = false;
bool showWireframe = false;
bool showDescriptor = false;
bool showFeaturePoints = false;
bool computedDescriptor = false;

void setupShaders()
{
    meshShader.setup(shaderPath, "Model.vert", "", "Model.frag");
    normalShader.setup(shaderPath, "Normal.vert", "Normal.geom", "Normal.frag");
    wireframeShader.setup(shaderPath, "Wireframe.vert", "", "Wireframe.frag");
    pointShader.setup(shaderPath, "Point.vert", "", "Point.frag");
}

void setupUniformBlocks()
{
    // 1) generate transform indices
    GLuint meshShaderIndex = glGetUniformBlockIndex(meshShader.program, "Transform");
    GLuint normalShaderIndex = glGetUniformBlockIndex(normalShader.program, "Transform");
    GLuint wireframeShaderIndex = glGetUniformBlockIndex(wireframeShader.program, "Transform");
    GLuint pointShaderIndex = glGetUniformBlockIndex(pointShader.program, "Transform");
    
    // bind
    glUniformBlockBinding(meshShader.program, meshShaderIndex, 0);
    glUniformBlockBinding(normalShader.program, normalShaderIndex, 0);
    glUniformBlockBinding(wireframeShader.program, wireframeShaderIndex, 0);
    glUniformBlockBinding(pointShader.program, pointShaderIndex, 0);
    
    // add transform data
    glGenBuffers(1, &transformUbo);
    glBindBuffer(GL_UNIFORM_BUFFER, transformUbo);
    glBindBufferBase(GL_UNIFORM_BUFFER, 0, transformUbo);
    glBufferData(GL_UNIFORM_BUFFER, 3*sizeof(Eigen::Matrix4f), NULL, GL_DYNAMIC_DRAW);
    glBindBuffer(GL_UNIFORM_BUFFER, 0);
    
    // 2) generate light index
    meshShaderIndex = glGetUniformBlockIndex(meshShader.program, "Light");
    
    // bind
    glUniformBlockBinding(meshShader.program, meshShaderIndex, 1);
    
    // add light data
    glGenBuffers(1, &lightUbo);
    glBindBuffer(GL_UNIFORM_BUFFER, lightUbo);
    glBufferData(GL_UNIFORM_BUFFER, 2*sizeof(Eigen::Vector4f), NULL, GL_STATIC_DRAW); // std140 alignment
    glBindBufferBase(GL_UNIFORM_BUFFER, 1, lightUbo);
    
    glBufferSubData(GL_UNIFORM_BUFFER, 0, sizeof(Eigen::Vector4f), lightPosition.data());
    glBufferSubData(GL_UNIFORM_BUFFER, sizeof(Eigen::Vector4f), sizeof(Eigen::Vector4f), lightColor.data());
    glBindBuffer(GL_UNIFORM_BUFFER, 0);
}

void printInstructions()
{
    std::cerr << "1: compute hks\n"
              << "2: compute wks\n"
              << "3: toggle descriptor\n"
              << "4: toggle feature points\n"
              << "5: generate patches\n"
              << "6: toggle normals\n"
              << "7: toggle wireframe\n"
              << "→/←: change descriptor level\n"
              << "w/s: move in/out\n"
              << "a/d: move left/right\n"
              << "e/q: move up/down\n"
              << "escape: exit program\n"
              << std::endl;
}

void setFeaturePoints()
{
    for (int i = 0; i < 2; i++) {
        for (VertexCIter v = meshes[i].vertices.begin(); v != meshes[i].vertices.end(); v++) {
            if (v->isFeature(t)) featureMaps[i][v->index] = true;
            else featureMaps[i][v->index] = false;
        }
    }
}

void setColor(int i, bool useDescriptor = false)
{
    colors[i] = std::vector<Eigen::Vector3f>(meshes[i].vertices.size());
    for (VertexCIter v = meshes[i].vertices.begin(); v != meshes[i].vertices.end(); v++) {
        if (useDescriptor) colors[i][v->index] = Eigen::Vector3f(v->descriptor(t), 0.0, 0.0);
        else colors[i][v->index] = defaultColor;
    }
}

void updateColor()
{
    for (int i = 0; i < (int)meshes.size(); i++) {
        if (showDescriptor) setColor(i, true);
        else setColor(i);
        glMeshes[i].update(colors[i]);
    }
}

void generatePatchColors(std::vector<Eigen::Vector3f>& patchColors)
{
    for (int i = 0; i < (int)patchColors.size(); i++) {
        patchColors[i].x() = (double)rand() / RAND_MAX;
        patchColors[i].y() = (double)rand() / RAND_MAX;
        patchColors[i].z() = (double)rand() / RAND_MAX;
    }
}

bool terminate(const std::vector<std::queue<VertexCIter>>& queues)
{
    bool allEmpty = true;
    for (int q = 0; q < (int)queues.size(); q++) {
        if (!queues[q].empty()) allEmpty = false;
    }
    
    return allEmpty;
}

void setPatchColors()
{
    for (int m = 0; m < 2; m++) {
        int p = 0;
        for (VertexCIter v = meshes[m].vertices.begin(); v != meshes[m].vertices.end(); v++) {
            if (featureMaps[m][v->index]) p++;
        }
        
        std::unordered_map<int, bool> visited;
        std::vector<std::queue<VertexCIter>> queues(p);
        std::vector<int> levels(p, 0);
        std::vector<Eigen::Vector3f> patchColors(p);
        generatePatchColors(patchColors);
        VertexCIter end = meshes[m].vertices.end();
        
        // initialize 
        int i = 0;
        for (VertexCIter v = meshes[m].vertices.begin(); v != meshes[m].vertices.end(); v++) {
            if (featureMaps[m][v->index]) {
                visited[v->index] = true;
                colors[m][v->index] = patchColors[i];
                queues[i].push(v);
                queues[i].push(end);
                i++;
            }
        }
        
        // perform bfs
        i = 0;
        while (!terminate(queues)) {
            while (queues[i].empty()) i = (i+1) % p;
            VertexCIter v = queues[i].front();
            queues[i].pop();
            
            if (v == end) {
                levels[i]++;
                queues[i].push(end);
                if (queues[i].front() == end) queues[i].pop();
                i = (i+1) % p;
                
            } else {
                HalfEdgeCIter h = v->he;
                do {
                    VertexCIter vn = h->flip->vertex;
                    if (!visited[vn->index]) {
                        colors[m][vn->index] = patchColors[i];
                        
                        queues[i].push(vn);
                        visited[vn->index] = true;
                    }
                    
                    h = h->flip->next;
                } while (h != v->he);
            }
        }
        
        // update
        glMeshes[m].update(colors[m]);
    }
}

void init()
{
    // enable depth testing and multisampling
    glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
    glEnable(GL_DEPTH_TEST);
    glEnable(GL_MULTISAMPLE);
    
    // setup shaders and blocks
    setupShaders();
    setupUniformBlocks();
    
    // read meshes
    success = true;
    for (int i = 0; i < (int)meshes.size(); i++) {
        if (!meshes[i].read(paths[i])) success = false;
        if (success) {
            setColor(i);
            glMeshes[i].setup(colors[i]);
            
            for (VertexCIter v = meshes[i].vertices.begin(); v != meshes[i].vertices.end(); v++) {
                glPoints[i].push_back(GLPoint(v->position.cast<float>()));
                glPoints[i][v->index].setup();
            }
        }
    }
    
    // print instructions
    printInstructions();
}

void uploadCameraTransforms()
{
    // set camera transformations
    glBindBuffer(GL_UNIFORM_BUFFER, transformUbo);
    glBufferSubData(GL_UNIFORM_BUFFER, 0, sizeof(glm::mat4),
                    glm::value_ptr(camera.projectionMatrix(gridX, gridY)));
    glBufferSubData(GL_UNIFORM_BUFFER, sizeof(glm::mat4), sizeof(glm::mat4),
                    glm::value_ptr(camera.viewMatrix()));
    glBindBuffer(GL_UNIFORM_BUFFER, 0);
    
    // set view position
    meshShader.use();
    glUniform3f(glGetUniformLocation(meshShader.program, "viewPosition"),
                camera.pos.x(), camera.pos.y(), camera.pos.z());
}

void uploadMeshTransform(const Eigen::Matrix4f& transform)
{
    // set transform
    glBindBuffer(GL_UNIFORM_BUFFER, transformUbo);
    glBufferSubData(GL_UNIFORM_BUFFER, 2*sizeof(Eigen::Matrix4f), sizeof(Eigen::Matrix4f),
                    transform.data());
    glBindBuffer(GL_UNIFORM_BUFFER, 0);
}

void draw(int m)
{
    // draw mesh
    glMeshes[m].draw(meshShader);
    if (showNormals) glMeshes[m].draw(normalShader);
    if (showWireframe) {
        glEnable(GL_BLEND);
        glBlendFunc(GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);
        glDepthMask(GL_FALSE);
        glMeshes[m].draw(wireframeShader);
        glDepthMask(GL_TRUE);
        glDisable(GL_BLEND);
    }
    
    // draw points
    if (computedDescriptor && showFeaturePoints) {
        for (VertexCIter v = meshes[m].vertices.begin(); v != meshes[m].vertices.end(); v++) {
            if (featureMaps[m][v->index]) {
                glPoints[m][v->index].draw(pointShader);
            }
        }
    }
}

void display()
{
    float elapsedTime = glutGet(GLUT_ELAPSED_TIME);
    dt = (elapsedTime - lastTime) / 1000.0;
    lastTime = elapsedTime;
    
    if (success) {
        // upload camera transforms
        uploadCameraTransforms();
        
        // clear
        glClearColor(0.1, 0.1, 0.1, 1.0);
        glClearDepth(1.0);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
        glDepthFunc(GL_LEQUAL);
        
        // draw
        Eigen::Matrix4f transform = Eigen::Matrix4f::Identity(); transform(0, 3) = -1.0;
        uploadMeshTransform(transform);
        draw(0);
        
        transform(0, 3) = 1.0;
        uploadMeshTransform(transform);
        draw(1);
        
        // swap
        glutSwapBuffers();
    }
}

void idle()
{
    glutPostRedisplay();
}

void reset()
{
    for (int i = 0; i < (int)glMeshes.size(); i++) {
        glMeshes[i].reset();
        for (VertexCIter v = meshes[i].vertices.begin(); v != meshes[i].vertices.end(); v++) {
            glPoints[i][v->index].reset();
        }
    }
    
    meshShader.reset();
    normalShader.reset();
    wireframeShader.reset();
    pointShader.reset();
    glDeleteBuffers(1, &transformUbo);
    glDeleteBuffers(1, &lightUbo);
}

void keyboardPressed(unsigned char key, int x, int y)
{
    keys[key] = true;
    
    if (keys[ESCAPE]) {
        reset();
        exit(0);
        
    } else if (keys[DIGIT_OFFSET + 1]) {
        for (int i = 0; i < (int)meshes.size(); i++) meshes[i].computeDescriptor(HKS);
        computedDescriptor = true;
        setFeaturePoints();
        if (showDescriptor) updateColor();
        
    } else if (keys[DIGIT_OFFSET + 2]) {
        for (int i = 0; i < (int)meshes.size(); i++) meshes[i].computeDescriptor(WKS);
        computedDescriptor = true;
        setFeaturePoints();
        if (showDescriptor) updateColor();
        
    } else if (keys[DIGIT_OFFSET + 3]) {
        showDescriptor = !showDescriptor;
        if (showDescriptor && !computedDescriptor) showDescriptor = false;
        else {
            updateColor();
            
            std::string title = "Mesh Correspondence";
            if (showDescriptor) title += ", t: " + std::to_string(t);
            glutSetWindowTitle(title.c_str());
        }
        
    } else if (keys[DIGIT_OFFSET + 4]) {
        showFeaturePoints = !showFeaturePoints;
        
    } else if (keys[DIGIT_OFFSET + 5]) {
        if (computedDescriptor) setPatchColors();
        
    } else if (keys[DIGIT_OFFSET + 6]) {
        showNormals = !showNormals;
        
    } else if (keys[DIGIT_OFFSET + 7]) {
        showWireframe = !showWireframe;
        
    } else if (keys['a']) {
        camera.processKeyboard(LEFT, dt);
        
    } else if (keys['d']) {
        camera.processKeyboard(RIGHT, dt);
        
    } else if (keys['w']) {
        camera.processKeyboard(FORWARD, dt);
        
    } else if (keys['s']) {
        camera.processKeyboard(BACKWARD, dt);
        
    } else if (keys['e']) {
        camera.processKeyboard(UP, dt);
        
    } else if (keys['q']) {
        camera.processKeyboard(DOWN, dt);
    }
}

void keyboardReleased(unsigned char key, int x, int y)
{
    if (key != ESCAPE) keys[key] = false;
}

void special(int i, int x, int y)
{
    switch (i) {
        case GLUT_KEY_LEFT:
            t--;
            if (t < 0) t = (int)meshes[0].vertices[0].descriptor.size() - 1;
            break;
        case GLUT_KEY_RIGHT:
            t++;
            int n = std::max(0, (int)meshes[0].vertices[0].descriptor.size() - 1);
            if (t > n) t = 0;
            break;
    }
    if (computedDescriptor) setFeaturePoints();
    if (showDescriptor) updateColor();
    
    std::string title = "Mesh Correspondence, t: " + std::to_string(t);
    glutSetWindowTitle(title.c_str());
}

void mouse(int x, int y)
{
    if (firstMouse) {
        lastX = x;
        lastY = y;
        firstMouse = false;
    }
    
    float dx = x - lastX;
    float dy = lastY - y;
    
    lastX = x;
    lastY = y;
    
    camera.processMouse(dx, dy);
}

int main(int argc, char** argv)
{
    if (argc != 4) {
        std::cout << "Usage: " << argv[0] << " OBJ_PATH_1 OBJ_PATH_2 SHADER_PATH" << std::endl;
        return 0;
    }
    
    paths.push_back(argv[1]);
    paths.push_back(argv[2]);
    shaderPath = argv[3];
    
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE | GLUT_DEPTH | GLUT_3_2_CORE_PROFILE | GLUT_MULTISAMPLE);
    glutInitWindowSize(gridX, gridY);
    glutCreateWindow("Mesh Correspondence");
    
    init();
    
    glutDisplayFunc(display);
    glutIdleFunc(idle);
    glutKeyboardFunc(keyboardPressed);
    glutKeyboardUpFunc(keyboardReleased);
    glutSpecialFunc(special);
    glutMotionFunc(mouse);
    glutMainLoop();
    
    return 0;
}
