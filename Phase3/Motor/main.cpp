#define _CRT_SECURE_NO_WARNINGS

#include "./tinyxml/tinyxml.h"
#include "./tinyxml/tinystr.h"
#include "Vertex.h"
#include "Group.h"
#include <vector>
#include <fstream>
#include <iostream>
#include <string>
#include <GL/glew.h>
#include <GL/glut.h>
#include <math.h>
#define _USE_MATH_DEFINES
#include <math.h>
#include <ctype.h>
#include <map>

using namespace std;

list<Group> groups;

float height = 2.0f;
float x = 0.0f;
float z = 0.0f;


float camX = 0, camY = 202.379883, camZ = 525.505249;
int startX,startY, tracking = 0;

int alpha = 0, beta = 22, r = 560;

bool orientation, points;

Group initGroup()
{
	Translate t;
	t.time = t.x = t.y = t.z = 0;

	Rotate r;
	r.time = 0;
	r.angle = 0;
	r.axisX = 0;
	r.axisY = 0;
	r.axisZ = 0;

	Scale s;
	s.x = s.y = s.z = 1;

	Colour c;
	c.r = c.b = c.g = 1;

	Group g;
	g.orientation = false;
	g.points = false;
	g.t = t;
	g.r = r;
	g.s = s;
	g.c = c;
	g.orbitRadius = 0.0;
	return g;
}


void loadGroups(TiXmlNode *currentranslateNode, list<Group>* groups)
{
	for (TiXmlNode* groupNode = currentranslateNode->FirstChild("group"); groupNode; groupNode = groupNode->NextSiblingElement())
	{
		Group g = initGroup();

		TiXmlNode* translateNode = groupNode->FirstChild("translate");
		if (translateNode != nullptr)
		{
			TiXmlElement* tElem = translateNode->ToElement();
			const char *time,*tx, *ty, *tz;
			time = tElem->Attribute("time");
			tx = tElem->Attribute("X");
			ty = tElem->Attribute("Y");
			tz = tElem->Attribute("Z");

			if (time)
				g.t.time = atof(time);
			if (tx)
				g.t.x = atof(tx);
			if (ty)
				g.t.y = atof(ty);
			if (tz)
				g.t.z = atof(tz);

			TiXmlNode* pointNode = tElem->FirstChild("point");
			if (pointNode != nullptr)
			{
				for (TiXmlElement* pointElem = tElem->FirstChild("point")->ToElement(); pointElem; pointElem = pointElem->NextSiblingElement())
				{
					vertex v;
					v.x = 0;
					v.y = 0;
					v.z = 0;
					const char *x, *y, *z;
					x = pointElem->Attribute("X");
					y = pointElem->Attribute("Y");
					z = pointElem->Attribute("Z");

					if (x)
						v.x = atof(x);
					if (y)
						v.y = atof(y);
					if (z)
						v.z = atof(z);

					g.t.controlPoints.push_back(v);

					
				}
			}
		}

		TiXmlNode* rotateNode = groupNode->FirstChild("rotate");
		if (rotateNode != nullptr)
		{
			TiXmlElement* rElem = rotateNode->ToElement();
			const char *angle, *time, *ax, *ay, *az;
			angle = rElem->Attribute("angle");
			time = rElem->Attribute("time");
			ax = rElem->Attribute("axisX");
			ay = rElem->Attribute("axisY");
			az = rElem->Attribute("axisZ");

			if (time)
				g.r.time = atof(time);
			if (angle)
				g.r.angle = atof(angle);
			if (ax)
				g.r.axisX = atoi(ax);
			if (ay)
				g.r.axisY = atoi(ay);
			if (az)
				g.r.axisZ = atoi(az);
		}

		TiXmlNode* scaleNode = groupNode->FirstChild("scale");
		if (scaleNode != nullptr)
		{
			TiXmlElement* sElem = scaleNode->ToElement();
			const char *ax, *ay, *az;
			ax = sElem->Attribute("X");
			ay = sElem->Attribute("Y");
			az = sElem->Attribute("Z");

			if (ax)
				g.s.x = atof(ax);
			if (ay)
				g.s.y = atof(ay);
			if (az)
				g.s.z = atof(az);
		}


		TiXmlNode* colourNode = groupNode->FirstChild("colour");
		if (colourNode != nullptr)
		{
			TiXmlElement* cElem = colourNode->ToElement();
			const char *R, *G, *B;
			R = cElem->Attribute("R");
			G = cElem->Attribute("G");
			B = cElem->Attribute("B");

			if (R)
				g.c.r = atof(R);
			if (G)
				g.c.g = atof(G);
			if (B)
				g.c.b = atof(B);
		}
		TiXmlNode* orbitNode = groupNode->FirstChild("orbit");
		if (orbitNode != nullptr)
		{
			TiXmlElement* oElem = orbitNode->ToElement();
			const char *r;
			r = oElem->Attribute("radius");
			if (r)
				g.orbitRadius = atof(r);
		}

		TiXmlNode* orientationNode = groupNode->FirstChild("orientation");
		if (orientationNode != nullptr)
		{
			TiXmlElement* orElem = orientationNode->ToElement();
			const char *f;
			int o;
			f = orElem->Attribute("FLAG");
			if (f)
				o = atof(f);
			if (o == 0)
				g.orientation = true;
		}

		TiXmlNode* pointNode = groupNode->FirstChild("points");
		if (pointNode != nullptr)
		{
			TiXmlElement* poElem = pointNode->ToElement();
			const char *f;
			int o;
			f = poElem->Attribute("FLAG");
			if (f)
				o = atof(f);
			if (o == 0)
				g.points = true;
		}

		int nbuffers = 0;
		TiXmlNode* modelsNode = groupNode->FirstChild("models");
		if (modelsNode != nullptr)
		{
			
			TiXmlElement* mElem = modelsNode->ToElement();
			TiXmlNode* modelNode = mElem->FirstChild("model");

			if (modelNode != nullptr)
			{
				for (TiXmlElement* modElem = mElem->FirstChild("model")->ToElement(); modElem; modElem = modElem->NextSiblingElement())
				{
					const char *file = modElem->Attribute("file");
					string name = "";
					name += file;

					string filename("..\\..\\Solids\\");
					filename += file;
					fstream f;
					f.open(filename);

					glEnableClientState(GL_VERTEX_ARRAY);

					if (f.is_open())
					{
						vector <float> solid;
						string line;
						int numVertices;
						getline(f, line);
						numVertices = stoi(line.c_str());

						while (getline(f, line))
						{
							float x, y, z;
							sscanf(line.c_str(), "%f %f %f\n", &x, &y, &z);

							solid.push_back(x);
							solid.push_back(y);
							solid.push_back(z);
							
						}
						f.close();

						Model m;
						m.name = name;
						m.numberPoints[nbuffers] = solid.size();
						glGenBuffers(1, m.vertices);
						glBindBuffer(GL_ARRAY_BUFFER, m.vertices[nbuffers++]);
						glBufferData(GL_ARRAY_BUFFER, solid.size() * sizeof(float), &solid[0], GL_STATIC_DRAW);

						
						g.m.push_back(m);

					}
				}
			}
		}
		if (groupNode != nullptr)
			loadGroups(groupNode, &g.g);

		groups->push_back(g);
	}
}

list<Group> load(const char* pFilename)
{
	list<Group> groupsI;
	TiXmlDocument xmlDoc(pFilename);
	if (xmlDoc.LoadFile())
	{
		TiXmlNode* root = TiXmlHandle(xmlDoc.RootElement()).ToNode();
		loadGroups(root, &groupsI);
	}
	return groupsI;

}

void changeSize(int w, int h) {

	// Prevent a divide by zero, when window is too short
	// (you cant make a window with zero width).
	if (h == 0)
		h = 1;

	// compute window's aspect ratio 
	float ratio = w * 1.0 / h;

	// Reset the coordinate system before modifying
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();

	// Set the viewport to be the entire window
	glViewport(0, 0, w, h);

	// Set the correct perspective
	gluPerspective(45, ratio, 1, 1000);

	// return to the model view matrix mode
	glMatrixMode(GL_MODELVIEW);
}

void drawTriangles(GLuint solid, int numberPoints)
{
	glBindBuffer(GL_ARRAY_BUFFER, solid);
	glVertexPointer(3, GL_FLOAT, 0, 0);


	if (points == false)
		glDrawArrays(GL_TRIANGLES, 0, numberPoints);
	else
		glDrawArrays(GL_POINTS, 0, numberPoints/3);
}


void applyRotation(float time, int x, int y, int z)
{

	float angle = (360 * glutGet(GLUT_ELAPSED_TIME))/(100*time);
	glRotatef(angle, x, y, z);
}


void buildRotMatrix(float *x, float *y, float *z, float *m) {

	m[0] = x[0]; m[1] = x[1]; m[2] = x[2]; m[3] = 0;
	m[4] = y[0]; m[5] = y[1]; m[6] = y[2]; m[7] = 0;
	m[8] = z[0]; m[9] = z[1]; m[10] = z[2]; m[11] = 0;
	m[12] = 0; m[13] = 0; m[14] = 0; m[15] = 1;
}


void cross(float *a, float *b, float *res) {

	res[0] = a[1] * b[2] - a[2] * b[1];
	res[1] = a[2] * b[0] - a[0] * b[2];
	res[2] = a[0] * b[1] - a[1] * b[0];
}


void normalize(float *a) {

	float l = sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2]);
	a[0] = a[0] / l;
	a[1] = a[1] / l;
	a[2] = a[2] / l;
}

void multMatrixVector(float *m, float *v, float *res) {

	for (int j = 0; j < 4; ++j) {
		res[j] = 0;
		for (int k = 0; k < 4; ++k) {
			res[j] += v[k] * m[j * 4 + k];
		}
	}

}


void getCatmullRomPoint(float t, float *p0, float *p1, float *p2, float *p3, float *res, float *deriv) {

	// catmull-rom matrix
	float m[4][4] = { { -0.5f,  1.5f, -1.5f,  0.5f },
	{ 1.0f, -2.5f,  2.0f, -0.5f },
	{ -0.5f,  0.0f,  0.5f,  0.0f },
	{ 0.0f,  1.0f,  0.0f,  0.0f } };

	// reset res and deriv
	res[0] = 0.0; res[1] = 0.0; res[2] = 0.0;
	deriv[0] = 0.0; deriv[1] = 0.0; deriv[2] = 0.0;

	// Compute point A=M*P

	float pointX[4];
	pointX[0] = p0[0];
	pointX[1] = p1[0];
	pointX[2] = p2[0];
	pointX[3] = p3[0];

	float pointY[4];
	pointY[0] = p0[1];
	pointY[1] = p1[1];
	pointY[2] = p2[1];
	pointY[3] = p3[1];


	float pointZ[4];
	pointZ[0] = p0[2];
	pointZ[1] = p1[2];
	pointZ[2] = p2[2];
	pointZ[3] = p3[2];

	float *ax = new float[4];
	float *ay = new float[4];
	float *az = new float[4];

	multMatrixVector(*m, pointX, ax);
	multMatrixVector(*m, pointY, ay);
	multMatrixVector(*m, pointZ, az);

	// Compute point res = T *A

	res[0] = pow(t, 3) * ax[0] + pow(t, 2) * ax[1] + t * ax[2] + ax[3];
	res[1] = pow(t, 3) * ay[0] + pow(t, 2) * ay[1] + t * ay[2] + ay[3];
	res[2] = pow(t, 3) * az[0] + pow(t, 2) * az[1] + t * az[2] + az[3];

	// compute deriv = T' * A



	deriv[0] = 3 * pow(t,2) * ax[0] + 2 * t * ax[1] +  ax[2] ;
	deriv[1] = 3 * pow(t, 2) * ay[0] + 2 * t * ay[1] +  ay[2] ;
	deriv[2] = 3 * pow(t, 2) * az[0] + 2 * t * az[1] +  az[2] ;



}

// given  global t, returns the point in the curve
void getGlobalCatmullRomPoint(float gt, float *res, float *deriv, vector<vertex> cPoints) {

	const int size = (int)cPoints.size();

	float t = gt * size; // this is the real global t
	int index = floor(t);  // which segment
	t = t - index; // where within  the segment

	

	// indices store the points
	int indices[4];
	indices[0] = (index + size - 1) % size;
	indices[1] = (indices[0] + 1) % size;
	indices[2] = (indices[1] + 1) % size;
	indices[3] = (indices[2] + 1) % size;


	float** p = new float*[size];
	for (int i = 0; i < size; i++)
		p[i] = new float[3];

	for (int i = 0; i < size; i++)
	{
		
		p[i][0] = cPoints[i].x;
		p[i][1] = cPoints[i].y;
		p[i][2] = cPoints[i].z;
	}
	
	getCatmullRomPoint(t, p[indices[0]], p[indices[1]], p[indices[2]], p[indices[3]], res, deriv);
}



void applyTranslation(vector<vertex> cPoints, float time)
{
	float res[3], div[3];
	

	int deltaTime = glutGet(GLUT_ELAPSED_TIME);
	
	float t = deltaTime / (1000 * time);

	getGlobalCatmullRomPoint(t, res, div, cPoints);
	glTranslatef(res[0], res[1], res[2]);
	
	if (orientation == true)
	{
		float up[3] = { 0, 1, 0 };
		float r[3];

		cross(div, up, r);
		normalize(r);

		cross(r, div, up);

		normalize(up);
		normalize(div);

		float m[4][4] =
		{
			{ div[0], up[0], r[0], res[0] },
			{ div[1], up[1], r[1], res[1] },
			{ div[2], up[2], r[2], res[2] },
			{ 0, 0, 0, 1 }
		};

		buildRotMatrix(div, up, r, *m);

		glMultMatrixf(*m);
	}

}

void drawGroupElements(list<Group> g)
{
	glPointSize(1.1);

	for (list<Group>::iterator itg = g.begin(); itg != g.end(); itg++)
	{
		glPushMatrix();

		orientation = itg->orientation;

		points = itg->points;

		glTranslatef(itg->t.x, itg->t.y, itg->t.z);
		
		glScalef(itg->s.x, itg->s.y, itg->s.z);
		glColor3f(itg->c.r, itg->c.g, itg->c.b);

		if (itg->r.time > 0)
			applyRotation(itg->r.time, itg->r.axisX, itg->r.axisY, itg->r.axisZ);
		
		if (itg->t.time > 0)
			applyTranslation(itg->t.controlPoints, itg->t.time);


		if(itg->r.angle > 0)
			glRotatef(itg->r.angle, itg->r.axisX, itg->r.axisY, itg->r.axisZ);

		int nbuffers = 0;

		for (list<Model>::iterator itm = itg->m.begin(); itm != itg->m.end(); itm++)
		{

			drawTriangles(itm->vertices[nbuffers], itm->numberPoints[nbuffers]);

			nbuffers++;
		}
		drawGroupElements(itg->g);
		glPopMatrix();
	}
}

/*
Função que desenha as primitivas gráficas.
*/
void renderScene(void)
{

	float pos[4] = { -1.0, 1.0, 1.0, 0.0 };

	glClearColor(0.0f, 0.0f, 0.0f, 0.0f);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	glLoadIdentity();
	gluLookAt(camX, camY, camZ,
		0.0, 0.0, 0.0,
		0.0f, 1.0f, 0.0f);


	glEnable(GL_CULL_FACE);
	glFrontFace(GL_CCW);

	glPolygonMode(GL_FRONT, GL_LINE);

	drawGroupElements(groups);


	glutSwapBuffers();
}


void processMouseButtons(int button, int state, int xx, int yy) {

	if (state == GLUT_DOWN) {
		startX = xx;
		startY = yy;
		if (button == GLUT_LEFT_BUTTON)
			tracking = 1;
		else if (button == GLUT_RIGHT_BUTTON)
			tracking = 2;
		else
			tracking = 0;
	}
	else if (state == GLUT_UP) {
		if (tracking == 1) {
			alpha += (xx - startX);
			beta += (yy - startY);
		}
		else if (tracking == 2) {

			r -= yy - startY;
			if (r < 3)
				r = 3.0;
		}
		tracking = 0;
	}
}


void processMouseMotion(int xx, int yy) {

	int deltaX, deltaY;
	int alphaAux, betaAux;
	int rAux;

	if (!tracking)
		return;

	deltaX = xx - startX;
	deltaY = yy - startY;

	if (tracking == 1) {


		alphaAux = alpha + deltaX;
		betaAux = beta + deltaY;

		if (betaAux > 85.0)
			betaAux = 85.0;
		else if (betaAux < -85.0)
			betaAux = -85.0;

		rAux = r;
	}
	else if (tracking == 2) {

		alphaAux = alpha;
		betaAux = beta;
		rAux = r - deltaY;
		if (rAux < 3)
			rAux = 3;
	}
	camX = rAux * sin(alphaAux * 3.14 / 180.0) * cos(betaAux * 3.14 / 180.0);
	camZ = rAux * cos(alphaAux * 3.14 / 180.0) * cos(betaAux * 3.14 / 180.0);
	camY = rAux * 							     sin(betaAux * 3.14 / 180.0);


}

/*
Função que inicia o programa e recebe como argumento o ficheiro XML que será lido.
*/
int main(int argc, char **argv) {

	orientation = false;
	
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGBA);
	glutInitWindowPosition(0, 0);
	glutInitWindowSize(1500, 700);
	glutCreateWindow("cg@di");

	glewInit();

	string filename("..\\..\\Solids\\");
	filename += argv[1];
	groups = load(filename.c_str());

	glutDisplayFunc(renderScene);
	glutIdleFunc(renderScene);
	glutReshapeFunc(changeSize);

	glutMouseFunc(processMouseButtons);
	glutMotionFunc(processMouseMotion);

	glEnable(GL_DEPTH_TEST);
	glEnable(GL_CULL_FACE);
	glClearColor(0.0f, 0.0f, 0.0f, 0.0f);

	glutMainLoop();

	return 1;
}

