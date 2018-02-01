#define _CRT_SECURE_NO_WARNINGS /*Fase 4 -> Motor*/

#include "./tinyxml/tinyxml.h"
#include "./tinyxml/tinystr.h"
#include "Vertex.h"
#include "Group.h"
#include "Lights.h"
#include <vector>
#include <fstream>
#include <iostream>
#include <string>
#include <GL/glew.h>
#include <GL/glut.h>
#include <IL/il.h>
#include <math.h>
#define _USE_MATH_DEFINES
#include <math.h>
#include <ctype.h>
#include <map>

using namespace std;

list<Group> groups;
vector<Light> lights;

float camX, camY, camZ;
int startX, startY, tracking = 0;

int alpha, beta, r;

bool orientation, points;

int numLights = 0;

int time, timebase, frame = 0, fps = 0;

void calculateCam()
{
	camX = r * sin(alpha * 3.14 / 180.0) * cos(beta * 3.14 / 180.0);
	camZ = r * cos(alpha * 3.14 / 180.0) * cos(beta * 3.14 / 180.0);
	camY = r * 							     sin(beta * 3.14 / 180.0);
}

void loadCamera(TiXmlNode *currentranslateNode)
{

	TiXmlNode* cameraNode = currentranslateNode->FirstChild("camera");
	if (cameraNode != nullptr)
	{

		TiXmlElement* cElem = cameraNode->ToElement();
		TiXmlNode* posNode = cElem->FirstChild("position");

		if (posNode != nullptr)
		{
			for (TiXmlElement* camElem = cElem->FirstChild("position")->ToElement(); camElem; camElem = camElem->NextSiblingElement())
			{

				const char *alph, *bet, *rad;

				alph = camElem->Attribute("alpha");
				bet = camElem->Attribute("beta");
				rad = camElem->Attribute("radius");


				if (alph)
					alpha = atof(alph);
				if (bet)
					beta = atof(bet);
				if (rad)
					r = atof(rad);

				calculateCam();
			}
		}
	}
}

void initLights()
{
	Light l;
	l.pos[0] = 0;
	l.pos[1] = 0;
	l.pos[2] = 0;
	l.pos[3] = 0;
	l.amb[0] = 0;
	l.amb[1] = 0;
	l.amb[2] = 0;
	l.amb[3] = 0;
	l.dif[0] = 0;
	l.dif[1] = 0;
	l.dif[2] = 0;
	l.dif[3] = 0;
	l.type = "";
	l.spotDir[0] = 0;
	l.spotDir[1] = 0;
	l.spotDir[2] = 0;
	l.angle = 0;
	l.exp = 0;
	for (int i = 0; i < 8; i++)
	{
		lights.push_back(l);
	}
}

void loadLights(TiXmlNode *currentranslateNode)
{
	initLights();

	TiXmlNode* lightsNode = currentranslateNode->FirstChild("lights");
	if (lightsNode != nullptr)
	{

		TiXmlElement* lElem = lightsNode->ToElement();
		TiXmlNode* lightNode = lElem->FirstChild("light");

		if (lightNode != nullptr)
		{
			for (TiXmlElement* lightElem = lElem->FirstChild("light")->ToElement(); lightElem; lightElem = lightElem->NextSiblingElement())
			{

				const char *type, *x, *y, *z, *dirx, *diry, *dirz, *angle, *exp,
					*ambR, *ambG, *ambB, *difR, *difG, *difB;


				type = lightElem->Attribute("type");
				x = lightElem->Attribute("posX");
				y = lightElem->Attribute("posY");
				z = lightElem->Attribute("posZ");
				dirx = lightElem->Attribute("dirX");
				diry = lightElem->Attribute("dirY");
				dirz = lightElem->Attribute("dirZ");
				angle = lightElem->Attribute("angle");
				exp = lightElem->Attribute("exp");
				ambR = lightElem->Attribute("ambR");
				ambG = lightElem->Attribute("ambG");
				ambB = lightElem->Attribute("ambB");
				difR = lightElem->Attribute("difR");
				difG = lightElem->Attribute("difG");
				difB = lightElem->Attribute("difB");

				if (ambB)
					lights[numLights].amb[2] = atof(ambB);
				if (ambR)
					lights[numLights].amb[0] = atof(ambR);
				if (ambG)
					lights[numLights].amb[1] = atof(ambG);
				if (difB)
					lights[numLights].dif[2] = atof(difB);
				if (difR)
					lights[numLights].dif[0] = atof(difR);
				if (difG)
					lights[numLights].dif[1] = atof(difG);
				if (x)
					lights[numLights].pos[0] = atof(x);
				if (y)
					lights[numLights].pos[1] = atof(y);
				if (z)
					lights[numLights].pos[2] = atof(z);
				if (dirx)
					lights[numLights].spotDir[0] = atof(dirx);
				if (diry)
					lights[numLights].spotDir[1] = atof(diry);
				if (dirz)
					lights[numLights].spotDir[2] = atof(dirz);
				if (angle)
					lights[numLights].angle = atof(angle);
				if (exp)
					lights[numLights].exp = atof(exp);

				string type1("POINT");
				string type2("DIRECTIONAL");
				string type3("SPOTLIGHT");

				if (type)
				{
					lights[numLights].type = "";
					lights[numLights].type += type;
					if (type1.compare(type) == 0 || type3.compare(type) == 0)
						lights[numLights].pos[3] = 1.0f;
					if (type2.compare(type) == 0)
						lights[numLights].pos[3] = 0.0f;
				}

				numLights++;
			}
		}
	}
}

int loadTexture(std::string s) {

	unsigned int t, tw, th;
	unsigned char *texData;
	unsigned int texID;

	ilInit();
	ilEnable(IL_ORIGIN_SET);
	ilOriginFunc(IL_ORIGIN_UPPER_LEFT);
	ilGenImages(1, &t);
	ilBindImage(t);
	ilLoadImage((ILstring)s.c_str());
	tw = ilGetInteger(IL_IMAGE_WIDTH);
	th = ilGetInteger(IL_IMAGE_HEIGHT);
	ilConvertImage(IL_RGBA, IL_UNSIGNED_BYTE);
	texData = ilGetData();

	glGenTextures(1, &texID);

	glBindTexture(GL_TEXTURE_2D, texID);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);

	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR_MIPMAP_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);

	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, tw, th, 0, GL_RGBA, GL_UNSIGNED_BYTE, texData);
	glGenerateMipmap(GL_TEXTURE_2D);

	glBindTexture(GL_TEXTURE_2D, 0);

	return texID;

}

Model initModel(int nbuffers)
{
	Model m;
	m.lc[nbuffers].ambR = 0;
	m.lc[nbuffers].ambG = 0;
	m.lc[nbuffers].ambB = 0;
	m.lc[nbuffers].difR = 0;
	m.lc[nbuffers].difG = 0;
	m.lc[nbuffers].difB = 0;
	m.lc[nbuffers].emiR = 0;
	m.lc[nbuffers].emiG = 0;
	m.lc[nbuffers].emiB = 0;
	m.texture = "";
	m.texID = 0;
	return m;
}

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
	r.rotationAngle = 360;

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
			const char *time, *tx, *ty, *tz;
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
			const char *rotationAngle, *angle, *time, *ax, *ay, *az;
			rotationAngle = rElem->Attribute("rotationAngle");
			angle = rElem->Attribute("angle");
			time = rElem->Attribute("time");
			ax = rElem->Attribute("axisX");
			ay = rElem->Attribute("axisY");
			az = rElem->Attribute("axisZ");

			if (rotationAngle)
				g.r.rotationAngle = atof(rotationAngle);
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
					Model m = initModel(nbuffers);

					const char *file = modElem->Attribute("file");
					const char *texture = modElem->Attribute("texture");
					const char *ax, *ay, *az, *dx, *dy, *dz, *ex, *ey, *ez;

					ax = modElem->Attribute("ambR");
					ay = modElem->Attribute("ambG");
					az = modElem->Attribute("ambB");
					dx = modElem->Attribute("difR");
					dy = modElem->Attribute("difG");
					dz = modElem->Attribute("difB");
					ex = modElem->Attribute("emiR");
					ey = modElem->Attribute("emiG");
					ez = modElem->Attribute("emiB");

					if (ax)
						m.lc[nbuffers].ambR = atof(ax);
					else m.lc[nbuffers].ambR = 0;
					if (ay)
						m.lc[nbuffers].ambG = atof(ay);
					else m.lc[nbuffers].ambG = 0;
					if (az)
						m.lc[nbuffers].ambB = atof(az);
					else m.lc[nbuffers].ambB = 0;

					if (dx)
						m.lc[nbuffers].difR = atof(dx);
					else m.lc[nbuffers].difR = 0;
					if (dy)
						m.lc[nbuffers].difG = atof(dy);
					else m.lc[nbuffers].difG = 0;
					if (dz)
						m.lc[nbuffers].difB = atof(dz);
					else m.lc[nbuffers].difB = 0;

					if (ex)
						m.lc[nbuffers].emiR = atof(ex);
					else m.lc[nbuffers].emiR = 0;
					if (ey)
						m.lc[nbuffers].emiG = atof(ey);
					else m.lc[nbuffers].emiG = 0;
					if (ez)
						m.lc[nbuffers].emiB = atof(ez);
					else m.lc[nbuffers].emiB = 0;

					if (texture)
						m.texture = texture;
					else
						m.texture = "";

					string name = "";
					name += file;

					string filename("..\\..\\Solids\\");
					filename += file;
					fstream f;
					f.open(filename);


					glEnableClientState(GL_VERTEX_ARRAY);
					glEnableClientState(GL_NORMAL_ARRAY);
					glEnableClientState(GL_TEXTURE_COORD_ARRAY);

					if (f.is_open())
					{
						vector <float> solid;
						vector<float> normal;
						vector<float> tex;
						string line;
						int numVertices;
						getline(f, line);
						numVertices = stoi(line.c_str());

						while (getline(f, line))
						{
							float x, y, z, nx, ny, nz, s, t;
							sscanf(line.c_str(), "%f %f %f %f %f %f %f %f\n", &x, &y, &z, &nx, &ny, &nz, &s, &t);

							solid.push_back(x);
							solid.push_back(y);
							solid.push_back(z);

							normal.push_back(nx);
							normal.push_back(ny);
							normal.push_back(nz);

							tex.push_back(s);
							tex.push_back(t);
						}
						f.close();

						m.numberPoints[nbuffers] = solid.size();

						glGenBuffers(3, m.buffers[nbuffers].b);


						glBindBuffer(GL_ARRAY_BUFFER, m.buffers[nbuffers].b[0]);
						glBufferData(GL_ARRAY_BUFFER, solid.size() * sizeof(float), &solid[0], GL_STATIC_DRAW);

						glBindBuffer(GL_ARRAY_BUFFER, m.buffers[nbuffers].b[1]);
						glBufferData(GL_ARRAY_BUFFER, normal.size() * sizeof(float), &normal[0], GL_STATIC_DRAW);

						glBindBuffer(GL_ARRAY_BUFFER, m.buffers[nbuffers++].b[2]);
						glBufferData(GL_ARRAY_BUFFER, tex.size() * sizeof(float), &tex[0], GL_STATIC_DRAW);

						string texName("..\\..\\Textures\\");
						texName += m.texture;

						int texID;
						if (texName.compare("..\\..\\Textures\\") != 0)
							texID = loadTexture(texName);
						else
							texID = 0;
						m.texID = texID;

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
		TiXmlNode *root, *root1, *root2;
		root = TiXmlHandle(xmlDoc.RootElement()).ToNode();
		root1 = TiXmlHandle(xmlDoc.RootElement()).ToNode();
		root2 = TiXmlHandle(xmlDoc.RootElement()).ToNode();
		loadGroups(root, &groupsI);

		loadLights(root1);

		loadCamera(root2);
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



void drawTriangles(GLuint solid, GLuint normal, GLuint tex, int numberPoints)
{
	glBindBuffer(GL_ARRAY_BUFFER, solid);
	glVertexPointer(3, GL_FLOAT, 0, 0);

	glBindBuffer(GL_ARRAY_BUFFER, normal);
	glNormalPointer(GL_FLOAT, 0, 0);

	glBindBuffer(GL_ARRAY_BUFFER, tex);
	glTexCoordPointer(2, GL_FLOAT, 0, 0);

	if (points == false)
		glDrawArrays(GL_TRIANGLES, 0, numberPoints);
	else
		glDrawArrays(GL_POINTS, 0, numberPoints / 3);
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



	deriv[0] = 3 * pow(t, 2) * ax[0] + 2 * t * ax[1] + ax[2];
	deriv[1] = 3 * pow(t, 2) * ay[0] + 2 * t * ay[1] + ay[2];
	deriv[2] = 3 * pow(t, 2) * az[0] + 2 * t * az[1] + az[2];



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

void applyRotation(float time, int x, int y, int z)
{
	float angle = (360 * glutGet(GLUT_ELAPSED_TIME)) / (100 * time);
	glRotatef(angle, x, y, z);
}

int flag = 1;
void applyRotation(float time, int ang, int x, int y, int z)
{
	float angle = (ang * glutGet(GLUT_ELAPSED_TIME)) / (100 * time);
	int an = (int)angle % (ang * 2);
	if (an > ang)
		angle = -an;
	glRotatef(angle, x, y, z);
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
		{
			if (itg->r.rotationAngle == 360)
				applyRotation(itg->r.time, itg->r.axisX, itg->r.axisY, itg->r.axisZ);
			else
				applyRotation(itg->r.time, itg->r.rotationAngle, itg->r.axisX, itg->r.axisY, itg->r.axisZ);
		}

		if (itg->t.time > 0)
			applyTranslation(itg->t.controlPoints, itg->t.time);


		if (itg->r.angle)
			glRotatef(itg->r.angle, itg->r.axisX, itg->r.axisY, itg->r.axisZ);

		int nbuffers = 0;

		for (list<Model>::iterator itm = itg->m.begin(); itm != itg->m.end(); itm++)
		{
			float amb[4];
			float dif[4];
			float emi[4];

			amb[0] = itm->lc[nbuffers].ambR;
			amb[1] = itm->lc[nbuffers].ambG;
			amb[2] = itm->lc[nbuffers].ambB;
			amb[3] = 0.0f;
			glMaterialfv(GL_FRONT, GL_AMBIENT, amb);

			dif[0] = itm->lc[nbuffers].difR;
			dif[1] = itm->lc[nbuffers].difG;
			dif[2] = itm->lc[nbuffers].difB;
			dif[3] = 0.0f;
			glMaterialfv(GL_FRONT, GL_DIFFUSE, dif);

			emi[0] = itm->lc[nbuffers].emiR;
			emi[1] = itm->lc[nbuffers].emiG;
			emi[2] = itm->lc[nbuffers].emiB;
			emi[3] = 1.0f;
			glMaterialfv(GL_FRONT, GL_EMISSION, emi);

			glBindTexture(GL_TEXTURE_2D, itm->texID);

			drawTriangles(itm->buffers[nbuffers].b[0], itm->buffers[nbuffers].b[1], itm->buffers[nbuffers].b[2], itm->numberPoints[nbuffers]);

			glBindTexture(GL_TEXTURE_2D, 0);

			nbuffers++;
		}
		drawGroupElements(itg->g);
		glPopMatrix();
	}
}

void enableLighting()
{
	string type3("SPOTLIGHT");

	if (numLights > 0)
		glEnable(GL_LIGHTING);

	for (int i = 0; i < numLights; i++)
	{
		glLightfv(GL_LIGHT0 + i, GL_POSITION, lights[i].pos);
		glLightfv(GL_LIGHT0 + i, GL_AMBIENT, lights[i].amb);
		glLightfv(GL_LIGHT0 + i, GL_DIFFUSE, lights[i].dif);
		if (type3.compare(lights[i].type) == 0)
		{
			glLightfv(GL_LIGHT0, GL_SPOT_DIRECTION, lights[i].spotDir);
			glLightfv(GL_LIGHT0, GL_SPOT_CUTOFF, &lights[i].angle);
			glLightfv(GL_LIGHT0, GL_SPOT_EXPONENT, &lights[i].exp);
		}
		glEnable(GL_LIGHT0 + i);
	}
}


/*
Função que desenha as primitivas gráficas.
*/
void renderScene(void)
{

	glClearColor(0.0f, 0.0f, 0.0f, 0.0f);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	glEnable(GL_CULL_FACE);

	glLoadIdentity();
	gluLookAt(camX, camY, camZ,
		0.0, 0.0, 0.0,
		0.0f, 1.0f, 0.0f);

	glEnable(GL_TEXTURE_2D);

	enableLighting();

	glPolygonMode(GL_FRONT, GL_FILL);

	drawGroupElements(groups);

	time = glutGet(GLUT_ELAPSED_TIME);
	frame++;
	time = glutGet(GLUT_ELAPSED_TIME);
	if (time - timebase > 1000)
	{
		fps = frame*1000.0 / (time - timebase);
		timebase = time;
		frame = 0;
	}
	string s = "cg@di   FPS:";
	string f = to_string(fps);
	s += f;
	glutSetWindowTitle(s.c_str());

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

	startX, startY, tracking = 0;

	alpha = 0; beta = 0; r = 20;
	calculateCam();

	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGBA);
	glutInitWindowPosition(0, 0);
	glutInitWindowSize(1500, 700);
	glutCreateWindow("cg@di");


	glewInit();


	glEnable(GL_RESCALE_NORMAL);
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