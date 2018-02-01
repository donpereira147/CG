// Generator.cpp : Defines the entry point for the console application.
//
#define _CRT_SECURE_NO_WARNINGS

#include "./tinyxml/tinyxml.h"
#include "./tinyxml/tinystr.h"
#include "Vertex.h"
#include "Group.h"
#include <vector>
#include <fstream>
#include <iostream>
#include <string>
#include <GL/glut.h> 
#include <math.h>
#define _USE_MATH_DEFINES
#include <math.h>
#include <ctype.h>
#include <map>
#include "Movement.h"

#define ASTEROIDS 100

using namespace std;

list<Group> groups;

float earth_degrees, moon_degrees,sun_degrees, mercury_degrees, venus_degrees, mars_degrees, jupiter_degrees,
saturn_degrees, uranus_degrees, neptune_degrees, pluto_degrees;


vector<Movement> orbits;

float height = 2.0f;
float x = 0.0f;
float z = 0.0f;


float camX = 00, camY = 30, camZ = 400;
int startX, startY, tracking = 0;

int alpha = 0, beta = 0, r = 400;

float rBeltI, rBeltS;

float earthRadius;

/*
	Inicializa a estrutura Group
*/

Group initGroup()
{
	Translate t;
	t.x = t.y = t.z = 0;
	
	Rotate r ;
	r.angle = 0;
	r.axisX = 0;
	r.axisY = 0;
	r.axisZ = 0;

	
	Scale s;
	s.x = s.y = s.z = 1;

	Colour c;
	c.r = c.b = c.g = 1;

	Group g;
	g.t = t;
	g.r = r;
	g.s = s;
	g.c = c;
	g.orbitRadius = 0.0;
	return g;
}

/*
	Faz a leitura do ficheiro XML e armazena os dados na lista de Groups
*/

void loadGroups(TiXmlNode *currentranslateNode, list<Group>* groups)
{
	for (TiXmlNode* groupNode = currentranslateNode->FirstChild("group"); groupNode; groupNode = groupNode->NextSiblingElement())
	{
		Group g = initGroup();

		TiXmlNode* earthRadiuscaleNode = groupNode->FirstChild("earth_radius");
		if (earthRadiuscaleNode != nullptr)
		{
			TiXmlElement* erElem = earthRadiuscaleNode->ToElement();
			const char *earth_radius = erElem->Attribute("R");
			if (earth_radius)
				earthRadius = atof(earth_radius);
		}

		TiXmlNode* translateNode = groupNode->FirstChild("translate");
		if (translateNode != nullptr)
		{
			TiXmlElement* tElem = translateNode->ToElement();
			const char *tx, *ty, *tz;
			tx = tElem->Attribute("X");
			ty = tElem->Attribute("Y");
			tz = tElem->Attribute("Z");

			if (tx)
				g.t.x = atof(tx);
			if (ty)
				g.t.y = atof(ty);
			if (tz)
				g.t.z = atof(tz);
		}

		TiXmlNode* rotateNode = groupNode->FirstChild("rotate");
		if (rotateNode != nullptr)
		{
			TiXmlElement* rElem = rotateNode->ToElement();
			const char *angle, *ax, *ay, *az;
			angle = rElem->Attribute("angle");
			ax = rElem->Attribute("axisX");
			ay = rElem->Attribute("axisY");
			az = rElem->Attribute("axisZ");

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

					if (f.is_open())
					{
						vector <vertex> solid;
						string line;
						vertex v;
						int numVertices;
						getline(f, line);
						numVertices = stoi(line.c_str());

						while (getline(f, line))
						{
							float x, y, z;
							sscanf(line.c_str(), "%f %f %f\n", &x, &y, &z);

							v.x = x;
							v.y = y;
							v.z = z;

							solid.push_back(v);

						}
						f.close();
						Model m;
						m.vertices = solid;
						m.name = name;
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

/*
	Chama a função que é responsável pela leitura do ficheiro XML
*/

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

/*
	Desenha as primitivas gráficas
*/

void drawTriangles(list<Model>::iterator itm)
{
	glBegin(GL_TRIANGLES);
	for (int i = 0; i < (int) itm->vertices.size(); i++)
	{
		glVertex3f(itm->vertices[i].x, itm->vertices[i].y, itm->vertices[i].z);
	}
	glEnd();
}

/*
	Desenha o planeta Terra
*/

void drawPlanetEarth(list<Model>::iterator itm)
{
	srand(30);
	int random;
	int g, b;
	glBegin(GL_TRIANGLES);
	for (int i = 0; i < (int) itm->vertices.size(); i++)
	{
		g = b = 0;
		random = rand() % 3;
		if ((random % 2)== 0)
			b = 1;
		else g = 1;
		glColor3f(0, g, b);
		glVertex3f(itm->vertices[i].x, itm->vertices[i].y, itm->vertices[i].z);
	}
	glEnd();
	
}

/*
	Desenha as órbitas dos planetas
*/

void drawOrbit(float r)
{
	float angleSlice = 2 * M_PI / 200;
	float angle = 0;

	glBegin(GL_POINTS);
	for (int i = 0; i < 200; i++)
	{
		glVertex3f(r*cos(angle), 0, r*sin(angle));
		angle += angleSlice;
	}
	glEnd();

}

/*
	Verifica se o asteroide está contido na cintura
*/

bool asteroidIn(float x, float z)
{
	float hip = sqrt(pow(x, 2) + pow(z, 2));
	return (hip > rBeltI && hip < rBeltS);
}

/*
	Desenha a cintura de asteroides
*/

void drawAsteroidBelt(list<Model>::iterator itm)
{
	srand(1000);
	int nasteroids = 0;
	float alpha, r, x, z;

	

	for (nasteroids = 0; nasteroids < ASTEROIDS; nasteroids++)
	{
		r = ((float)rand()/RAND_MAX) * (rBeltS - rBeltI);
		alpha = ((float)rand()/RAND_MAX) *(2*M_PI);

		x = sin(alpha) * (r+rBeltI);
		z = cos(alpha) * (r+ rBeltI);

		if (asteroidIn(x, z))
		{
			glPushMatrix();
			glTranslatef(x, 0, z);
			drawTriangles(itm);
			glPopMatrix();
		}
	}
}

/*
	Desenha os elementos do Sistema Solar contidos na lista de grupos
*/

void drawGroupElements(list<Group> g)
{
	glPointSize(1.1);
	string sun("sun.3d");
	string mercury("mercury.3d");
	string venus("venus.3d");
	string earth("earth.3d");
	string moon("moon.3d");
	string mars("mars.3d");
	string jupiter("jupiter.3d");
	string saturn("saturn.3d");
	string uranus("uranus.3d");
	string neptune("neptune.3d");
	string pluto("pluto.3d");
	string asteroids("asteroid.3d");

	bool planetEarth, asteroid;

	for (list<Group>::iterator itg = g.begin(); itg != g.end(); itg++)
	{
		glPushMatrix();
		glTranslatef(itg->t.x, itg->t.y, itg->t.z);
		glRotatef(itg->r.angle, itg->r.axisX, itg->r.axisY, itg->r.axisZ);
		glScalef(itg->s.x, itg->s.y, itg->s.z);
		glColor3f(itg->c.r, itg->c.g, itg->c.b);

		if(itg->orbitRadius != 0)
			drawOrbit(itg->orbitRadius);

		for (list<Model>::iterator itm = itg->m.begin(); itm != itg->m.end(); itm++)
		{
			planetEarth = false;
			asteroid = false;

			if (sun.compare(itm->name) == 0)
			{
				glRotatef(sun_degrees, 0, 1, 0);
			}
			else if (mercury.compare(itm->name) == 0)
			{
				glRotatef(mercury_degrees, 0, 1, 0);
			}
			else if (venus.compare(itm->name) == 0)
			{
				glRotatef(venus_degrees, 0, 1, 0);
			}
			else if (earth.compare(itm->name) == 0)
			{
				planetEarth = true;
				glRotatef(earth_degrees, 0, 1, 0);
			}
			else if (moon.compare(itm->name) == 0)
			{
				glRotatef(moon_degrees, 0, 1, 0);
			}
			else if (mars.compare(itm->name) == 0)
			{
				glRotatef(mars_degrees, 0, 1, 0);
				rBeltI = (itg->t.x + 5) + earthRadius*itg->s.x;
			}
			else if (jupiter.compare(itm->name) == 0)
			{
				glRotatef(jupiter_degrees, 0, 1, 0);
				rBeltS = (itg->t.x-5) - earthRadius*itg->s.x;
			}
			else if (saturn.compare(itm->name) == 0)
			{
				glRotatef(saturn_degrees, 0, 1, 0);
			}
			else if (uranus.compare(itm->name) == 0)
			{
				glRotatef(uranus_degrees, 0, 1, 0);
			}
			else if (neptune.compare(itm->name) == 0)
			{
				glRotatef(neptune_degrees, 0, 1, 0);
			}
			else if (pluto.compare(itm->name) == 0)
			{
				glRotatef(pluto_degrees, 0, 1, 0);
			}
			else if (asteroids.compare(itm->name) == 0)
			{
				asteroid = true;
			}

			if (planetEarth)
				drawPlanetEarth(itm);
			else if(asteroid)
				drawAsteroidBelt(itm);
			else
				drawTriangles(itm);
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

/*
	Movimento de rotação dos planetas
*/

void rotate(unsigned char key, int x, int y)
{
	switch (key)
	{
		case 'a':
			sun_degrees += 10;
			mercury_degrees += 10;
			earth_degrees += 5;
			venus_degrees += 7;
			mars_degrees += 4;
			jupiter_degrees += 2.5;
			saturn_degrees += 1.7;
			uranus_degrees += 1.3;
			neptune_degrees += 1.2;
			pluto_degrees += 1.1;
			moon_degrees += 12;
			break;
		case 'q':
			mercury_degrees += 10;
			break;
		case 'w':
			venus_degrees += 7;
			break;
		case 'e':
			earth_degrees += 5;
			break;
		case 'r':
			mars_degrees += 4;
			break;
		case 't':
			jupiter_degrees += 2.5;
			break;
		case 'y':
			saturn_degrees += 1.7;
			break;
		case 'u':
			uranus_degrees += 1.3;
			break;
		case 'i':
			neptune_degrees += 1.2;
			break;
		case 'o':
			pluto_degrees += 1.1;
			break;
		case 'p':
			moon_degrees += 12;
			break;
		default: break;
	}
	glutPostRedisplay();
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

	rBeltI = rBeltS = 0.0;
	sun_degrees = mercury_degrees = venus_degrees = earth_degrees =
			mars_degrees = jupiter_degrees = saturn_degrees = uranus_degrees = pluto_degrees = neptune_degrees =  0;
	string filename("..\\..\\Solids\\");
	filename += argv[1];
	groups = load(filename.c_str());

	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGBA);
	glutInitWindowPosition(0, 0);
	glutInitWindowSize(1500, 700);
	glutCreateWindow("cg@di");

	glutDisplayFunc(renderScene);
	glutIdleFunc(renderScene);
	glutReshapeFunc(changeSize);

	glutKeyboardFunc(rotate);
	glutMouseFunc(processMouseButtons);
	glutMotionFunc(processMouseMotion);

	glEnable(GL_DEPTH_TEST);
	glEnable(GL_CULL_FACE);
	glClearColor(0.0f, 0.0f, 0.0f, 0.0f);

	glutMainLoop();

	return 1;
}

