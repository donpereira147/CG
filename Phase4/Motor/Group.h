#pragma once

#include "Vertex.h"
#include <vector>
#include <list>
#include <GL/glew.h>


using namespace std;

struct Translate
{
	float time;
	float x, y, z;
	vector<vertex> controlPoints;
};

struct Rotate
{
	float time;
	float angle, rotationAngle;
	int axisX, axisY, axisZ;
};

struct Colour
{
	float r, g, b;
};

struct Scale
{
	float x, y, z;
};

struct Buffers {
	GLuint b[3];
};


struct LightComp{
	float ambR, ambG, ambB;
	float difR, difG, difB;
	float emiR, emiG, emiB;
};

struct Model
{
	string texture;
	int texID;
	int numberPoints[20];
	Buffers buffers[20];
	LightComp lc[20];
};


struct Group
{
	Rotate r;
	Translate t;
	Scale s;
	Colour c;
	float orbitRadius;
	bool orientation;
	bool points;
	list<Model> m;
	list<Group> g;
};