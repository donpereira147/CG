#pragma once
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
	float angle;
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

struct Model
{
	string name;
	int numberPoints[20];
	GLuint vertices[20];
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