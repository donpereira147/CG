#pragma once
#pragma once

#include "Vertex.h"
#include <vector>
#include <list>


using namespace std;

struct Translate
{
	float x, y, z;
};

struct Rotate
{
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
	vector<vertex> vertices;
};



struct Group
{
	Rotate r;
	Translate t;
	Scale s;
	Colour c;
	float orbitRadius;
	list<Model> m;
	list<Group> g;
};