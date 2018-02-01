#pragma once

#include <string>
#include <vector>

using namespace std;

struct Light {
	string type;
	float pos[4];
	float amb[4];
	float dif[4];
	float spotDir[3];
	GLfloat angle;
	GLfloat exp;
};
