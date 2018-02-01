#define _CRT_SECURE_NO_WARNINGS /*FASE 3->Motor*/
#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include "Vertex.h"
#define _USE_MATH_DEFINES
#include <math.h>
#include <ctype.h>
#include <stdlib.h>

using namespace std;

/*
	Função que escreve os vértices de uma primitiva gráfica num ficheiro cujo nome está especificado no argumento.
*/
void toFile(vector<vertex> vertices, string name) {
	FILE *file;
	string filename("..\\..\\Solids\\");
	filename += name;
	file = fopen(filename.c_str(), "w");
	fprintf(file, "%d\n", vertices.size());
	for (int i = 0; i <(int)vertices.size(); i++)
	{
		fprintf(file, "%f %f %f\n", vertices[i].x, vertices[i].y, vertices[i].z);
	}
	fclose(file);

	cout << "Generated with success!\n";
}

/*
	Função que calcula os pontos de dois triângulos na face do lado direito da box.
*/
void createTrianglesXRight(vector<vertex> *v, float x, float yf, float zf, float ry, float rz) {
	vertex p;
	p.x = x;

	// right hand rule
	// first triangle
	p.y = yf;
	p.z = zf;
	v->push_back(p);

	p.y = yf;
	p.z = zf - rz;
	v->push_back(p);

	p.y = yf + ry;
	p.z = zf - rz;
	v->push_back(p);

	// second triangle
	p.y = yf;
	p.z = zf;
	v->push_back(p);

	p.y = yf + ry;
	p.z = zf - rz;
	v->push_back(p);

	p.y = yf + ry;
	p.z = zf;
	v->push_back(p);
}

/*
	Função que calcula os pontos de dois triângulos na face do lado esquerdo da box.
*/
void createTrianglesXLeft(vector<vertex> *v, float x, float yf, float zf, float ry, float rz) {
	vertex p;
	p.x = x;

	// right hand rule
	// first triangle
	p.y = yf;
	p.z = zf;
	v->push_back(p);

	p.y = yf + ry;
	p.z = zf - rz;
	v->push_back(p);

	p.y = yf;
	p.z = zf - rz;
	v->push_back(p);

	// second triangle
	p.y = yf;
	p.z = zf;
	v->push_back(p);

	p.y = yf + ry;
	p.z = zf;
	v->push_back(p);

	p.y = yf + ry;
	p.z = zf - rz;
	v->push_back(p);

}

/*
	Função que calcula os pontos de dois triângulos na face superior da box.
*/
void createTrianglesYUp(vector<vertex> *v, float xf, float y, float zf, float rx, float rz) {
	vertex p;
	p.y = y;

	// right hand rule
	// first triangle
	p.x = xf;
	p.z = zf;
	v->push_back(p);

	p.x = xf + rx;
	p.z = zf;
	v->push_back(p);

	p.x = xf + rx;
	p.z = zf - rz;
	v->push_back(p);

	// second triangle
	p.x = xf;
	p.z = zf;
	v->push_back(p);

	p.x = xf + rx;
	p.z = zf - rz;
	v->push_back(p);

	p.x = xf;
	p.z = zf - rz;
	v->push_back(p);
}

/*
	Função que calcula os pontos de dois triângulos na face inferior da box.
*/
void createTrianglesYDown(vector<vertex> *v, float xf, float y, float zf, float rx, float rz) {
	vertex p;
	p.y = y;

	// right hand rule
	// first triangle
	p.x = xf;
	p.z = zf;
	v->push_back(p);

	p.x = xf + rx;
	p.z = zf - rz;
	v->push_back(p);

	p.x = xf + rx;
	p.z = zf;
	v->push_back(p);

	// second triangle
	p.x = xf;
	p.z = zf;
	v->push_back(p);

	p.x = xf;
	p.z = zf - rz;
	v->push_back(p);

	p.x = xf + rx;
	p.z = zf - rz;
	v->push_back(p);
}

/*
	Função que calcula os pontos de dois triângulos na face da frente da box.
*/
void createTrianglesZFront(vector<vertex> *v, float xf, float yf, float z, float rx, float ry) {
	vertex p;
	p.z = z;

	// right hand rule
	// first triangle
	p.x = xf;
	p.y = yf;
	v->push_back(p);

	p.x = xf + rx;
	p.y = yf;
	v->push_back(p);

	p.x = xf + rx;
	p.y = yf + ry;
	v->push_back(p);

	// second triangle
	p.x = xf + rx;
	p.y = yf + ry;
	v->push_back(p);

	p.x = xf;
	p.y = yf + ry;
	v->push_back(p);

	p.x = xf;
	p.y = yf;
	v->push_back(p);
}

/*
	Função que calcula os pontos de dois triângulos na face de trás da box.
*/
void createTrianglesZBack(vector<vertex> *v, float xf, float yf, float z, float rx, float ry) {
	vertex p;
	p.z = z;

	// right hand rule
	// first triangle
	p.x = xf;
	p.y = yf;
	v->push_back(p);

	p.x = xf + rx;
	p.y = yf + ry;
	v->push_back(p);

	p.x = xf + rx;
	p.y = yf;
	v->push_back(p);

	// second triangle
	p.x = xf + rx;
	p.y = yf + ry;
	v->push_back(p);

	p.x = xf;
	p.y = yf;
	v->push_back(p);

	p.x = xf;
	p.y = yf + ry;
	v->push_back(p);




}

/*
	Função que calcula todos os pontos da face superior da box.
*/
void faceYUp(vector<vertex> *v, float x, float y, float z, float dx, float dz, int d) {
	float rx = dx / d;
	float rz = dz / d;
	float xf;
	float zf;

	int i, j;
	xf = x;
	for (i = 0; i < d; i++)
	{
		zf = z;
		for (j = 0; j < d; j++)
		{
			createTrianglesYUp(v, xf, y, zf, rx, rz);
			zf -= rz;
		}
		xf += rx;
	}
}

/*
	Função que calcula todos os pontos da face inferior da box.
*/
void faceYDown(vector<vertex> *v, float x, float y, float z, float dx, float dz, int d) {
	float rx = dx / d;
	float rz = dz / d;
	float xf;
	float zf;

	int i, j;
	xf = x;
	for (i = 0; i < d; i++)
	{
		zf = z;
		for (j = 0; j < d; j++)
		{
			createTrianglesYDown(v, xf, y, zf, rx, rz);
			zf -= rz;
		}
		xf += rx;
	}
}

/*
	Função que calcula todos os pontos da face de cima da box.
*/
void faceZFront(vector<vertex> *v, float x, float y, float z, float dx, float dy, int d) {
	float rx = dx / d;
	float ry = dy / d;
	float xf;
	float yf;
	
	yf = y;
	
	int i, j;

	for (i = 0; i < d; i++)
	{
		xf = x;
		for (j = 0; j < d; j++)
		{
			createTrianglesZFront(v, xf, yf, z, rx, ry);

			xf += rx;
		}
		yf += ry;
	}
}

/*
	Função que calcula todos os pontos da face de trás da box.
*/
void faceZBack(vector<vertex> *v, float x, float y, float z, float dx, float dy, int d) {
	float rx = dx / d;
	float ry = dy / d;
	float xf;
	float yf;

	yf = y;

	int i, j;

	for (i = 0; i < d; i++)
	{
		xf = x;
		for (j = 0; j < d; j++)
		{
			createTrianglesZBack(v, xf, yf, z, rx, ry);

			xf += rx;
		}
		yf += ry;
	}
}

/*
	Função que calcula todos os pontos da face do lado direito da box.
*/
void faceXRight(vector<vertex> *v, float x, float y, float z, float dy, float dz, int d) {
	float ry = dy / d;
	float rz = dz / d;
	float yf;
	float zf;
	int i, j;


	yf = y;
	for (i = 0; i < d; i++)
	{
		zf = z;
		for (j = 0; j < d; j++)
		{
			createTrianglesXRight(v, x, yf, zf, ry, rz);
			zf -= rz;
		}
		yf += ry;
	}

}

/*
	Função que calcula todos os pontos da face do lado esquerdo da box.
*/
void faceXLeft(vector<vertex> *v, float x, float y, float z, float dy, float dz, int d) {
	float ry = dy / d;
	float rz = dz / d;
	float yf;
	float zf;
	int i, j;


	yf = y;
	for (i = 0; i < d; i++)
	{
		zf = z;
		for (j = 0; j < d; j++)
		{
			createTrianglesXLeft(v, x, yf, zf, ry, rz);
			zf -= rz;
		}
		yf += ry;
	}
}

/*
	Função que calcula todos os pontos do plano.
*/
vector<vertex> plane(float coord) {
	vector<vertex> vertices;
	float coordF = coord / 2;
	vertex v;
	
	v.x = -coordF;
	v.y = 0;
	v.z = coordF;
	vertices.push_back(v);

	v.x = coordF;
	v.y = 0;
	v.z = coordF;
	vertices.push_back(v);


	v.x = coordF;
	v.y = 0;
	v.z = -coordF;
	vertices.push_back(v);
	
	v.x = -coordF;
	v.y = 0;
	v.z = coordF;
	vertices.push_back(v);

	v.x = coordF;
	v.y = 0;
	v.z = -coordF;
	vertices.push_back(v);

	v.x = -coordF;
	v.y = 0;
	v.z = -coordF;
	vertices.push_back(v);
	
	return vertices;
}

/*
	Função que calcula todos os pontos da box, utilizando funções auxiliares para cada face. Os pontos foram inicialmente gerados
no octante positivo e depois centrados na origem.
*/
vector<vertex> box(float x, float y, float z, int d) {
	vector<vertex> vertices;

	faceXRight(&vertices, x, 0, z, y, z, d);
	faceXLeft(&vertices, 0, 0, z, y, z, d);

	faceYUp(&vertices, 0, y, z, x, z, d);
	faceYDown(&vertices, 0, 0, z, x, z, d);


	faceZFront(&vertices, 0, 0, z, x, y, d);
	faceZBack(&vertices, 0, 0, 0, x, y, d);
	
	for (int i = 0; i <(int)vertices.size(); i++)
	{
		vertices[i].x -= x / 2;
		vertices[i].y -= y / 2;
		vertices[i].z -= z / 2;
	}
	
	return vertices;
}

/*
	Função que calcula todos os pontos da esfera.
*/
vector<vertex> sphere(float radius, int slices, int stacks) {

	vector<vertex> v;

	float alpha = 2 * (float)M_PI / slices;
	float rBeta = (float)M_PI / stacks;

	int i, j;

	float beta;

	vertex a, b, c, d;

	beta = (float)M_PI / 2;

	for (i = 0; i < stacks; i++)
	{
		
		for (j = 0; j < slices; j++)
		{
			a.x = radius * cos(beta) * sin(j * alpha);
			a.y = radius * sin(beta);
			a.z = radius * cos(beta) * cos(j * alpha);

			b.x = radius * cos(beta) * sin(j * alpha + alpha);
			b.y = radius * sin(beta);
			b.z = radius * cos(beta) * cos(j * alpha + alpha);

			c.x = radius * cos(beta - rBeta) * sin(j * alpha + alpha);
			c.y = radius * sin(beta - rBeta);
			c.z = radius * cos(beta - rBeta) * cos(j * alpha + alpha);

			d.x = radius * cos(beta - rBeta) * sin(j * alpha);
			d.y = radius * sin(beta - rBeta);
			d.z = radius * cos(beta - rBeta) * cos(j * alpha);

			v.push_back(a);
			v.push_back(b);
			v.push_back(c);

			v.push_back(a);
			v.push_back(c);
			v.push_back(d);

		}
		beta -= rBeta;
	}



	return v;
}

/*
	Função auxiliar que calcula todos os pontos da base do cone.
*/
vector<vertex> base(float radius, int slices)
{
	vector<vertex> v;
	vertex vertex, vOrigin;
	vertex.y = 0;
	int i;
	float angle;

	if (slices == 1)
	{
		angle = (2 * (float) M_PI)/3;
		for (i = 3; i > 0; i--)
		{
			vertex.x = radius * sin(i* angle);
			vertex.z = radius * cos(i* angle);
			v.push_back(vertex);
		}
	}
	else if (slices == 2)
	{
		angle = (float) M_PI / 2;
		for (i = 2; i >= 0; i-=2 )
		{
			vertex.x = radius*sin(i*angle);
			vertex.z = radius * cos(i* angle);
			v.push_back(vertex);

			vertex.x = radius*sin(i*angle - angle);
			vertex.z = radius * cos(i* angle - angle);
			v.push_back(vertex);

			vertex.x = radius*sin(i*angle - 2 * angle);
			vertex.z = radius * cos(i* angle - 2 * angle);
			v.push_back(vertex);
		}
	}
	else
	{
		angle = 2 * (float) M_PI / slices;

		vOrigin.x = 0;
		vOrigin.y = 0;
		vOrigin.z = 0;

		for (i = slices; i > 0; i--)
		{
			vertex.x = radius * sin(i * angle + angle);
			vertex.z = radius * cos(i * angle + angle);

			v.push_back(vertex);

			vertex.x = radius * sin(i*angle);
			vertex.z = radius * cos(i*angle);

			v.push_back(vertex);

			v.push_back(vOrigin);

		}
	}
	
	return v;
}

/*
	Função auxiliar que calcula todos os pontos laterais do cone.
*/
vector<vertex> lateral(vector<vertex> v, float radius, float height, int slices, int stacks)
{
	float rh = height / stacks;
	float angle = 2 * (float)M_PI / slices;
	float old_r = radius;
	int i, j, n_slices = slices;
	float new_r;

	vertex a, b, c, d;

	if (slices == 1)
	{

		angle = (2 * (float)M_PI) / 3;
		n_slices = 3;
	}
	if (slices == 2)
	{
		angle = (float)M_PI / 2;
		n_slices = 4;
	}
	for (i = 0; i < stacks - 1; i++)
	{
		for (j = 0; j < n_slices; j++)
		{
			a.x = old_r *  sin(j*angle);
			a.y = i*rh;
			a.z = old_r *  cos(j*angle);


			b.x = old_r * sin(j * angle + angle);
			b.y = i * rh;
			b.z = old_r * cos(j * angle + angle);

			new_r = (radius * (height - (i + 1) * rh)) / height;

			c.x = new_r * sin(j * angle + angle);
			c.y = (i + 1) * rh;
			c.z = new_r * cos(j * angle + angle);

			d.x = new_r *  sin(j*angle);
			d.y = (i + 1)*rh;
			d.z = new_r *  cos(j*angle);


			v.push_back(a);
			v.push_back(b);
			v.push_back(c);

			v.push_back(a);
			v.push_back(c);
			v.push_back(d);
		}
		old_r = new_r;

	}
	c.x = 0;
	c.y = height;
	c.z = 0;

	for (int i = 0; i < n_slices; i++)
	{
			
		a.x = old_r *  sin(i*angle);
		a.y = height - rh;
		a.z = old_r *  cos(i*angle);


		b.x = old_r * sin(i * angle + angle);
		b.y = height - rh;
		b.z = old_r * cos(i * angle + angle);

		v.push_back(a);
		v.push_back(b);
		v.push_back(c);
		
	}
	return v;
}

/*
	Função auxiliar que calcula todos os pontos do cone.
*/
vector<vertex> cone(float radius, float height, int slices, int stacks)
{
	vector<vertex> v= base(radius, slices);


	return lateral(v, radius, height, slices, stacks);

}

/*
	Função que calcula todos os pontos do anel 
*/

vector<vertex> ring(vector<vertex> vertices, float radiusB, float radiusS, int slices)
{
	vertex v;
	v.y = 0;

	

	float angle = 2 * (float)M_PI / slices;

	float ang1 = 0;
	float ang2 = ang1 - angle;
	for(int i = slices; i > 0; i--)
	{		
		ang2 = ang1 - angle;
		v.x = radiusB*sin(ang1);
		v.z = radiusB*cos(ang1);

		vertices.push_back(v);

		v.x = radiusB*sin(ang2);
		v.z = radiusB*cos(ang2);

		vertices.push_back(v);

		v.x = radiusS*sin(ang1);
		v.z = radiusS*cos(ang1);

		vertices.push_back(v);

		v.x = radiusB*sin(ang2);
		v.z = radiusB*cos(ang2);

		vertices.push_back(v);

		v.x = radiusS*sin(ang2);
		v.z = radiusS*cos(ang2);

		vertices.push_back(v);

		v.x = radiusS*sin(ang1);
		v.z = radiusS*cos(ang1);

		vertices.push_back(v);

		v.x = radiusB*sin(ang2);
		v.z = radiusB*cos(ang2);

		vertices.push_back(v);

		v.x = radiusB*sin(ang1);
		v.z = radiusB*cos(ang1);

		vertices.push_back(v);

		v.x = radiusS*sin(ang1);
		v.z = radiusS*cos(ang1);

		vertices.push_back(v);
		

		v.x = radiusB*sin(ang2);
		v.z = radiusB*cos(ang2);

		vertices.push_back(v);

		v.x = radiusS*sin(ang1);
		v.z = radiusS*cos(ang1);

		vertices.push_back(v);

		v.x = radiusS*sin(ang2);
		v.z = radiusS*cos(ang2);

		vertices.push_back(v);

		ang1 -= angle;
	}

	return vertices;
}
/*
	Função que se encarrega das verificações do input. Se algum parâmetro for incorreto, retorna falso.
*/
bool checkInput(int argc, char** argv) {
	bool ok = true;
	try {
		for (int i = 2; i < argc - 1; i++) {
			stof(argv[i]);
		}
	}
	catch (const invalid_argument&) {
		ok = false;
		cout << "The inserted values aren't valid.\n";
	}

	string model = argv[1];
	if (ok)
	{
		if (model.compare("plane") == 0) {
			if (argc != 4)
			{
				ok = false;
				cout << "Invalid number of arguments.\n";
			}
			else
			{
				try
				{
					if (stof(argv[2]) <= 0)
					{
						cout << "Invalid dimension.\n";
						ok = false;
					}
				}
				catch (const invalid_argument&)
				{
					ok = false;
					cout << "The inserted values aren't valid.\n";
				}
			}
		}
		else if (model.compare("box") == 0) {
			if (argc != 7)
			{
				ok = false;
				cout << "Invalid number of arguments.\n";
			}
			else
			{
				try
				{
					if (stof(argv[2]) <= 0 || stof(argv[3]) <= 0 || stof(argv[4]) <= 0 || stoi(argv[5]) <=0)
					{
						cout << "Invalid dimension.\n";
						ok = false;
					}
				}
				catch (const invalid_argument&)
				{
					ok = false;
					cout << "The inserted values aren't valid.\n";
				}
			}
		}
		else if (model.compare("sphere") == 0) {
			if (argc != 6)
			{
				ok = false;
				cout << "Invalid number of arguments.\n";
			}
			else
			{
				try
				{
					if (stof(argv[2]) <= 0 || stoi(argv[3]) < 2 || stoi(argv[4]) < 2)
					{
						cout << "Invalid dimension.\n";
						ok = false;
					}
				}
				catch (const invalid_argument&)
				{
					ok = false;
					cout << "The inserted values aren't valid.\n";
				}
			}
		}
		else if (model.compare("cone") == 0) {
			if (argc != 7)
			{
				ok = false;
				cout << "Invalid number of arguments.\n";
			}
			else
			{
				try
				{
					if (stof(argv[2]) <= 0 || stof(argv[3]) <= 0 || stoi(argv[4]) <= 0 || stoi(argv[5]) <= 0)
					{
						cout << "Invalid dimension.\n";
						ok = false;
					}
				}
				catch (const invalid_argument&)
				{
					ok = false;
					cout << "The inserted values aren't valid.\n";
				}
			}
		}
		else if (model.compare("ring") == 0) {
			if (argc != 6)
			{
				ok = false;
				cout << "Invalid number of arguments.\n";
			}
			else
			{
				try
				{
					if (stof(argv[2]) <= 0 || stof(argv[3]) <= 0 || stoi(argv[4]) <= 2)
					{
						cout << "Invalid dimension.\n";
						ok = false;
					}
				}
				catch (const invalid_argument&)
				{
					ok = false;
					cout << "The inserted values aren't valid.\n";
				}
			}
		}
		else {
			cout << "Inexistent graphical primitive.\n";
		}
	}
	return ok;
}

int main(int argc, char **argv) {
	string model;
	string file;

	if (argc > 2)
	{
		model = argv[1];

		if (checkInput(argc, argv) == true) {
			if (model.compare("plane") == 0) {
				float coord = stof(argv[2]);
				toFile(plane(coord), argv[argc - 1]);
			}
			else if (model.compare("box") == 0) {
				float x = stof(argv[2]);
				float y = stof(argv[3]);
				float z = stof(argv[4]);
				int div = stoi(argv[5]);
				vector<vertex> v = box(x, y, z, div);
				toFile(v, argv[argc - 1]);
			}
			else if (model.compare("sphere") == 0) {
				float radius = stof(argv[2]);
				int slices = stoi(argv[3]);
				int stacks = stoi(argv[4]);
				vector<vertex> v = sphere(radius, slices, stacks);
				toFile(v, argv[argc - 1]);
			}
			else if (model.compare("cone") == 0) {
				float radius = stof(argv[2]);
				float height = stof(argv[3]);
				int slices = stoi(argv[4]);
				int stacks = stoi(argv[5]);

				vector<vertex> v = cone(radius, height, slices, stacks);
				toFile(v, argv[argc - 1]);
			}
			else if (model.compare("ring") == 0)
			{
				float radiusB = stof(argv[2]);
				float radiusS = stof(argv[3]);
				int slices = stoi(argv[4]);

				float innerRing = (radiusB + radiusS) / 2;
				vector<vertex> v;
				v = ring(v, innerRing, radiusS, slices);
				v = ring(v,radiusB, innerRing, slices);
				toFile(v, argv[argc - 1]);
			}
		}
	}
	else {
		cout << "Not enough arguments.\n";
	}
	return 0;
}