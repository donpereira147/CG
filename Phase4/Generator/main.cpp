#define _CRT_SECURE_NO_WARNINGS /*FASE 4 -> Generator*/
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


int numPatches, numPoints;

float px[4][4];
float py[4][4];
float pz[4][4];

void toFile(vector<vertex> vertices, string name) {
	FILE *file;
	string filename("..\\..\\Solids\\");
	filename += name;
	file = fopen(filename.c_str(), "w");
	fprintf(file, "%d\n", vertices.size());
	int j = 0;
	for (int i = 0; i <(int)vertices.size(); i++)
	{
		fprintf(file, "%f %f %f\n", vertices[i].x, vertices[i].y, vertices[i].z);
	}
	fclose(file);

	cout << "Generated with success!\n";
}

void toFile(vector<vertex> vertices, vector<vertex> normals,  string name) {
	FILE *file;
	string filename("..\\..\\Solids\\");
	filename += name;
	file = fopen(filename.c_str(), "w");
	fprintf(file, "%d\n", vertices.size());
	int j = 0;
	for (int i = 0; i <(int)vertices.size(); i++)
	{
		fprintf(file, "%f %f %f %f %f %f\n", vertices[i].x, vertices[i].y, vertices[i].z, normals[i].x, normals[i].y, normals[i].z);
	}
	fclose(file);

	cout << "Generated with success!\n";
}


/*
Função que escreve os vértices de uma primitiva gráfica num ficheiro cujo nome está especificado no argumento.
*/

void toFile(vector<vertex> vertices, vector<vertex> normals, vector<vertex> tex, string name) {
	FILE *file;
	string filename("..\\..\\Solids\\");
	filename += name;
	file = fopen(filename.c_str(), "w");
	fprintf(file, "%d\n", vertices.size());
	int j = 0;
	for (int i = 0; i <(int)vertices.size(); i++)
	{
		fprintf(file, "%f %f %f %f %f %f %f %f\n", vertices[i].x, vertices[i].y, vertices[i].z, normals[i].x, normals[i].y, normals[i].z,
													tex[i].x, tex[i].y);
	}
	fclose(file);
	
	cout << "Generated with success!\n";
}

/*
Função que calcula todos os pontos do plano.
*/
vector<vertex> plane(float coord, vector<vertex> *n, vector<vertex> *tex)
{
	vector<vertex> vertices;
	float coordF = coord / 2;
	vertex v;

	vertex normal;
	normal.x = 0;
	normal.y = 1;
	normal.z = 0;

	vertex texCoord;
	//a
	v.x = -coordF;
	v.y = 0;
	v.z = coordF;

	texCoord.y = 1;
	texCoord.x = 0;
	tex->push_back(texCoord);
	vertices.push_back(v);
	n->push_back(normal);
	//b
	v.x = coordF;
	v.y = 0;
	v.z = coordF;
	texCoord.x = 1;
	texCoord.y = 1;
	tex->push_back(texCoord);
	vertices.push_back(v);
	n->push_back(normal);

	//c
	v.x = coordF;
	v.y = 0;
	v.z = -coordF;
	texCoord.x = 1;
	texCoord.y = 0;
	tex->push_back(texCoord);
	vertices.push_back(v);
	n->push_back(normal);

	//a
	v.x = -coordF;
	v.y = 0;
	v.z = coordF;
	texCoord.x = 0;
	texCoord.y = 1;
	tex->push_back(texCoord);
	vertices.push_back(v);
	n->push_back(normal);

	//c
	v.x = coordF;
	v.y = 0;
	v.z = -coordF;
	texCoord.x = 1;
	texCoord.y = 0;
	tex->push_back(texCoord);
	vertices.push_back(v);
	n->push_back(normal);

	//d
	v.x = -coordF;
	v.y = 0;
	v.z = -coordF;
	texCoord.x = 0;
	texCoord.y = 0;
	tex->push_back(texCoord);
	vertices.push_back(v);
	n->push_back(normal);

	
	normal.x = 0;
	normal.y = -1;
	normal.z = 0;

	v.x = -coordF;
	v.y = 0;
	v.z = coordF;
	vertices.push_back(v);
	n->push_back(normal);
	texCoord.y = 1;
	texCoord.x = 0;
	tex->push_back(texCoord);

	v.x = coordF;
	v.y = 0;
	v.z = -coordF;
	vertices.push_back(v);
	n->push_back(normal);
	texCoord.x = 1;
	texCoord.y = 0;
	tex->push_back(texCoord);

	v.x = coordF;
	v.y = 0;
	v.z = coordF;
	vertices.push_back(v);
	n->push_back(normal);
	texCoord.x = 1;
	texCoord.y = 1;
	tex->push_back(texCoord);

	v.x = coordF;
	v.y = 0;
	v.z = -coordF;
	vertices.push_back(v);
	n->push_back(normal);
	texCoord.x = 1;
	texCoord.y = 0;
	tex->push_back(texCoord);

	v.x = -coordF;
	v.y = 0;
	v.z = coordF;
	vertices.push_back(v);
	n->push_back(normal);
	texCoord.y = 1;
	texCoord.x = 0;
	tex->push_back(texCoord);

	v.x = -coordF;
	v.y = 0;
	v.z = -coordF;
	vertices.push_back(v);
	n->push_back(normal);
	texCoord.x = 0;
	texCoord.y = 0;
	tex->push_back(texCoord);
	
	return vertices;
}

void normalize(float *a) {

	float l = sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2]);
	a[0] = a[0] / l;
	a[1] = a[1] / l;
	a[2] = a[2] / l;
}

void cross(float *a, float *b, float *res) {

	res[0] = a[1] * b[2] - a[2] * b[1];
	res[1] = a[2] * b[0] - a[0] * b[2];
	res[2] = a[0] * b[1] - a[1] * b[0];
}

void multMatrixVector(float *m, float *v, float *res) {

	for (int j = 0; j < 4; ++j) {
		res[j] = 0;
		for (int k = 0; k < 4; ++k) {
			res[j] += v[k] * m[j * 4 + k];
		}
	}

}

float multVectorVector(float *v1, float *v2)
{
	float res = 0;
	for (int i = 0; i < 4; ++i)
		res += v1[i] * v2[i];
	return res;
}

vector<vertex> createTeapot(vector<vertex> *n, vector<vertex> *tex, char* file, int divU, int divV)
{
	string line;
	int* indexes;
	float* points;
	int patch, point;
	string filename("..\\..\\Solids\\");
	filename += file;
	fstream f;
	f.open(filename);
	int i = 0, j = 0;


	if (f.is_open())
	{
		getline(f, line);
		numPatches = stoi(line);
		indexes = (int*)malloc(16 * numPatches * sizeof(int));

		for (patch = 0; patch < numPatches; patch++)
		{
			getline(f, line);
			sscanf(line.c_str(),
				"%d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d",
				&(indexes[i]), &(indexes[i + 1]), &(indexes[i + 2]),
				&(indexes[i + 3]), &(indexes[i + 4]), &(indexes[i + 5]), &(indexes[i + 6]),
				&(indexes[i + 7]), &(indexes[i + 8]), &(indexes[i + 9]), &(indexes[i + 10]),
				&(indexes[i + 11]), &(indexes[i + 12]), &(indexes[i + 13]), &(indexes[i + 14]), &(indexes[i + 15]));
			i += 16;
		}

		getline(f, line);
		numPoints = stoi(line.c_str());

		points = (float*)malloc(numPoints * 3 * sizeof(float));

		for (point = 0; point < numPoints; point++)
		{
			getline(f, line);
			sscanf(line.c_str(), " %f, %f, %f", &(points[j]), &(points[j + 1]), &(points[j + 2]));
			j += 3;
		}
	}

	float m[4][4] = { { -1.0f,  3.0f, -3.0f,  1.0f },
	{ 3.0f, -6.0f,  3.0f, 0.0f },
	{ -3.0f,  3.0f,  0.0f,  0.0f },
	{ 1.0f,  0.0f,  0.0f,  0.0f } };


	float mt[4][4] = { { -1.0f,  3.0f, -3.0f,  1.0f },
	{ 3.0f, -6.0f,  3.0f, 0.0f },
	{ -3.0f,  3.0f,  0.0f,  0.0f },
	{ 1.0f,  0.0f,  0.0f,  0.0f } };


	//para cada patch, calcular a grelha

	vector<vertex> vertices;

	vertex vertexA, vertexB, vertexC, vertexD;

	float nUA[3], nUB[3], nUC[3], nUD[3];
	float nVA[3], nVB[3], nVC[3], nVD[3];

	vertex texA, texB, texC, texD;


	float res[3];

	vertex normalA, normalB, normalC, normalD;

	for (patch = 0; patch < numPatches; patch++)
	{
		int b = 0, c = 0;
		//olhando para os índices indicados na patch, guardar os pontos referentes a esses índices numa matriz P
		for (int a = 0; a < 16; a++)
		{
			int f = patch * 16 + a;
			int index = indexes[f];
			if (((a % 4) == 0) && (a != 0))
			{
				c = 0;
				b++;
			}
			px[c][b] = points[3 * index];
			py[c][b] = points[3 * index + 1];
			pz[c][b] = points[3 * index + 2];
			c++;
		}

		//desenhar patch


		float u, v;
		u = 1.0f / divU;
		v = 1.0f / divV;

		float uu, vv;
		uu = vv = 0.0f;
		for (int i = 0; i < divU; i++)
		{
			float vU[4] = { pow(uu, 3) , pow(uu, 2), uu, 1 };
			
			float vU2[4] = { pow(uu + u, 3), pow(uu + u, 2), uu + u, 1 };

			float dU[4] = { 3* pow(uu, 2) , 2*uu, 1, 0 };

			float dU2[4] = { 3* pow(uu + u, 2), 2*(uu + u), 1, 0 };


			vv = 0;

			for (int j = 0; j < divV; j++)
			{
				float vV[4] = {pow(vv,3), pow(vv, 2), vv, 1 };

				float vV2[4] = { pow(vv + v, 3), pow(vv + v, 2), vv + v, 1 };

				float dV[4] = { 3 * pow(vv, 2) , 2 * vv, 1, 0 };

				float dV2[4] = { 3 * pow(vv + v, 2), 2 * (vv + v), 1, 0 };


				//cálculo do ponto em X

				float *r1 = new float[4];
				float *r2 = new float[4];
				float *r3 = new float[4];


				// Primeiro ponto. Ex: a

				multMatrixVector(*mt, vV, r1);


				multMatrixVector(*px, r1, r2);
				multMatrixVector(*m, r2, r3);
				vertexA.x = multVectorVector(vU, r3);

				multMatrixVector(*py, r1, r2);
				multMatrixVector(*m, r2, r3);
				vertexA.y = multVectorVector(vU, r3);

				multMatrixVector(*pz, r1, r2);
				multMatrixVector(*m, r2, r3);
				vertexA.z = multVectorVector(vU, r3);

				vertices.push_back(vertexA);
				texA.x = uu;
				texA.y = vv;
				texA.z = 0;

				tex->push_back(texA);

				//Normal do 1ºponto de U

				multMatrixVector(*mt, vV, r1);
				multMatrixVector(*px, r1, r2);
				multMatrixVector(*m, r2, r3);
				nUA[0] = multVectorVector(dU, r3);

				multMatrixVector(*py, r1, r2);
				multMatrixVector(*m, r2, r3);
				nUA[1] = multVectorVector(dU, r3);

				multMatrixVector(*pz, r1, r2);
				multMatrixVector(*m, r2, r3);
				nUA[2] = multVectorVector(dU, r3);

				//Normal do 1ºponto de V

				multMatrixVector(*mt, dV, r1);
				
				multMatrixVector(*px, r1, r2);
				multMatrixVector(*m, r2, r3);
				nVA[0] = multVectorVector(vU, r3);

				multMatrixVector(*py, r1, r2);
				multMatrixVector(*m, r2, r3);
				nVA[1] = multVectorVector(vU, r3);

				multMatrixVector(*pz, r1, r2);
				multMatrixVector(*m, r2, r3);
				nVA[2] = multVectorVector(vU, r3);

				normalize(nVA);
				normalize(nUA);
				cross(nUA, nVA, res);
				normalize(res);
				normalA.x = res[0];
				normalA.y = res[1];
				normalA.z = res[2];

				n->push_back(normalA);


				// Segundo ponto. Ex: b

				multMatrixVector(*mt, vV, r1);
				
				multMatrixVector(*px, r1, r2);
				multMatrixVector(*m, r2, r3);
				vertexB.x = multVectorVector(r3, vU2);

				multMatrixVector(*py, r1, r2);
				multMatrixVector(*m, r2, r3);
				vertexB.y = multVectorVector(r3, vU2);

				multMatrixVector(*pz, r1, r2);
				multMatrixVector(*m, r2, r3);
				vertexB.z = multVectorVector(r3, vU2);

				vertices.push_back(vertexB);

				texB.x = uu+u;
				texB.y = vv;
				texB.z = 0;

				tex->push_back(texB);

				
				//Normal do 2ºponto de U

				multMatrixVector(*mt, vV2, r1);

				multMatrixVector(*px, r1, r2);
				multMatrixVector(*m, r2, r3);
				nUB[0] = multVectorVector(dU2, r3);

				multMatrixVector(*py, r1, r2);
				multMatrixVector(*m, r2, r3);
				nUB[1] = multVectorVector(dU2, r3);

				multMatrixVector(*pz, r1, r2);
				multMatrixVector(*m, r2, r3);
				nUB[2] = multVectorVector(dU2, r3);

				//Normal do 2ºponto de V

				multMatrixVector(*mt, dV, r1);

				multMatrixVector(*px, r1, r2);
				multMatrixVector(*m, r2, r3);
				nVB[0] = multVectorVector(vU2, r3);

				multMatrixVector(*py, r1, r2);
				multMatrixVector(*m, r2, r3);
				nVB[1] = multVectorVector(vU2, r3);

				multMatrixVector(*pz, r1, r2);
				multMatrixVector(*m, r2, r3);
				nVB[2] = multVectorVector(vU2, r3);

				normalize(nVB);
				normalize(nUB);
				cross(nUB, nVB, res);
				normalize(res);
				normalB.x = res[0];
				normalB.y = res[1];
				normalB.z = res[2];

				n->push_back(normalB);

				// Primeiro ponto. Ex: c

				multMatrixVector(*mt, vV2, r1);

				multMatrixVector(*px, r1, r2);
				multMatrixVector(*m, r2, r3);
				vertexC.x = multVectorVector(vU, r3);

				multMatrixVector(*py, r1, r2);
				multMatrixVector(*m, r2, r3);
				vertexC.y = multVectorVector(vU, r3);

				multMatrixVector(*pz, r1, r2);
				multMatrixVector(*m, r2, r3);
				vertexC.z = multVectorVector(r3, vU);

				vertices.push_back(vertexC);

				texC.x = uu;
				texC.y = vv + v;
				texC.z = 0;

				tex->push_back(texC);

				//Normal do 3ºponto de U

				multMatrixVector(*mt, vV2, r1);

				multMatrixVector(*px, r1, r2);
				multMatrixVector(*m, r2, r3);
				nUC[0] = multVectorVector(dU, r3);

				multMatrixVector(*py, r1, r2);
				multMatrixVector(*m, r2, r3);
				nUC[1] = multVectorVector(dU, r3);

				multMatrixVector(*pz, r1, r2);
				multMatrixVector(*m, r2, r3);
				nUC[2] = multVectorVector(dU, r3);

				//Normal do 3ºponto de V

				multMatrixVector(*mt, dV2, r1);

				multMatrixVector(*px, r1, r2);
				multMatrixVector(*m, r2, r3);
				nVC[0] = multVectorVector(vU, r3);

				multMatrixVector(*py, r1, r2);
				multMatrixVector(*m, r2, r3);
				nVC[1] = multVectorVector(vU, r3);

				multMatrixVector(*pz, r1, r2);
				multMatrixVector(*m, r2, r3);
				nVC[2] = multVectorVector(vU, r3);

				normalize(nVC);
				normalize(nUC);
				cross(nUC, nVC, res);
				normalize(res);
				normalC.x = res[0];
				normalC.y = res[1];
				normalC.z = res[2];

				n->push_back(normalC);


				vertices.push_back(vertexC);

				tex->push_back(texC);

				n->push_back(normalC);

				vertices.push_back(vertexB);
				
				n->push_back(normalB);

				tex->push_back(texB);

				// Primeiro ponto. Ex: d

				multMatrixVector(*mt, vV2, r1);

				multMatrixVector(*px, r1, r2);
				multMatrixVector(*m, r2, r3);
				vertexD.x = multVectorVector(vU2, r3);

				multMatrixVector(*py, r1, r2);
				multMatrixVector(*m, r2, r3);
				vertexD.y = multVectorVector(vU2, r3);

				multMatrixVector(*pz, r1, r2);
				multMatrixVector(*m, r2, r3);
				vertexD.z = multVectorVector(r3, vU2);

				vertices.push_back(vertexD);

				texD.x = uu+u;
				texD.y = vv+v;
				texD.z = 0;

				tex->push_back(texD);


				//Normal do 4ºponto de U

				multMatrixVector(*mt, vV2, r1);

				multMatrixVector(*px, r1, r2);
				multMatrixVector(*m, r2, r3);
				nUD[0] = multVectorVector(dU2, r3);

				multMatrixVector(*py, r1, r2);
				multMatrixVector(*m, r2, r3);
				nUD[1] = multVectorVector(dU2, r3);

				multMatrixVector(*pz, r1, r2);
				multMatrixVector(*m, r2, r3);
				nUD[2] = multVectorVector(dU2, r3);

				//Normal do 4ºponto de V

				multMatrixVector(*mt, dV2, r1);

				multMatrixVector(*px, r1, r2);
				multMatrixVector(*m, r2, r3);
				nVD[0] = multVectorVector(vU2, r3);

				multMatrixVector(*py, r1, r2);
				multMatrixVector(*m, r2, r3);
				nVD[1] = multVectorVector(vU2, r3);

				multMatrixVector(*pz, r1, r2);
				multMatrixVector(*m, r2, r3);
				nVD[2] = multVectorVector(vU2, r3);

				normalize(nVD);
				normalize(nUD);
				cross(nUD, nVD, res);
				normalize(res);
				normalD.x = res[0];
				normalD.y = res[1];
				normalD.z = res[2];

				n->push_back(normalD);

				// Primeiro ponto. Ex: a
				vertices.push_back(vertexA);
				n->push_back(normalA);
				tex->push_back(texA);

				//ponto C
				vertices.push_back(vertexC);
				n->push_back(normalC);
				tex->push_back(texC);

				//PONTO B
				vertices.push_back(vertexB);
				n->push_back(normalB);
				tex->push_back(texB);

				//ponto b
				vertices.push_back(vertexB);
				n->push_back(normalB);
				tex->push_back(texB);

				//ponto C
				vertices.push_back(vertexC);
				n->push_back(normalC);
				tex->push_back(texC);

				// Primeiro ponto. Ex: d
				vertices.push_back(vertexD);
				n->push_back(normalD);
				tex->push_back(texD);
				
				vv += v;
			}
			uu += u;
		}

	}

	return vertices;
}
/*
Função que calcula os pontos de dois triângulos na face do lado direito da box.
*/
void createTrianglesXRight(vector<vertex> *v, vector<vertex> *n, float x, float yf, float zf, float ry, float rz) 
{
	vertex normal;
	normal.x = 1;
	normal.y = 0;
	normal.z = 0;
	vertex p;
	p.x = x;

	// right hand rule
	// first triangle
	p.y = yf;
	p.z = zf;
	v->push_back(p);
	n->push_back(normal);

	p.y = yf - ry;
	p.z = zf;
	v->push_back(p);
	n->push_back(normal);

	p.y = yf;
	p.z = zf - rz;
	v->push_back(p);
	n->push_back(normal);

	// second triangle
	p.y = yf;
	p.z = zf - rz;
	v->push_back(p);
	n->push_back(normal);

	p.y = yf - ry;
	p.z = zf;
	v->push_back(p);
	n->push_back(normal);

	p.y = yf - ry;
	p.z = zf - rz;
	v->push_back(p);
	n->push_back(normal);
}

/*
Função que calcula os pontos de dois triângulos na face do lado esquerdo da box.
*/
void createTrianglesXLeft(vector<vertex> *v, vector<vertex> *n, float x, float yf, float zf, float ry, float rz) 
{
	vertex normal;
	normal.x = -1;
	normal.y = 0;
	normal.z = 0;
	vertex p;
	p.x = x;

	// first triangle
	p.y = yf;
	p.z = zf;
	v->push_back(p);
	n->push_back(normal);

	p.y = yf - ry;
	p.z = zf;
	v->push_back(p);
	n->push_back(normal);

	p.y = yf;
	p.z = zf + rz;
	v->push_back(p);
	n->push_back(normal);

	// second triangle
	p.y = yf;
	p.z = zf + rz;
	v->push_back(p);
	n->push_back(normal);

	p.y = yf - ry;
	p.z = zf;
	v->push_back(p);
	n->push_back(normal);

	p.y = yf - ry;
	p.z = zf + rz;
	v->push_back(p);
	n->push_back(normal);

	

}

/*
Função que calcula os pontos de dois triângulos na face superior da box.
*/
void createTrianglesYUp(vector<vertex> *v, vector<vertex> *n, float xf, float y, float zf, float rx, float rz) 
{
	vertex normal;
	normal.x = 0;
	normal.y = 1;
	normal.z = 0;

	vertex p;
	p.y = y;

	// right hand rule
	// first triangle
	p.x = xf;
	p.z = zf;
	v->push_back(p);
	n->push_back(normal);

	p.x = xf;
	p.z = zf + rz;
	v->push_back(p);
	n->push_back(normal);

	p.x = xf + rx;
	p.z = zf;
	v->push_back(p);
	n->push_back(normal);

	// second triangle
	p.x = xf + rx;
	p.z = zf;
	v->push_back(p);
	n->push_back(normal);

	p.x = xf ;
	p.z = zf + rz;
	v->push_back(p);
	n->push_back(normal);

	p.x = xf + rx;
	p.z = zf + rz;
	v->push_back(p);
	n->push_back(normal);
}

/*
Função que calcula os pontos de dois triângulos na face inferior da box.
*/
void createTrianglesYDown(vector<vertex> *v, vector<vertex> *n, float xf, float y, float zf, float rx, float rz) 
{
	vertex normal;
	normal.x = 0;
	normal.y = -1;
	normal.z = 0;
	vertex p;
	p.y = y;

	// right hand rule
	// first triangle
	p.x = xf;
	p.z = zf;
	v->push_back(p);
	n->push_back(normal);

	p.x = xf;
	p.z = zf - rz;
	v->push_back(p);
	n->push_back(normal);

	p.x = xf + rx;
	p.z = zf;
	v->push_back(p);
	n->push_back(normal);

	// second triangle
	p.x = xf + rx;
	p.z = zf;
	v->push_back(p);
	n->push_back(normal);

	p.x = xf;
	p.z = zf - rz;
	v->push_back(p);
	n->push_back(normal);

	p.x = xf + rx;
	p.z = zf - rz;
	v->push_back(p);
	n->push_back(normal);
}

/*
Função que calcula os pontos de dois triângulos na face da frente da box.
*/
void createTrianglesZFront(vector<vertex> *v, vector<vertex> *n, float xf, float yf, float z, float rx, float ry) 
{
	vertex normal;
	normal.x = 0;
	normal.y = 0;
	normal.z = 1;
	vertex p;
	p.z = z;

	// right hand rule
	// first triangle
	p.x = xf;
	p.y = yf;
	v->push_back(p);
	n->push_back(normal);

	p.x = xf;
	p.y = yf - ry;
	v->push_back(p);
	n->push_back(normal);

	p.x = xf + rx;
	p.y = yf;
	v->push_back(p);
	n->push_back(normal);

	

	// second triangle
	p.x = xf + rx;
	p.y = yf;
	v->push_back(p);
	n->push_back(normal);

	p.x = xf;
	p.y = yf - ry;
	v->push_back(p);
	n->push_back(normal);

	p.x = xf + rx;
	p.y = yf - ry;
	v->push_back(p);
	n->push_back(normal);
}

/*
Função que calcula os pontos de dois triângulos na face de trás da box.
*/
void createTrianglesZBack(vector<vertex> *v, vector<vertex> *n, float xf, float yf, float z, float rx, float ry)
{
	vertex normal;
	normal.x = 0;
	normal.y = 0;
	normal.z = -1;
	vertex p;
	p.z = z;

	// right hand rule
	// first triangle
	p.x = xf;
	p.y = yf;
	v->push_back(p);
	n->push_back(normal);

	p.x = xf;
	p.y = yf - ry;
	v->push_back(p);
	n->push_back(normal);

	p.x = xf - rx;
	p.y = yf;
	v->push_back(p);
	n->push_back(normal);

	// second triangle
	p.x = xf - rx;
	p.y = yf;
	v->push_back(p);
	n->push_back(normal);

	p.x = xf;
	p.y = yf - ry;
	v->push_back(p);
	n->push_back(normal);

	p.x = xf - rx;
	p.y = yf - ry;
	v->push_back(p);
	n->push_back(normal);

	
}

/*
Função que calcula todos os pontos da face superior da box.
*/
void faceYUp(vector<vertex> *v, vector<vertex> *n, vector<vertex> *tex, float x, float y, float z, float dx, float dz, int d) {
	float rx = dx / d;
	float rz = dz / d;
	float xf;
	float zf;

	vertex t;
	t.z = 0;

	float rt = 1.0f / d;

	int i, j;
	xf = x;
	for (i = 0; i < d; i++)
	{
		zf = z;
		for (j = 0; j < d; j++)
		{
			
			t.x = i * rt;
			t.y = j	* rt;
			tex->push_back(t);

			t.x = i * rt;
			t.y = (j + 1)	* rt;
			tex->push_back(t);

			t.x = (i + 1) * rt;
			t.y = j	* rt;
			tex->push_back(t);


			t.x = (i + 1) * rt;
			t.y = j	* rt;
			tex->push_back(t);

			t.x = i * rt;
			t.y = (j + 1)	* rt;
			tex->push_back(t);

			t.x = (i + 1) * rt;
			t.y = (j + 1)	* rt;
			tex->push_back(t);


			createTrianglesYUp(v, n, xf, y, zf, rx, rz);
			zf += rz;
		}
		xf += rx;
	}
}

/*
Função que calcula todos os pontos da face inferior da box.
*/
void faceYDown(vector<vertex> *v, vector<vertex> *n, vector<vertex> *tex, float x, float y, float z, float dx, float dz, int d) {
	float rx = dx / d;
	float rz = dz / d;
	float xf;
	float zf;

	float rt = 1.0f / d;

	vertex t;
	t.z = 0;

	int i, j;
	xf = x;
	for (i = 0; i < d; i++)
	{
		zf = z;
		for (j = 0; j < d; j++)
		{
			t.x = i * rt;
			t.y = j	* rt;
			tex->push_back(t);

			t.x = i * rt;
			t.y = (j + 1)	* rt;
			tex->push_back(t);

			t.x = (i + 1) * rt;
			t.y = j	* rt;
			tex->push_back(t);


			t.x = (i + 1) * rt;
			t.y = j	* rt;
			tex->push_back(t);

			t.x = i * rt;
			t.y = (j + 1)	* rt;
			tex->push_back(t);

			t.x = (i + 1) * rt;
			t.y = (j + 1)	* rt;
			tex->push_back(t);

			createTrianglesYDown(v, n, xf, y, zf, rx, rz);
			zf -= rz;
		}
		xf += rx;
	}
}

/*
Função que calcula todos os pontos da face de cima da box.
*/
void faceZFront(vector<vertex> *v, vector<vertex> *n, vector<vertex> *tex, float x, float y, float z, float dx, float dy, int d) {
	float rx = dx / d;
	float ry = dy / d;
	float xf;
	float yf;

	yf = y;

	int i, j;

	float rt = 1.0f / d;

	vertex t;
	t.z = 0;

	for (i = 0; i < d; i++)
	{
		xf = x;
		for (j = 0; j < d; j++)
		{
			t.x = j * rt;
			t.y = i	* rt;
			tex->push_back(t);

			t.x = j * rt;
			t.y = (i + 1)	* rt;
			tex->push_back(t);

			t.x = (j + 1) * rt;
			t.y = i	* rt;
			tex->push_back(t);


			t.x = (j + 1) * rt;
			t.y = i	* rt;
			tex->push_back(t);

			t.x = j * rt;
			t.y = (i + 1)	* rt;
			tex->push_back(t);

			t.x = (j + 1) * rt;
			t.y = (i + 1)	* rt;
			tex->push_back(t);

			createTrianglesZFront(v, n, xf, yf, z, rx, ry);

			xf += rx;
		}
		yf -= ry;
	}
}

/*
Função que calcula todos os pontos da face de trás da box.
*/
void faceZBack(vector<vertex> *v, vector<vertex> *n, vector<vertex> *tex, float x, float y, float z, float dx, float dy, int d) {
	float rx = dx / d;
	float ry = dy / d;
	float xf;
	float yf;

	yf = y;

	float rt = 1.0f / d;

	vertex t;
	t.z = 0;

	int i, j;

	for (i = 0; i < d; i++)
	{
		xf = x;
		for (j = 0; j < d; j++)
		{
			t.x = j * rt;
			t.y = i	* rt;
			tex->push_back(t);

			t.x = j * rt;
			t.y = (i + 1)	* rt;
			tex->push_back(t);

			t.x = (j + 1) * rt;
			t.y = i	* rt;
			tex->push_back(t);


			t.x = (j + 1) * rt;
			t.y = i	* rt;
			tex->push_back(t);

			t.x = j * rt;
			t.y = (i + 1)	* rt;
			tex->push_back(t);

			t.x = (j + 1) * rt;
			t.y = (i + 1)	* rt;
			tex->push_back(t);


			createTrianglesZBack(v, n, xf, yf, z, rx, ry);

			xf -= rx;
		}
		yf -= ry; 
	}
}

/*
Função que calcula todos os pontos da face do lado direito da box.
*/
void faceXRight(vector<vertex> *v, vector<vertex> *n, vector<vertex> *tex,  float x, float y, float z, float dy, float dz, int d) {
	float ry = dy / d;
	float rz = dz / d;
	float yf;
	float zf;
	int i, j;

	vertex t;
	t.z = 0;

	float rt = 1.0f / d;

	yf = y;
	for (i = 0; i < d; i++)
	{
		zf = z;
		for (j = 0; j < d; j++)
		{
			t.x = j * rt;
			t.y = i	* rt;
			tex->push_back(t);

			t.x = j * rt;
			t.y = (i + 1)	* rt;
			tex->push_back(t);

			t.x = (j + 1) * rt;
			t.y = i	* rt;
			tex->push_back(t);


			t.x = (j + 1) * rt;
			t.y = i	* rt;
			tex->push_back(t);

			t.x = j * rt;
			t.y = (i + 1)	* rt;
			tex->push_back(t);

			t.x = (j + 1) * rt;
			t.y = (i + 1)	* rt;
			tex->push_back(t);

			createTrianglesXRight(v, n, x, yf, zf, ry, rz);
			zf -= rz;
		}
		yf -= ry;
	}

}

/*
Função que calcula todos os pontos da face do lado esquerdo da box.
*/
void faceXLeft(vector<vertex> *v, vector<vertex> *n, vector<vertex> *tex, float x, float y, float z, float dy, float dz, int d)
{
	
	float ry = dy / d;
	float rz = dz / d;
	float yf;
	float zf;
	int i, j;

	float rt = 1.0f / d;

	vertex t;
	t.z = 0;

	yf = y;
	for (i = 0; i < d; i++)
	{
		zf = z;
		for (j = 0; j < d; j++)
		{
			
			t.x = j * rt;
			t.y = i	* rt;
			tex->push_back(t);

			t.x = j * rt;
			t.y = (i + 1)	* rt;
			tex->push_back(t);

			t.x = (j + 1) * rt;
			t.y = i	* rt;
			tex->push_back(t);


			t.x = (j + 1) * rt;
			t.y = i	* rt;
			tex->push_back(t);

			t.x = j * rt;
			t.y = (i + 1)	* rt;
			tex->push_back(t);

			t.x = (j + 1) * rt;
			t.y = (i + 1)	* rt;
			tex->push_back(t);

			createTrianglesXLeft(v, n, x, yf, zf, ry, rz);
			zf += rz;
		}
		yf -= ry; 
	}
}

/*
Função que calcula todos os pontos da box, utilizando funções auxiliares para cada face. Os pontos foram inicialmente gerados
no octante positivo e depois centrados na origem.
*/
vector<vertex> box(float x, float y, float z, int d, vector<vertex> *n, vector<vertex> *tex) {
	vector<vertex> vertices;

	faceXRight(&vertices, n, tex, x, y, z, y, z, d);
	faceXLeft(&vertices, n, tex, 0, y, 0, y, z, d);

	faceYUp(&vertices, n, tex, 0, y, 0, x, z, d);
	faceYDown(&vertices, n, tex, 0, 0, z, x, z, d);


	faceZFront(&vertices, n, tex, 0, y, z, x, y, d);
	faceZBack(&vertices, n, tex, x, y, 0, x, y, d);

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

vector<vertex> sphere(float radius, int slices, int stacks, vector<vertex> *n, vector<vertex> *tex) {

	vector<vertex> v;
	

	float alpha = 2 * (float)M_PI / slices;
	float rBeta = (float)M_PI / stacks;

	int i, j;

	float beta;

	vertex a, b, c, d;
	vertex na, nb, nc, nd;

	beta = (float)M_PI / 2;

	float rstacks = 1.0f / stacks;
	float rslices = 1.0f / slices;

	vertex t;

	for (i = 0; i < stacks; i++)
	{

		for (j = 0; j < slices; j++)
		{
			a.x = radius * cos(beta) * sin(j * alpha);
			a.y = radius * sin(beta);
			a.z = radius * cos(beta) * cos(j * alpha);

			na.x = cos(beta) * sin(j * alpha);
			na.y = sin(beta);
			na.z = cos(beta) * cos(j * alpha);
			

			b.x = radius * cos(beta) * sin(j * alpha + alpha);
			b.y = radius * sin(beta);
			b.z = radius * cos(beta) * cos(j * alpha + alpha);

			nb.x =  cos(beta) * sin(j * alpha + alpha);
			nb.y =  sin(beta);
			nb.z =  cos(beta) * cos(j * alpha + alpha);


			c.x = radius * cos(beta - rBeta) * sin(j * alpha + alpha);
			c.y = radius * sin(beta - rBeta);
			c.z = radius * cos(beta - rBeta) * cos(j * alpha + alpha);

			nc.x =  cos(beta - rBeta) * sin(j * alpha + alpha);
			nc.y =  sin(beta - rBeta);
			nc.z =  cos(beta - rBeta) * cos(j * alpha + alpha);


			d.x = radius * cos(beta - rBeta) * sin(j * alpha);
			d.y = radius * sin(beta - rBeta);
			d.z = radius * cos(beta - rBeta) * cos(j * alpha);

			nd.x =  cos(beta - rBeta) * sin(j * alpha);
			nd.y =  sin(beta - rBeta);
			nd.z =  cos(beta - rBeta) * cos(j * alpha);

			v.push_back(b);
			v.push_back(a);
			v.push_back(c);

			n->push_back(nb);
			n->push_back(na);
			n->push_back(nc);
			
			t.y = i*rstacks;
			t.x = (j + 1)*rslices;
			t.z = 0;
			
			tex->push_back(t);
			
			t.y = i*rstacks;
			t.x = j*rslices;
			
			tex->push_back(t);

			t.y = (i + 1)*rstacks;
			t.x = (j + 1)*rslices;

			tex->push_back(t);
			
			v.push_back(c);
			v.push_back(a);
			v.push_back(d);

			n->push_back(nc);
			n->push_back(na);
			n->push_back(nd);

			t.y = (i + 1)*rstacks;
			t.x = (j + 1)*rslices;

			tex->push_back(t);

			t.y = i*rstacks;
			t.x = j*rslices;

			tex->push_back(t);

			t.y = (i + 1)*rstacks;
			t.x = j*rslices;

			tex->push_back(t);
		}
		beta -= rBeta;
	}



	return v;
}

/*
Função auxiliar que calcula todos os pontos da base do cone.
*/
vector<vertex> base(vector<vertex> *n, vector<vertex> *tex, float radius, int slices)
{
	vertex normal;
	normal.x = 0;
	normal.y = -1;
	normal.z = 0;

	vertex t;
	t.z = 0;

	vector<vertex> v;
	vertex vertex, vOrigin;
	vertex.y = 0;

	int i;
	float angle;

	float height = 1.0f-0.1875f;
	float width = 0.4375f;


	if (slices == 1)
	{

		angle = (2 * (float)M_PI) / 3;

		for (i = 0; i < 3; i++)
		{
			vertex.x = radius * sin(i* angle);
			vertex.z = radius * cos(i* angle);
			v.push_back(vertex);

			n->push_back(normal);

			t.x = width + 0.1875f*sin(i*angle);
			t.y = height + 0.1875f*cos(i*angle);

			tex->push_back(t);

		}
	}
	else if (slices == 2)
	{

		angle = (float)M_PI / 2;
		for (i = 2; i >= 0; i -= 2)
		{
			vertex.x = radius*sin(i*angle);
			vertex.z = radius * cos(i* angle);
			v.push_back(vertex);
			n->push_back(normal);

			t.x = width + height*sin(i*angle);
			t.y = height + height*cos(i*angle);

			tex->push_back(t);


			vertex.x = radius*sin(i*angle - angle);
			vertex.z = radius * cos(i* angle - angle);
			v.push_back(vertex);
			n->push_back(normal);

			t.x = width + height*sin(i*angle-angle);
			t.y = height + height*cos(i*angle-angle);

			tex->push_back(t);

			vertex.x = radius*sin(i*angle - 2 * angle);
			vertex.z = radius * cos(i* angle - 2 * angle);
			v.push_back(vertex);
			n->push_back(normal);

			t.x = width + height*sin(i*angle - 2 * angle);
			t.y = height + height*cos(i*angle - 2 * angle);

			tex->push_back(t);
		}
	}
	else
	{
		angle = 2 * (float)M_PI / slices;

		vOrigin.x = 0;
		vOrigin.y = 0;
		vOrigin.z = 0;

		for (i = slices; i > 0; i--)
		{
			vertex.x = radius * sin(i * angle + angle);
			vertex.z = radius * cos(i * angle + angle);

			v.push_back(vertex);
			n->push_back(normal);

			t.x = width + 0.1875f*sin(i*angle + angle);
			t.y = height + 0.1875f*cos(i*angle + angle);

			tex->push_back(t);

			vertex.x = radius * sin(i*angle);
			vertex.z = radius * cos(i*angle);

			v.push_back(vertex);
			n->push_back(normal);

			t.x = width + 0.1875f*sin(i*angle);
			t.y = height + 0.1875f*cos(i*angle);

			tex->push_back(t);

			v.push_back(vOrigin);
			n->push_back(normal);

			t.x = width;
			t.y = height;
			tex->push_back(t);

		}
	}

	return v;
}



/*
Função auxiliar que calcula todos os pontos laterais do cone.
*/
vector<vertex> lateral(vector<vertex> v, vector<vertex> *n, vector<vertex> *tex, float radius, float height, int slices, int stacks)
{
	float rh = height / stacks;
	float angle = 2 * (float)M_PI / slices;
	float old_r = (radius * rh) / height;
	int i, j, n_slices = slices;
	float new_r;

	vertex a, b, c, d, normalA, normalB, normalC, normalD;

	float rstacks = 0.625f / stacks;
	float rslices = 1.0f / slices;

	vertex t;

	float istacks = 0.625;

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
	for (i = 1; i < stacks; i++)
	{
		for (j = 0; j < n_slices; j++)
		{
			d.x = old_r *  sin(j*angle);
			d.y = height - i*rh;
			d.z = old_r *  cos(j*angle);

			c.x = old_r * sin(j * angle + angle);
			c.y = height - i * rh;
			c.z = old_r * cos(j * angle + angle);

			new_r = (radius * (/*height -*/ (i + 1) * rh)) / height;

			b.x = new_r * sin(j * angle + angle);
			b.y = height - (i + 1) * rh;
			b.z = new_r * cos(j * angle + angle);

			a.x = new_r *  sin(j*angle);
			a.y = height - (i + 1)*rh;
			a.z = new_r *  cos(j*angle);

			//Normais

			float res[3];

			//a
			
			float vab[3];
			vab[0] = b.x - a.x;
			vab[1] = b.y - a.y;
			vab[2] = b.z - a.z;

			float vad[3];
			vad[0] = d.x - a.x;
			vad[1] = d.y - a.y;
			vad[2] = d.z - a.z;
			
			cross( vab, vad, res);
			normalize(res);
			normalA.x = res[0];
			normalA.y = res[1];
			normalA.z = res[2];

			//b
			float vba[3];
			vba[0] = a.x - b.x;
			vba[1] = a.y - b.y;
			vba[2] = a.z - b.z;

			float vbd[3];
			vbd[0] = d.x - b.x;
			vbd[1] = d.y - b.y;
			vbd[2] = d.z - b.z;

			cross(vbd, vba, res);
			normalize(res);
			normalB.x = res[0];
			normalB.y = res[1];
			normalB.z = res[2];

			//d
			float vda[3];
			vda[0] = a.x - d.x;
			vda[1] = a.y - d.y;
			vda[2] = a.z - d.z;

			float vdb[3];
			vdb[0] = b.x - d.x;
			vdb[1] = b.y - d.y;
			vdb[2] = b.z - d.z;

			cross(vda, vdb, res);
			normalize(res);
			normalD.x = res[0];
			normalD.y = res[1];
			normalD.z = res[2];
			
			
			v.push_back(d);
			v.push_back(a);
			v.push_back(b);
			

			n->push_back(normalD);
			n->push_back(normalA);
			n->push_back(normalB);
			

			
			t.x = j*rslices;
			t.y = i*rstacks;
			t.z = 0;
			tex->push_back(t);
			

			t.x = j * rslices;
			t.y = (i + 1) * rstacks;
			tex->push_back(t);

			t.x = (j + 1)*rslices;
			t.y = ( i+1) * rstacks;
			tex->push_back(t);

			

			//2ºtriângulo
			
			//d
			float vdc[3];
			vdc[0] = c.x - d.x;
			vdc[1] = c.y - d.y;
			vdc[2] = c.z - d.z;

			vdb[0] = b.x - d.x;
			vdb[1] = b.y - d.y;
			vdb[2] = b.z - d.z;

			cross(vdb, vdc, res);
			normalize(res);
			normalD.x = res[0];
			normalD.y = res[1];
			normalD.z = res[2];

			//b
			float vbc[3];
			vbc[0] = c.x - b.x;
			vbc[1] = c.y - b.y;
			vbc[2] = c.z - b.z;

			vbd[0] = d.x - b.x;
			vbd[1] = d.y - b.y;
			vbd[2] = d.z - b.z;


			cross(vbc, vbd, res);
			normalize(res);
			normalB.x = res[0];
			normalB.y = res[1];
			normalB.z = res[2];

			//c
			float vcd[3];
			vcd[0] = d.x - c.x;
			vcd[1] = d.y - c.y;
			vcd[2] = d.z - c.z;

			float vcb[3];
			vcb[0] = b.x - c.x;
			vcb[1] = b.y - c.y;
			vcb[2] = b.z - c.z;

			cross(vcd, vcb, res);
			normalize(res);
			normalC.x = res[0];
			normalC.y = res[1];
			normalC.z = res[2];
			
			v.push_back(d);
			v.push_back(b);
			v.push_back(c);
			

			n->push_back(normalD);
			n->push_back(normalB);
			n->push_back(normalC);
			

			t.x = j*rslices;
			t.y = i*rstacks;

			tex->push_back(t);

			t.x = (j+1)*rslices;
			t.y = (i+1)*rstacks;
		
			tex->push_back(t);

			t.x = (j + 1)*rslices;
			t.y= i * rstacks;

			tex->push_back(t);


		}
		old_r = new_r;
	}

	old_r = (radius * rh) / height;
	c.x = 0;
	c.y = height;
	c.z = 0;
	float res[3];
	
	for (int k = 0; k < n_slices; k++)
	{

		a.x = old_r *  sin(k*angle);
		a.y = height - rh;
		a.z = old_r *  cos(k*angle);

		b.x = old_r * sin(k * angle + angle);
		b.y = height - rh;
		b.z = old_r * cos(k * angle + angle);
		
		//a
		float vab[3];
		vab[0] = b.x - a.x;
		vab[1] = b.y - a.y;
		vab[2] = b.z - a.z;

		float vac[3];
		vac[0] = c.x - a.x;
		vac[1] = c.y - a.y;
		vac[2] = c.z - a.z;

		cross(vab, vac, res);
		normalize(res);
		normalA.x = res[0];
		normalA.y = res[1];
		normalA.z = res[2];

		//b
		float vba[3];
		vba[0] = a.x - b.x;
		vba[1] = a.y - b.y;
		vba[2] = a.z - b.z;

		float vbc[3];
		vbc[0] = c.x - b.x;
		vbc[1] = c.y - b.y;
		vbc[2] = c.z - b.z;

		cross( vbc, vba, res);
		normalize(res);
		normalB.x = res[0];
		normalB.y = res[1];
		normalB.z = res[2];



		//c
		float vca[3];
		vca[0] = a.x - c.x;
		vca[1] = a.y - c.y;
		vca[2] = a.z - c.z;

		float vcb[3];
		vcb[0] = b.x - c.x;
		vcb[1] = b.y - c.y;
		vcb[2] = b.z - c.z;

		cross(vca, vcb, res);
		normalize(res);
		normalC.x = res[0];
		normalC.y = res[1];
		normalC.z = res[2];

		n->push_back(normalA);
		n->push_back(normalB);
		n->push_back(normalC);

		v.push_back(a);
		v.push_back(b);
		v.push_back(c);


		t.x = k*rslices;
		t.y = rstacks;

		tex->push_back(t);

		t.x = (k + 1)*rslices;
		t.y = rstacks;

		tex->push_back(t);

		t.x = k*rslices;
		t.y = 0;

		tex->push_back(t);
	}


	return v;
}

/*
Função auxiliar que calcula todos os pontos do cone.
*/
vector<vertex> cone(vector<vertex> *n, vector<vertex> *tex, float radius, float height, int slices, int stacks)
{
	vector<vertex> v = base(n, tex, radius, slices);

	return lateral(v, n, tex, radius, height, slices, stacks);

}

/*
Função que calcula todos os pontos do anel
*/

vector<vertex> ring(vector<vertex> vertices, vector<vertex> *n, vector<vertex> *tex, float radiusB, float radiusS, int slices)
{
	vertex up, down;
	up.x = down.x = 0;
	up.z = down.z = 0;
	up.y = 1;
	down.y = -1;

	vertex v;
	v.y = 0;

	vertex t;
	t.z = 0;

	float angle = 2 * (float)M_PI / slices;

	float ang1 = 0;
	float ang2 = ang1 + angle;

	float rslices = 1.0f / slices;

	for (int i = 0; i < slices; i++)
	{
		ang2 = ang1 + angle;


		//b
		v.x = radiusB*sin(ang1);
		v.z = radiusB*cos(ang1);

		vertices.push_back(v);
		n->push_back(up);

		t.x = 0;
		t.y = 1;
		tex->push_back(t);

		//d
		v.x = radiusS*sin(ang1);
		v.z = radiusS*cos(ang1);

		vertices.push_back(v);
		n->push_back(up);

		t.x = 1;
		t.y = 1;
		tex->push_back(t);

		//a
		v.x = radiusB*sin(ang2);
		v.z = radiusB*cos(ang2);

		vertices.push_back(v);
		n->push_back(up);
		t.x = 0;
		t.y = 0;
		tex->push_back(t);


		//a
		v.x = radiusB*sin(ang2);
		v.z = radiusB*cos(ang2);

		vertices.push_back(v);
		n->push_back(up);
		t.x = 0;
		t.y = 0;
		tex->push_back(t);

		//d
		v.x = radiusS*sin(ang1);
		v.z = radiusS*cos(ang1);

		vertices.push_back(v);
		n->push_back(up);
		t.x = 1;
		t.y = 1;
		tex->push_back(t);

		//c

		v.x = radiusS*sin(ang2);
		v.z = radiusS*cos(ang2);

		vertices.push_back(v);
		n->push_back(up);
		t.x = 1;
		t.y = 0;
		tex->push_back(t);

		
		//PARTE DE BAIXO

		//d
		v.x = radiusS*sin(ang1);
		v.z = radiusS*cos(ang1);

		vertices.push_back(v);
		n->push_back(down);

		t.x = 1;
		t.y = 1;
		tex->push_back(t);

		//b
		v.x = radiusB*sin(ang1);
		v.z = radiusB*cos(ang1);

		vertices.push_back(v);
		n->push_back(down);

		t.x = 0;
		t.y = 1;
		tex->push_back(t);

		//a
		v.x = radiusB*sin(ang2);
		v.z = radiusB*cos(ang2);

		vertices.push_back(v);
		n->push_back(down);
		t.x = 0;
		t.y = 0;
		tex->push_back(t);


		//d
		v.x = radiusS*sin(ang1);
		v.z = radiusS*cos(ang1);

		vertices.push_back(v);
		n->push_back(down);
		t.x = 1;
		t.y = 1;
		tex->push_back(t);

		//a
		v.x = radiusB*sin(ang2);
		v.z = radiusB*cos(ang2);

		vertices.push_back(v);
		n->push_back(down);
		t.x = 0;
		t.y = 0;
		tex->push_back(t);

		//c

		v.x = radiusS*sin(ang2);
		v.z = radiusS*cos(ang2);

		vertices.push_back(v);
		n->push_back(down);
		t.x = 1;
		t.y = 0;
		tex->push_back(t);

		ang1 += angle;
	}
	return vertices;
}

vector<vertex> createOrbit(float radius, int numPoints)
{
	vector<vertex> vertices;
	float angleSlice = 2 * (float) M_PI / 200;
	float angle = 0;

	vertex v;
	for (int i = 0; i < numPoints ; i++)
	{
		v.x = radius*cos(angle);
		v.y = 0;
		v.z = radius*sin(angle);
		vertices.push_back(v);
		angle += angleSlice;
	}
	return vertices;
}
/*
Função que se encarrega das verificações do input. Se algum parâmetro for incorreto, retorna falso.
*/
bool checkInput(int argc, char** argv) {
	bool ok = true;
	string model = argv[1];
	try {
		if (!(model.compare("teapot") == 0))
		{
			for (int i = 2; i < argc - 1; i++) {
				stof(argv[i]);
			}
		}
		else
		{
			for (int i = 2; i < argc - 2; i++) {
				stof(argv[i]);
			}
		}
	}
	catch (const invalid_argument&) {
		ok = false;
		cout << "The inserted values aren't valid.\n";
	}

	
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
					if (stof(argv[2]) <= 0 || stof(argv[3]) <= 0 || stof(argv[4]) <= 0 || stoi(argv[5]) <= 0)
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
					if (stof(argv[2]) <= 0 || stof(argv[3]) <= 0 || stoi(argv[4]) <= 0)
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
		else if (model.compare("teapot") == 0) {
			if (argc != 6)
			{
				ok = false;
				cout << "Invalid number of arguments.\n";
			}
			else
			{
				try
				{
					if (stoi(argv[2]) <= 0 || stoi(argv[3]) <= 0)
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
		else if (model.compare("orbit") == 0) {
			if (argc != 5)
			{
				ok = false;
				cout << "Invalid number of arguments.\n";
			}
			else
			{
				try
				{
					if (stof(argv[2]) <= 0 || stoi(argv[3]) <= 0)
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

int main(int argc, char **argv)
{
	string model;
	string file;
	
	if (argc >= 2)
	{
		model = argv[1];

		if (checkInput(argc, argv) == true) 
		{
			
			if (model.compare("plane") == 0)
			{
				float coord = stof(argv[2]);
				vector<vertex> n;
				vector<vertex> tex;
				vector<vertex> v = plane(coord, &n, &tex);
				toFile(v, n, tex, argv[argc - 1]);
			}
			else if (model.compare("box") == 0) 
			{
				float x = stof(argv[2]);
				float y = stof(argv[3]);
				float z = stof(argv[4]);
				int div = stoi(argv[5]);
				vector<vertex> n;
				vector<vertex> tex;
				vector<vertex> v = box(x, y, z, div, &n, &tex);
				toFile(v, n, tex, argv[argc - 1]);
			}
			else if (model.compare("sphere") == 0) 
			{
				float radius = stof(argv[2]);
				int slices = stoi(argv[3]);
				int stacks = stoi(argv[4]);

				vector<vertex> n;
				vector<vertex> tex;
				vector<vertex> v = sphere(radius, slices, stacks, &n, &tex);
				toFile(v, n, tex, argv[argc - 1]);
			}
			else if (model.compare("cone") == 0) {
				float radius = stof(argv[2]);
				float height = stof(argv[3]);
				int slices = stoi(argv[4]);
				int stacks = stoi(argv[5]);

				vector<vertex> n;
				vector<vertex> tex;
				vector<vertex> v = cone(&n, &tex, radius, height, slices, stacks);
				toFile(v, n, tex, argv[argc - 1]);
			}
			else if (model.compare("ring") == 0)
			{
				float radiusB = stof(argv[2]);
				float radiusS = stof(argv[3]);
				int slices = stoi(argv[4]);

				vector<vertex> n;
				vector<vertex> v;
				vector<vertex> tex;
				v = ring(v, &n, &tex, radiusB, radiusS, slices);
				toFile(v, n, tex, argv[argc - 1]);
			}
			else if (model.compare("teapot") == 0)
			{
				int divU = stoi(argv[2]);
				int divV = stoi(argv[3]);
				vector<vertex> v;
				vector<vertex> n;
				vector<vertex> tex;

				v = createTeapot(&n, &tex, argv[4], divU, divV);
				toFile(v, n, tex, argv[argc - 1]);
			}
			else if (model.compare("orbit") == 0)
			{
				float radius = stof(argv[2]);
				int numPoints = stoi(argv[3]);
				vector<vertex> v;

				v = createOrbit(radius, numPoints);
				toFile(v, argv[argc - 1]);
			}
		}
	}
	else {
		cout << "Not enough arguments.\n";
	}
	return 0;
}