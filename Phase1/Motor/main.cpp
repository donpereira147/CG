// Generator.cpp : Defines the entry point for the console application.
//
#define _CRT_SECURE_NO_WARNINGS

#include "./tinyxml/tinyxml.h"
#include "./tinyxml/tinystr.h"
#include "Vertex.h"
#include <vector>
#include <fstream>
#include <iostream>
#include <string>
#include <GL/glut.h> 
#include <math.h>


using namespace std;

/*
	Variáveis globais que serão usadas para determinar as rotações dependendo do grau e do eixo selecionado(oX, oY, oZ).
*/
float graus, xx, yy, zz;

/*
	Variáveis que permitem alterar o modo de preenchimento e a cor das prmiitivas gráficas.
*/
int mode, colour;
/*
	Este vector armazena vários vectors que, por sua vez, armazenam os pontos das primitivas gráficas.
*/
vector<vector<vertex>> solids;

/*
	A função load lê o nome do ficheiro XML (que se encontra na pasta conjunta Solids) dado como parâmetro e armazena num vector solid todos os pontos lidos, que
	vai ser armazenado num vector solids, que armazena os pontos das várias primitivas gráficas.
*/
void load(const char* filename)
{
	TiXmlDocument doc(filename);
	if (!doc.LoadFile()) return;
	TiXmlElement* root = doc.RootElement();


	for (TiXmlElement* modelNode = root->FirstChild("model")->ToElement(); modelNode; modelNode = modelNode->NextSiblingElement())
	{

		const char *file = modelNode->Attribute("file");

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

			solids.push_back(solid);

		}
	}
}

void changeSize(int w, int h) 
{
	if (h == 0)
		h = 1;
	float ratio = w * 1.0 / h;
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	glViewport(0, 0, w, h);
	gluPerspective(45.0f, ratio, 1.0f, 1000.0f);
	glMatrixMode(GL_MODELVIEW);
}


/*
	Função que desenha as primitivas gráficas.
*/
void renderScene(void)
{
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	glLoadIdentity();
	gluLookAt(7.0, 5.0, 5.0,
		0.0, 0.0, 0.0,
		0.0f, 1.0f, 0.0f);
	
	if ((mode % 3) == 0)
	{
		glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
	}
	if ((mode % 3) == 1)
	{
		glPolygonMode(GL_FRONT_AND_BACK, GL_POINT);
	}
	if ((mode % 3) == 2)
	{
		glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	}
	glRotatef(graus, xx, yy, zz);

	glEnable(GL_CULL_FACE);
	glFrontFace(GL_CCW);

	if (colour == 0)
	{
		glColor3f(1.0f, 1.0f, 1.0f);
	}
	if (colour == 1)
	{
		glColor3f(1.0f, 0.0f, 0.0f);
	}
	if (colour == 2)
	{
		glColor3f(0.0f, 1.0f, 0.0f);
	}
	if (colour == 3)
	{
		glColor3f(0.0f, 0.0f, 1.0f);
	}
	glBegin(GL_TRIANGLES);
	
	
	for (int i = 0; i < (int) solids.size(); i++)
	{
		vector<vertex> solid = solids[i];

		for (int j = 0; j < (int) solids[i].size(); j++)
		{
			float x = solid[j].x;
			float y = solid[j].y;
			float z = solid[j].z;

			glVertex3f(x, y, z);

		}
	}
	glEnd();
	

	glutSwapBuffers();
}

/*
	Função que altera o modo de apresentação (linha, preenchido, ou por pontos) das primitivas gráficas.
*/
void mouse(int button, int state, int x, int y)
{
	mode++;
	glutPostRedisplay();
}

/*
	Função que altera a cor (vermelho, azul, verde) das primitivas gráficas.
*/
void color(int key_code, int x, int y)
{
	if (key_code == GLUT_KEY_UP)
	{
		colour = 1;
	}
	if (key_code == GLUT_KEY_LEFT)
	{
		colour = 2;
	}
	if (key_code == GLUT_KEY_RIGHT)
	{
		colour = 3;
	}
	if (key_code == GLUT_KEY_DOWN)
	{
		colour = 0;
	}
	glutPostRedisplay();
}

/*
	A função rotate encarrega-se de, dependendo da tecla premida, fazer a rotação da primitiva gráfica em causa. 
*/
void rotate(unsigned char key, int x, int y)
{
	bool posRotate = true;
	if (key == 'w')
	{
		xx = 1;
		yy = 0;
		zz = 0;
	}
		
	if (key == 's')
	{
		xx = 0;
		yy = 1;
		zz = 0;
	}
	if (key == 'x')
	{
		xx = 0;
		yy = 0;
		zz = 1;
	}
	if (key == 'q')
	{
		posRotate = false;
		xx = 1;
		yy = 0;
		zz = 0;
	}

	if (key == 'a')
	{
		posRotate = false;
		xx = 0;
		yy = 1;
		zz = 0;
	}
	if (key == 'z')
	{
		posRotate = false;
		xx = 0;
		yy = 0;
		zz = 1;
	}
	if(posRotate)
	{
		graus += 2;
	}
	else
	{
		graus -= 2;
	}
	glutPostRedisplay();
}


/*
	Função que inicia o programa e recebe como argumento o ficheiro XML que será lido.
*/
int main(int argc, char **argv) {
	 graus = 0;
	 mode = 0;
	 colour = 0;

	string filename("..\\..\\Solids\\");
	filename += argv[1];
	load(filename.c_str());
	 

	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGBA);
	glutInitWindowPosition(0, 0);
	glutInitWindowSize(1000, 1000);
	glutCreateWindow("cg@di");



	glutDisplayFunc(renderScene);
	glutReshapeFunc(changeSize);

	glutKeyboardFunc(rotate);
	glutMouseFunc(mouse);
	glutSpecialFunc(color);


	glEnable(GL_DEPTH_TEST);
	glEnable(GL_CULL_FACE);
	glClearColor(0.0f, 0.0f, 0.0f, 0.0f);

	glutMainLoop();
	
	return 1;
}

