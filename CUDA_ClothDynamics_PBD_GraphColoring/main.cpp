#include <Windows.h>
#include <stdio.h>

#include "ClothDynamics.h"

int _frame = 0;
int _width = 800;
int _height = 600;
float _zoom = -2.5f;
float _rotx = 0;
float _roty = 0.001f;
float _tx = 0;
float _ty = 0;
int _lastx = 0;
int _lasty = 0;
unsigned char _buttons[3] = { 0 };
bool _simulation = false;
char _FPS_str[100];

Mesh* _cloth_mesh = nullptr, * _obstacle0_mesh = nullptr, * _obstacle1_mesh = nullptr;
ClothDynamics* _cloth = nullptr;

#define SCREEN_CAPTURE

void DrawText(float x, float y, const char* text, void* font = NULL)
{
	glColor3f(0, 0, 0);
	glDisable(GL_DEPTH_TEST);

	glMatrixMode(GL_PROJECTION);
	glPushMatrix();
	glLoadIdentity();
	glOrtho(0.0, (double)_width, 0.0, (double)_height, -1.0, 1.0);

	glMatrixMode(GL_MODELVIEW);
	glPushMatrix();
	glLoadIdentity();

	if (font == NULL) {
		font = GLUT_BITMAP_9_BY_15;
	}

	size_t len = strlen(text);

	glRasterPos2f(x, y);
	for (const char* letter = text; letter < text + len; letter++) {
		if (*letter == '\n') {
			y -= 12.0f;
			glRasterPos2f(x, y);
		}
		glutBitmapCharacter(font, *letter);
	}

	glPopMatrix();
	glMatrixMode(GL_PROJECTION);
	glPopMatrix();
	glMatrixMode(GL_MODELVIEW);
	glEnable(GL_DEPTH_TEST);
}

void Init(void)
{
	glEnable(GL_DEPTH_TEST);
	//_cloth_mesh = new Mesh("../obj/LR_cloth.obj");
	//_cloth_mesh = new Mesh("../obj/MR_cloth.obj");
	_cloth_mesh = new Mesh("../obj/HR_cloth.obj");
	_cloth = new ClothDynamics(_cloth_mesh, 1.0);
}

void FPS(void)
{
	static float framesPerSecond = 0.0f;
	static float lastTime = 0.0f;
	float currentTime = GetTickCount() * 0.001f;
	++framesPerSecond;
	if (currentTime - lastTime > 1.0f) {
		lastTime = currentTime;
		sprintf(_FPS_str, "FPS : %d", (int)framesPerSecond);
		framesPerSecond = 0;
	}
}
void Darw(void)
{
	glEnable(GL_LIGHTING);
	glEnable(GL_LIGHT0);
	char text[100];

	if(_cloth)
		_cloth->draw();

	glDisable(GL_LIGHTING);
	DrawText(10.0f, 580.0f, "Projective Dynamics, ACM TOG 2014");

	glDisable(GL_LIGHTING);
	if (_cloth)
		sprintf(text, "Number of triangles : %d", _cloth_mesh->_fs.size() / 3u);
	DrawText(10.0f, 560.0f, text);
	DrawText(10.0f, 540.0f, _FPS_str);
	sprintf(text, "Frame : %d", _frame);
	DrawText(10.0f, 520.0f, text);
}
void Capture(int endFrame)
{
	if (_frame == 0 || _frame % 2 == 0) {
		static int index = 0;
		char filename[100];
		sprintf_s(filename, "../capture/capture-%d.bmp", index);
		BITMAPFILEHEADER bf;
		BITMAPINFOHEADER bi;
		unsigned char* image = (unsigned char*)malloc(sizeof(unsigned char) * _width * _height * 3);
		FILE* file;
		fopen_s(&file, filename, "wb");
		if (image != NULL) {
			if (file != NULL) {
				glReadPixels(0, 0, _width, _height, 0x80E0, GL_UNSIGNED_BYTE, image);
				memset(&bf, 0, sizeof(bf));
				memset(&bi, 0, sizeof(bi));
				bf.bfType = 'MB';
				bf.bfSize = sizeof(bf) + sizeof(bi) + _width * _height * 3;
				bf.bfOffBits = sizeof(bf) + sizeof(bi);
				bi.biSize = sizeof(bi);
				bi.biWidth = _width;
				bi.biHeight = _height;
				bi.biPlanes = 1;
				bi.biBitCount = 24;
				bi.biSizeImage = _width * _height * 3;
				fwrite(&bf, sizeof(bf), 1, file);
				fwrite(&bi, sizeof(bi), 1, file);
				fwrite(image, sizeof(unsigned char), _height * _width * 3, file);
				fclose(file);
			}
			free(image);
		}
		index++;
		if (index == endFrame) {
			exit(0);
		}
	}
}
void Update(void)
{
	if (_simulation) {
#ifdef SCREEN_CAPTURE
		Capture(617);
#endif
		if (_cloth)
			_cloth->simulation();
		if (_frame == 500) {
		//	exit(0);
		}
		_frame++;
	}
	::glutPostRedisplay();
}
void SetView(double zoom, double tx, double ty, double rotx, double roty)
{
	_zoom = zoom;
	_tx = tx;
	_ty = ty;
	_rotx = rotx;
	_roty = roty;
}
void SetView(void)
{
	//printf("%f, %f, %f, %f, %f\n", _zoom, _tx, _ty, _rotx, _roty);
	//SetView(-4.099999, 0.810000, 0.570000, 0.000000, 0.001000);
	//SetView(-3.2, 0.2, 0.1, 5.0, 12.500999);
	//SetView(-4.800001, 0.55, 0.9, 0.0, 0.0001);
	//SetView(-2.050000, -0.240000, 0.900000, 10.000000, 75.500999);
	//SetView(-2.550000, -0.150000, 0.630000, 12.000000, 68.000999);
	//SetView(-4.050001, -1.110000, 0.570000, 13.000000, 7.500999);

	//----< Self Collision >-------------------
	//SetView(-4.150000, -1.169999, 0.630000, 13.500000, 11.000999);
	//SetView(-4.750002, 0.750000, 0.570000, 11.000000, 222.001007);
	//SetView(-3.450000, -0.300000, 0.240000, 26.000000, 42.501007);

	//SetView(-4.099999, 0.000000, 0.000000, 0.000000, 0.001000);

	//----< Boundary Collision >-------------------
	//SetView(-3.000001, -0.300000, 1.170000, 13.000000, -71.998993);
	//SetView(-2.700003, -0.960000, 1.080000, 12.500000, -13.999001);
	//SetView (-2.549999, -0.120000, 1.589999, 23.500000, -19.999001);
	//SetView(-1.350001, -0.060000, 1.230000, 44.000000, -12.999001);

	//----< Drop Boundary Collision >-----------------
	//SetView(-1.000001, 0.090000, 1.440000, 33.500000, 12.000999);

	//----< Obstacle Collision >-----------------
	//SetView(-1.450001, 0.030000, 1.260000, 35.000000, 43.000999);
	//SetView(-2.000000, 0.030000, 1.080000, 38.500000, -18.999001);
}

void Display(void)
{
	glClearColor(0.8980392156862745f, 0.9490196078431373f, 1.0f, 1.0f);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glEnable(GL_POLYGON_OFFSET_FILL);
	glPolygonOffset(1.1f, 4.0f);
	glLoadIdentity();

	glTranslatef(0, 0, _zoom);
	glTranslatef(_tx, _ty, 0);
	glRotatef(_rotx, 1, 0, 0);
	glRotatef(_roty, 0, 1, 0);

	SetView();

	//glTranslatef(-0.5f, -0.5f, -0.5f);
	Darw();
	FPS();
	glutSwapBuffers();
}

void Reshape(int w, int h)
{
	if (w == 0) {
		h = 1;
	}
	glViewport(0, 0, w, h);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluPerspective(45, (float)w / h, 0.1, 100);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
}

void Motion(int x, int y)
{
	int diffx = x - _lastx;
	int diffy = y - _lasty;
	_lastx = x;
	_lasty = y;

	if (_buttons[2]) {
		_zoom += (float)0.05f * diffx;
	}
	else if (_buttons[0]) {
		_rotx += (float)0.5f * diffy;
		_roty += (float)0.5f * diffx;
	}
	else if (_buttons[1]) {
		_tx += (float)0.03f * diffx;
		_ty -= (float)0.03f * diffy;
	}

	if (_simulation) {
		SetView();
	}
	glutPostRedisplay();
}

void Mouse(int button, int state, int x, int y)
{
	_lastx = x;
	_lasty = y;
	switch (button)
	{
	case GLUT_LEFT_BUTTON:
		_buttons[0] = ((GLUT_DOWN == state) ? 1 : 0);
		break;
	case GLUT_MIDDLE_BUTTON:
		_buttons[1] = ((GLUT_DOWN == state) ? 1 : 0);
		break;
	case GLUT_RIGHT_BUTTON:
		_buttons[2] = ((GLUT_DOWN == state) ? 1 : 0);
		break;
	default:
		break;
	}
	glutPostRedisplay();
}

void SpecialInput(int key, int x, int y)
{
	glutPostRedisplay();
}

void Keyboard(unsigned char key, int x, int y)
{
	switch (key)
	{
	case 'q':
	case 'Q':
		exit(0);
	case ' ':
		_simulation = !_simulation;
		break;
	case 'r':
	case 'R':
		_cloth->reset();
		break;
	case 'c':
	case 'C':
		printf("%f, %f, %f, %f, %f\n", _zoom, _tx, _ty, _rotx, _roty);
		break;
	case 'd':
	case 'D':
		_cloth->_bvh->_test++;
		//CollisionSolver::Debug();
		break;
	}
	glutPostRedisplay();
}

int main(int argc, char** argv)
{

	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA | GLUT_DEPTH);
	glutInitWindowSize(_width, _height);
	glutInitWindowPosition(100, 100);
	glutCreateWindow("Cloth Simulator");
	glutDisplayFunc(Display);
	glutReshapeFunc(Reshape);
	glutIdleFunc(Update);
	glutMouseFunc(Mouse);
	glutMotionFunc(Motion);
	glutKeyboardFunc(Keyboard);
	glutSpecialFunc(SpecialInput);
	Init();
	glutMainLoop();
}