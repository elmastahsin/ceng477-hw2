#ifndef _SCENE_H_
#define _SCENE_H_
#include <vector>
#include "Vec3WithColor.h"
#include "Vec4.h"
#include "Vec4WithColor.h" // Bu satırı ekle!
#include "Color.h"
#include "Rotation.h"
#include "Scaling.h"
#include "Translation.h"
#include "Camera.h"
#include "Mesh.h"
#include "Instance.h"
#include "Triangle.h"


class Scene
{
public:
	Color backgroundColor;
	bool cullingEnabled;

	std::vector<std::vector<Color>> image;
	std::vector<std::vector<double>> depth;
	std::vector<Camera *> cameras;
	std::vector<Vec3WithColor *> vertices;
	std::vector<Scaling *> scalings;
	std::vector<Rotation *> rotations;
	std::vector<Translation *> translations;
	std::vector<Mesh *> meshes;
	std::vector<Instance *> instances;

	Scene(const char *xmlPath);

	void assignColorToPixel(int i, int j, Color c);
	void initializeImage(Camera *camera);
	int makeBetweenZeroAnd255(double value);
	void writeImageToPPMFile(Camera *camera);
	void forwardRenderingPipeline(Camera *camera);
	bool lbClipping(Vec4WithColor &v0, Vec4WithColor &v1); // Scene:: kaldırıldı
	void lineRasterizer(Vec4WithColor &p0, Vec4WithColor &p1, Camera *camera);
	void drawTriangleRasterization(Triangle *triangle, Camera *camera);
};

#endif