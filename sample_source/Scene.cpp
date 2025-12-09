#include <iostream>
#include <fstream>
#include <cstdio>
#include <cstdlib>
#include <iomanip>
#include <cstring>
#include <string>
#include <vector>
#include <cmath>

#include "tinyxml2.h"
#include "Triangle.h"
#include "Helpers.h"
#include "Scene.h"
#include "Vec4WithColor.h"
#include "Vec4.h"

using namespace tinyxml2;
using namespace std;

/*
	Parses XML file
*/
Scene::Scene(const char *xmlPath)
{
	const char *str;
	XMLDocument xmlDoc;
	XMLElement *xmlElement;

	xmlDoc.LoadFile(xmlPath);

	XMLNode *rootNode = xmlDoc.FirstChild();

	// read background color
	xmlElement = rootNode->FirstChildElement("BackgroundColor");
	str = xmlElement->GetText();
	sscanf(str, "%lf %lf %lf", &backgroundColor.r, &backgroundColor.g, &backgroundColor.b);

	// read culling
	xmlElement = rootNode->FirstChildElement("Culling");
	if (xmlElement != NULL)
	{
		str = xmlElement->GetText();

		if (strcmp(str, "enabled") == 0)
		{
			this->cullingEnabled = true;
		}
		else
		{
			this->cullingEnabled = false;
		}
	}

	// read cameras
	xmlElement = rootNode->FirstChildElement("Cameras");
	XMLElement *camElement = xmlElement->FirstChildElement("Camera");
	XMLElement *camFieldElement;
	while (camElement != NULL)
	{
		Camera *camera = new Camera();

		camElement->QueryIntAttribute("id", &camera->cameraId);

		// read projection type
		str = camElement->Attribute("type");

		if (strcmp(str, "orthographic") == 0)
		{
			camera->projectionType = ORTOGRAPHIC_PROJECTION;
		}
		else
		{
			camera->projectionType = PERSPECTIVE_PROJECTION;
		}

		camFieldElement = camElement->FirstChildElement("Position");
		str = camFieldElement->GetText();
		sscanf(str, "%lf %lf %lf", &camera->position.x, &camera->position.y, &camera->position.z);

		camFieldElement = camElement->FirstChildElement("Gaze");
		str = camFieldElement->GetText();
		sscanf(str, "%lf %lf %lf", &camera->gaze.x, &camera->gaze.y, &camera->gaze.z);

		camFieldElement = camElement->FirstChildElement("Up");
		str = camFieldElement->GetText();
		sscanf(str, "%lf %lf %lf", &camera->v.x, &camera->v.y, &camera->v.z);

		camera->gaze = normalizeVec3(camera->gaze);
		camera->u = crossProductVec3(camera->gaze, camera->v);
		camera->u = normalizeVec3(camera->u);

		camera->w = inverseVec3(camera->gaze);
		camera->v = crossProductVec3(camera->u, camera->gaze);
		camera->v = normalizeVec3(camera->v);

		camFieldElement = camElement->FirstChildElement("ImagePlane");
		str = camFieldElement->GetText();
		sscanf(str, "%lf %lf %lf %lf %lf %lf %d %d",
			   &camera->left, &camera->right, &camera->bottom, &camera->top,
			   &camera->near, &camera->far, &camera->horRes, &camera->verRes);

		camFieldElement = camElement->FirstChildElement("OutputName");
		str = camFieldElement->GetText();
		camera->outputFilename = string(str);

		this->cameras.push_back(camera);

		camElement = camElement->NextSiblingElement("Camera");
	}

	// read vertices
	xmlElement = rootNode->FirstChildElement("Vertices");
	XMLElement *vertexElement = xmlElement->FirstChildElement("Vertex");
	int vertexId = 1;

	while (vertexElement != NULL)
	{
		Vec3WithColor *vertex = new Vec3WithColor();
		vertex->vertexId = vertexId;

		str = vertexElement->Attribute("position");
		sscanf(str, "%lf %lf %lf", &vertex->x, &vertex->y, &vertex->z);

		str = vertexElement->Attribute("color");
		sscanf(str, "%lf %lf %lf", &vertex->color.r, &vertex->color.g, &vertex->color.b);

		this->vertices.push_back(vertex);

		vertexElement = vertexElement->NextSiblingElement("Vertex");

		vertexId++;
	}

	// read translations
	xmlElement = rootNode->FirstChildElement("Translations");
	XMLElement *translationElement = xmlElement->FirstChildElement("Translation");
	while (translationElement != NULL)
	{
		Translation *translation = new Translation();

		translationElement->QueryIntAttribute("id", &translation->translationId);

		str = translationElement->Attribute("value");
		sscanf(str, "%lf %lf %lf", &translation->tx, &translation->ty, &translation->tz);

		this->translations.push_back(translation);

		translationElement = translationElement->NextSiblingElement("Translation");
	}

	// read scalings
	xmlElement = rootNode->FirstChildElement("Scalings");
	XMLElement *scalingElement = xmlElement->FirstChildElement("Scaling");
	while (scalingElement != NULL)
	{
		Scaling *scaling = new Scaling();

		scalingElement->QueryIntAttribute("id", &scaling->scalingId);
		str = scalingElement->Attribute("value");
		sscanf(str, "%lf %lf %lf", &scaling->sx, &scaling->sy, &scaling->sz);

		this->scalings.push_back(scaling);

		scalingElement = scalingElement->NextSiblingElement("Scaling");
	}

	// read rotations
	xmlElement = rootNode->FirstChildElement("Rotations");
	XMLElement *rotationElement = xmlElement->FirstChildElement("Rotation");
	while (rotationElement != NULL)
	{
		Rotation *rotation = new Rotation();

		rotationElement->QueryIntAttribute("id", &rotation->rotationId);
		str = rotationElement->Attribute("value");
		sscanf(str, "%lf %lf %lf %lf", &rotation->angle, &rotation->ux, &rotation->uy, &rotation->uz);

		this->rotations.push_back(rotation);

		rotationElement = rotationElement->NextSiblingElement("Rotation");
	}

	// read meshes
	xmlElement = rootNode->FirstChildElement("Meshes");

	XMLElement *meshElement = xmlElement->FirstChildElement("Mesh");
	while (meshElement != NULL)
	{
		Mesh *mesh = new Mesh();

		meshElement->QueryIntAttribute("id", &mesh->meshId);

		// read mesh faces
		char *row;
		char *cloneStr;
		int vertexId1, vertexId2, vertexId3;
		str = meshElement->GetText();
		cloneStr = strdup(str);

		row = strtok(cloneStr, "\n");
		while (row != NULL)
		{
			int result = sscanf(row, "%d %d %d", &vertexId1, &vertexId2, &vertexId3);

			if (result != EOF)
			{
				Vec3WithColor v1 = *(this->vertices[vertexId1 - 1]);
				Vec3WithColor v2 = *(this->vertices[vertexId2 - 1]);
				Vec3WithColor v3 = *(this->vertices[vertexId3 - 1]);

				mesh->triangles.push_back(Triangle(v1, v2, v3));
			}
			row = strtok(NULL, "\n");
		}
		mesh->numberOfTriangles = mesh->triangles.size();
		this->meshes.push_back(mesh);

		meshElement = meshElement->NextSiblingElement("Mesh");
	}

	// read instances
	xmlElement = rootNode->FirstChildElement("Instances");

	XMLElement *instanceElement = xmlElement->FirstChildElement("Instance");
	while (instanceElement != NULL)
	{
		Instance *instance = new Instance();
		int meshId;

		instanceElement->QueryIntAttribute("id", &instance->instanceId);
		instanceElement->QueryIntAttribute("meshId", &meshId);

		instance->mesh = *(this->meshes[meshId - 1]);

		// read projection type
		str = instanceElement->Attribute("type");

		if (strcmp(str, "wireframe") == 0)
		{
			instance->instanceType = WIREFRAME_INSTANCE;
		}
		else
		{
			instance->instanceType = SOLID_INSTANCE;
		}

		// read instance transformations
		XMLElement *instanceTransformationsElement = instanceElement->FirstChildElement("Transformations");
		XMLElement *instanceTransformationElement = instanceTransformationsElement->FirstChildElement("Transformation");

		while (instanceTransformationElement != NULL)
		{
			char transformationType;
			int transformationId;

			str = instanceTransformationElement->GetText();
			sscanf(str, "%c %d", &transformationType, &transformationId);

			instance->transformationTypes.push_back(transformationType);
			instance->transformationIds.push_back(transformationId);

			instanceTransformationElement = instanceTransformationElement->NextSiblingElement("Transformation");
		}

		instance->numberOfTransformations = instance->transformationIds.size();
		this->instances.push_back(instance);

		instanceElement = instanceElement->NextSiblingElement("Instance");
	}
}

void Scene::assignColorToPixel(int i, int j, Color c)
{
	this->image[i][j].r = c.r;
	this->image[i][j].g = c.g;
	this->image[i][j].b = c.b;
}

/*
	Initializes image with background color
*/
void Scene::initializeImage(Camera *camera)
{
	if (this->image.empty())
	{
		for (int i = 0; i < camera->horRes; i++)
		{
			vector<Color> rowOfColors;
			vector<double> rowOfDepths;

			for (int j = 0; j < camera->verRes; j++)
			{
				rowOfColors.push_back(this->backgroundColor);
				rowOfDepths.push_back(1.01);
			}

			this->image.push_back(rowOfColors);
			this->depth.push_back(rowOfDepths);
		}
	}
	else
	{
		for (int i = 0; i < camera->horRes; i++)
		{
			for (int j = 0; j < camera->verRes; j++)
			{
				assignColorToPixel(i, j, this->backgroundColor);

				this->depth[i][j] = 1.01;
				this->depth[i][j] = 1.01;
				this->depth[i][j] = 1.01;
			}
		}
	}
}

/*
	If given value is less than 0, converts value to 0.
	If given value is more than 255, converts value to 255.
	Otherwise returns value itself.
*/
int Scene::makeBetweenZeroAnd255(double value)
{
	if (value >= 255.0)
		return 255;
	if (value <= 0.0)
		return 0;
	return (int)(value);
}

/*
	Writes contents of image (vector<vector<Color>>) into a PPM file.
*/
void Scene::writeImageToPPMFile(Camera *camera)
{
	ofstream fout;

	fout.open(camera->outputFilename.c_str());

	fout << "P3" << endl;
	fout << "# " << camera->outputFilename << endl;
	fout << camera->horRes << " " << camera->verRes << endl;
	fout << "255" << endl;

	for (int j = camera->verRes - 1; j >= 0; j--)
	{
		for (int i = 0; i < camera->horRes; i++)
		{
			fout << makeBetweenZeroAnd255(this->image[i][j].r) << " "
				 << makeBetweenZeroAnd255(this->image[i][j].g) << " "
				 << makeBetweenZeroAnd255(this->image[i][j].b) << " ";
		}
		fout << endl;
	}
	fout.close();
}

/*
	Transformations, clipping, culling, rasterization are done here.
*/
void Scene::forwardRenderingPipeline(Camera *camera)
{
	// TODO: Implement this function

	//camera matrixi oluştur
	Matrix4 camMatrix;
	Matrix4 projMatrix;
	Matrix4 vpMatrix; 

	vector<Triangle *> newTriangles;
	
	camMatrix = cameraTransformationMatrix(camera);

	//projection matrix 
	if(camera->projectionType == PERSPECTIVE_PROJECTION){
		projMatrix = perspectiveMatrix(camera);
	}
	else{
		projMatrix = orthographicMatrix(camera);
	}


	vpMatrix = viewportMatrix(camera);

	//her mesh için matriceleri oluştur 
	for(int i=0; i < this->instances.size(); i++){
		newTriangles.clear(); //her instance için yeni üçgenler listesi oluştur
		Instance currentInstance = *(this->instances[i]);
		Matrix4 modelMatrix = getIdentityMatrix();

		//model matrixi oluştur
		for(int j=0; j < currentInstance.numberOfTransformations; j++){
			char tType = currentInstance.transformationTypes[j];
			int tId = currentInstance.transformationIds[j];

			Matrix4 transformation_matrix;
			if(tType == 't'){
				Translation *t = this->translations[tId - 1];
				transformation_matrix = translationMatrix(t);
			}
			else if(tType == 's'){
				Scaling *s = this->scalings[tId - 1];
				transformation_matrix = scalingMatrix(s);
			}
			else if(tType == 'r'){
				Rotation *r = this->rotations[tId - 1];
				transformation_matrix = rotationMatrix(r);
			}
			
			modelMatrix = multiplyMatrixWithMatrix(transformation_matrix, modelMatrix);
		}

		Mesh currentMesh = currentInstance.mesh;

		for(int k=0; k < currentMesh.numberOfTriangles; k++){
			Triangle *currentTriangle = new Triangle(currentMesh.triangles[k]);

			Vec4WithColor v1(currentTriangle->v1.x, currentTriangle->v1.y, currentTriangle->v1.z, 1.0, currentTriangle->v1.color);
			Vec4WithColor v2(currentTriangle->v2.x, currentTriangle->v2.y, currentTriangle->v2.z, 1.0, currentTriangle->v2.color);
			Vec4WithColor v3(currentTriangle->v3.x, currentTriangle->v3.y, currentTriangle->v3.z, 1.0, currentTriangle->v3.color);

			//model matrix uygula
			v1 = multiplyMatrixWithVec4WithColor(modelMatrix, v1);
			v2 = multiplyMatrixWithVec4WithColor(modelMatrix, v2);	
			v3 = multiplyMatrixWithVec4WithColor(modelMatrix, v3);

			//view matrix uygula
			v1 = multiplyMatrixWithVec4WithColor(camMatrix, v1);
			v2 = multiplyMatrixWithVec4WithColor(camMatrix, v2);
			v3 = multiplyMatrixWithVec4WithColor(camMatrix, v3);

			//projection matrix uygula
			v1 = multiplyMatrixWithVec4WithColor(projMatrix, v1);
			v2 = multiplyMatrixWithVec4WithColor(projMatrix, v2);
			v3 = multiplyMatrixWithVec4WithColor(projMatrix, v3);

			//vertice perspective divide
			if(v1.t != 0.0) {
				v1.x = v1.x / v1.t;
				v1.y = v1.y / v1.t;
				v1.z = v1.z / v1.t;
				v1.t = 1.0;
			}

			if(v2.t != 0.0){
				v2.x = v2.x / v2.t;
				v2.y = v2.y / v2.t;
				v2.z = v2.z / v2.t;
				v2.t = 1.0;
			}

			if(v3.t != 0.0){
				v3.x = v3.x / v3.t;
				v3.y = v3.y / v3.t;
				v3.z = v3.z / v3.t;
				v3.t = 1.0;
			}

			/*
			//viewport matrix uygula
			v1 = multiplyMatrixWithVec4WithColor(vpMatrix, v1);
			v2 = multiplyMatrixWithVec4WithColor(vpMatrix, v2);
			v3 = multiplyMatrixWithVec4WithColor(vpMatrix, v3);

			//yeni üçgeni ekle, current triangle koordinatları güncelle
			currentTriangle->v1.x = v1.x;
			currentTriangle->v1.y = v1.y;
			currentTriangle->v1.z = v1.z;
			currentTriangle->v1.color = v1.color;
			currentTriangle->v2.x = v2.x;
			currentTriangle->v2.y = v2.y;
			currentTriangle->v2.z = v2.z;
			currentTriangle->v2.color = v2.color;
			currentTriangle->v3.x = v3.x;
			currentTriangle->v3.y = v3.y;
			currentTriangle->v3.z = v3.z;
			currentTriangle->v3.color = v3.color;
			newTriangles.push_back(currentTriangle);
			*/

						
			// --- wireframe clipping in NDC (optional, for WIREFRAME_INSTANCE) ---
			if (currentInstance.instanceType == WIREFRAME_INSTANCE)
			{
				// For each edge, make a copy and clip it
				Vec4WithColor e1v0 = v1;
				Vec4WithColor e1v1 = v2;
				if (lbClipping(e1v0, e1v1))
				{
					// Convert to screen: apply viewport
					e1v0 = multiplyMatrixWithVec4WithColor(vpMatrix, e1v0);
					e1v1 = multiplyMatrixWithVec4WithColor(vpMatrix, e1v1);

					// Now you can draw this edge (midpoint/Bresenham) using e1v0/e1v1
					// e.g. drawLineMidpoint(e1v0, e1v1);
				}

				Vec4WithColor e2v0 = v2;
				Vec4WithColor e2v1 = v3;
				if (lbClipping(e2v0, e2v1))
				{
					e2v0 = multiplyMatrixWithVec4WithColor(vpMatrix, e2v0);
					e2v1 = multiplyMatrixWithVec4WithColor(vpMatrix, e2v1);
					// drawLineMidpoint(e2v0, e2v1);
				}

				Vec4WithColor e3v0 = v3;
				Vec4WithColor e3v1 = v1;
				if (lbClipping(e3v0, e3v1))
				{
					e3v0 = multiplyMatrixWithVec4WithColor(vpMatrix, e3v0);
					e3v1 = multiplyMatrixWithVec4WithColor(vpMatrix, e3v1);
					// drawLineMidpoint(e3v0, e3v1);
				}

				// For wireframe, you *don't* need to push this triangle into newTriangles.
			}
			else
			{
				// SOLID instance: no clipping here (as your HW says)
				// Just do viewport and then store in Triangle / newTriangles as before.

				v1 = multiplyMatrixWithVec4WithColor(vpMatrix, v1);
				v2 = multiplyMatrixWithVec4WithColor(vpMatrix, v2);
				v3 = multiplyMatrixWithVec4WithColor(vpMatrix, v3);

				currentTriangle->v1.x = v1.x; currentTriangle->v1.y = v1.y; currentTriangle->v1.z = v1.z;
				currentTriangle->v1.color = v1.color;
				currentTriangle->v2.x = v2.x; currentTriangle->v2.y = v2.y; currentTriangle->v2.z = v2.z;
				currentTriangle->v2.color = v2.color;
				currentTriangle->v3.x = v3.x; currentTriangle->v3.y = v3.y; currentTriangle->v3.z = v3.z;
				currentTriangle->v3.color = v3.color;

				newTriangles.push_back(currentTriangle);
			}

		}

		//backface culling 
		if(this->cullingEnabled){
			//yeni normalleri hesapla

			for(int i=0; i < newTriangles.size(); i++){
				Triangle *tri = newTriangles[i];

				Vec3 a(tri->v1.x, tri->v1.y, tri->v1.z);
				Vec3 b(tri->v2.x, tri->v2.y, tri->v2.z);
				Vec3 c(tri->v3.x, tri->v3.y, tri->v3.z);

				Vec3 normal = tri->triangleNormal(a, b, c);

				Vec3 v = a;
				double nv = dotProductVec3(normal, v);

				//dot product pozitifse üçgeni sil
				if(nv > 0){
					newTriangles.erase(newTriangles.begin() + i);
					i--; 
				}
			}
		}

/*
		//liang-barsky clipping eğer wireframe ise 
		if(currentInstance.instanceType == WIREFRAME_INSTANCE){
			//her üçgen için
			for(int i=0; i < newTriangles.size(); i++){
				Triangle *tri = newTriangles[i];

			}
		}
		else{
			//solid instance ise clipping yok
			
		
		}
*/

	//rasterization

	//draw triangles

	}

	
}

bool visible(double den, double num, double &tE, double &tL)
{
	double t;
	if(num>0) return false;

	if(den>0)
	{		
		t = num/den;
		if(t>tL) return false;
		if(t>tE) tE = t;
	}
	else
	{
		t = num/den;
		if(t<tE) return false;
		if(t<tL) tL = t;
	}
	
	return true;
}

/*
bool lbClipping(double xmin, double xmax, double ymin, double ymax, 
					double &x0, double &y0, double &x1, double &y1)
{
	double dx = x1 - x0;
	double dy = y1 - y0;

	double tE = 0.0;
	double tL = 1.0;

	//sol
	if(!lbVisible(dx, xmin - x0, tE, tL)) return false;
	//sağ
	if(!lbVisible(-dx, x0 - xmax, tE, tL)) return false;
	//alt
	if(!lbVisible(dy, ymin - y0, tE, tL)) return false;
	//üst
	if(!lbVisible(-dy, y0 - ymax, tE, tL)) return false;

	if(tL < 1.0)
	{
		x1 = x0 + dx * tL;
        y1 = y0 + dy * tL;
	}
	if(tE > 0.0)
	{
        x0 = x0 + dx * tE;
        y0 = y0 + dy * tE;
	}

	return true;
}
*/
bool Scene::lbClipping(Vec4WithColor &v0, Vec4WithColor &v1)
{
    // Parametric range
    double tE = 0.0;   // entering
    double tL = 1.0;   // leaving

    // Copy original endpoints (so we can safely interpolate)
    Vec4WithColor p0 = v0;
    Vec4WithColor p1 = v1;

    double x0 = p0.x, y0 = p0.y, z0 = p0.z;
    double x1 = p1.x, y1 = p1.y, z1 = p1.z;

    double dx = x1 - x0;
    double dy = y1 - y0;
    double dz = z1 - z0;

    // NDC clipping box
    const double xmin = -1.0, xmax = 1.0;
    const double ymin = -1.0, ymax = 1.0;
    const double zmin = -1.0, zmax = 1.0;

    // Liang–Barsky tests 
    if (!visible(dx,   xmin - x0, tE, tL)) return false; // left
    if (!visible(-dx,  x0   - xmax, tE, tL)) return false; // right

    if (!visible(dy,   ymin - y0, tE, tL)) return false; // bottom
    if (!visible(-dy,  y0   - ymax, tE, tL)) return false; // top

    if (!visible(dz,   zmin - z0, tE, tL)) return false; // front
    if (!visible(-dz,  z0   - zmax, tE, tL)) return false; // back

    if (tL < 1.0)
    {
        x1 = x0 + dx * tL;
        y1 = y0 + dy * tL;
        z1 = z0 + dz * tL;
    }

    if (tE > 0.0)
    {
        x0 = x0 + dx * tE;
        y0 = y0 + dy * tE;
        z0 = z0 + dz * tE;
    }

    Color c0 = p0.color;
    Color c1 = p1.color;

    Color newC0 = c0;
    Color newC1 = c1;

    if (tE > 0.0)
    {
        newC0.r = c0.r + tE * (c1.r - c0.r);
        newC0.g = c0.g + tE * (c1.g - c0.g);
        newC0.b = c0.b + tE * (c1.b - c0.b);
    }

    if (tL < 1.0)
    {
        newC1.r = c0.r + tL * (c1.r - c0.r);
        newC1.g = c0.g + tL * (c1.g - c0.g);
        newC1.b = c0.b + tL * (c1.b - c0.b);
    }

    v0.x = x0; v0.y = y0; v0.z = z0; v0.color = newC0;
    v1.x = x1; v1.y = y1; v1.z = z1; v1.color = newC1;

    return true;
}
