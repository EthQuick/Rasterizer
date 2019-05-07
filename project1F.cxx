/* Ethan Quick | CIS 441 | Winter 2019 | Project 1F
 * This file contains starter code provided in class
 * though the bulk of it is code written by me to make it functional.
 */
#include <iostream>
#include <vtkDataSet.h>
#include <vtkImageData.h>
#include <vtkPNGWriter.h>
#include <vtkPointData.h>
#include <vtkPolyData.h>
#include <vtkPolyDataReader.h>
#include <vtkPoints.h>
#include <vtkUnsignedCharArray.h>
#include <vtkFloatArray.h>
#include <vtkCellArray.h>
#include <vtkDoubleArray.h>
#include <math.h>
#include <string>
#include <stdio.h>

#define NORMALS

using std::cerr;
using std::endl;

enum T_Type { TOP, BOTTOM, OTHER};
int PRI_INT = 0;
int tri_id = 0;

struct LightingParameters
{
    LightingParameters(void)
    {
         lightDir[0] = -0.6;
         lightDir[1] = 0;
         lightDir[2] = -0.8;
         Ka = 0.3;
         Kd = 0.7;
         Ks = 2.3;
         alpha = 2.5;
    };
  

    double lightDir[3]; // The direction of the light source
    double Ka;           // The coefficient for ambient lighting.
    double Kd;           // The coefficient for diffuse lighting.
    double Ks;           // The coefficient for specular lighting.
    double alpha;        // The exponent term for specular lighting.
};

LightingParameters lp; //make it global so I dont need to pass it around

double ceil_441(double f)
{
    return ceil(f-0.00001);
}

double floor_441(double f)
{
    return floor(f+0.00001);
}

class Matrix
{
  public:
    double          A[4][4] = {0};  // A[i][j] means row i, column j

    void            TransformPoint(const double *ptIn, double *ptOut);
    static Matrix   ComposeMatrices(const Matrix &, const Matrix &);
    void            Print(ostream &o);
};

void
Matrix::Print(ostream &o)
{
    for (int i = 0 ; i < 4 ; i++)
    {
        char str[256];
        sprintf(str, "(%.7f %.7f %.7f %.7f)\n", A[i][0], A[i][1], A[i][2], A[i][3]);
        o << str;
    }
}

Matrix
Matrix::ComposeMatrices(const Matrix &M1, const Matrix &M2)
{
    Matrix rv;
    for (int i = 0 ; i < 4 ; i++)
        for (int j = 0 ; j < 4 ; j++)
        {
            rv.A[i][j] = 0;
            for (int k = 0 ; k < 4 ; k++)
                rv.A[i][j] += M1.A[i][k]*M2.A[k][j];
        }

    return rv;
}

void
Matrix::TransformPoint(const double *ptIn, double *ptOut)
{
    ptOut[0] = ptIn[0]*A[0][0]
             + ptIn[1]*A[1][0]
             + ptIn[2]*A[2][0]
             + ptIn[3]*A[3][0];
    ptOut[1] = ptIn[0]*A[0][1]
             + ptIn[1]*A[1][1]
             + ptIn[2]*A[2][1]
             + ptIn[3]*A[3][1];
    ptOut[2] = ptIn[0]*A[0][2]
             + ptIn[1]*A[1][2]
             + ptIn[2]*A[2][2]
             + ptIn[3]*A[3][2];
    ptOut[3] = ptIn[0]*A[0][3]
             + ptIn[1]*A[1][3]
             + ptIn[2]*A[2][3]
             + ptIn[3]*A[3][3];
}

class Camera
{
  public:
    double          near, far;
    double          angle;
    double          position[3];
    double          focus[3];
    double          up[3];

    Matrix          ViewTransform() {
                        Matrix res;
                        double temp1 = 1/tan(angle/2);
                        
                        res.A[0][0] = temp1;
                        res.A[1][1] = temp1;
                        res.A[2][2] = (far+near)/(far-near);
                        res.A[3][2] = (2*far*near)/(far - near);
                        res.A[2][3] = -1;

                        return res;
                    }

    Matrix          CameraTransform() {
                        //create the camera frame
                        Matrix cf; //camera frame
                        //create w
                        cf.A[0][2] = position[0] - focus[0];
                        cf.A[1][2] = position[1] - focus[1];
                        cf.A[2][2] = position[2] - focus[2];
                        //create u
                        cf.A[0][0] = up[1]*cf.A[2][2] - up[2]*cf.A[1][2];
                        cf.A[1][0] = up[2]*cf.A[0][2] - up[0]*cf.A[2][2];
                        cf.A[2][0] = up[0]*cf.A[1][2] - up[1]*cf.A[0][2];
                        //create v
                        cf.A[0][1] = cf.A[1][2]*cf.A[2][0] - cf.A[2][2]*cf.A[1][0];
                        cf.A[1][1] = cf.A[2][2]*cf.A[0][0] - cf.A[0][2]*cf.A[2][0];
                        cf.A[2][1] = cf.A[0][2]*cf.A[1][0] - cf.A[1][2]*cf.A[0][0];
                        //magnitude of u
                        double magu = sqrt(cf.A[0][0]*cf.A[0][0] + cf.A[1][0]*cf.A[1][0] + cf.A[2][0]*cf.A[2][0]); 
                        //for v                        
                        double magv = sqrt(cf.A[0][1]*cf.A[0][1] + cf.A[1][1]*cf.A[1][1] + cf.A[2][1]*cf.A[2][1]);
                        //for w
                        double magw = sqrt(cf.A[0][2]*cf.A[0][2] + cf.A[1][2]*cf.A[1][2] + cf.A[2][2]*cf.A[2][2]);
                        //normalize
                        for(int i = 0; i < 3; i++){
                            cf.A[i][0] = cf.A[i][0] / magu;
                            cf.A[i][1] = cf.A[i][1] / magv;
                            cf.A[i][2] = cf.A[i][2] / magw;
                        }
                        double O[3] = {0};
                        double t[3];
                        for(int i = 0; i < 3; i++){
                            t[i] = O[i] - position[i];
                        }
                        //dot the last row of cf with t
                        cf.A[3][0] = cf.A[0][0]*t[0] + cf.A[1][0]*t[1] + cf.A[2][0]*t[2]; 
                        cf.A[3][1] = cf.A[0][1]*t[0] + cf.A[1][1]*t[1] + cf.A[2][1]*t[2];
                        cf.A[3][2] = cf.A[0][2]*t[0] + cf.A[1][2]*t[1] + cf.A[2][2]*t[2];
                        cf.A[3][3] = 1;
                        
                        return cf;
                    }
    Matrix          DeviceTransform(double n, double m) {
                        Matrix res;
                        double temp1 = n/2;
                        double temp2 = m/2;
                        res.A[0][0] = temp1;
                        res.A[1][1] = temp2;
                        res.A[2][2] = 1;
                        res.A[3][3] = 1;
                        res.A[3][0] = temp1;
                        res.A[3][1] = temp2;

                        return res;
                    }
};

double SineParameterize(int curFrame, int nFrames, int ramp)
{
    int nNonRamp = nFrames-2*ramp;
    double height = 1./(nNonRamp + 4*ramp/M_PI);
    if (curFrame < ramp)
    {
        double factor = 2*height*ramp/M_PI;
        double eval = cos(M_PI/2*((double)curFrame)/ramp);
        return (1.-eval)*factor;
    }
    else if (curFrame > nFrames-ramp)
    {
        int amount_left = nFrames-curFrame;
        double factor = 2*height*ramp/M_PI;
        double eval =cos(M_PI/2*((double)amount_left/ramp));
        return 1. - (1-eval)*factor;
    }
    double amount_in_quad = ((double)curFrame-ramp);
    double quad_part = amount_in_quad*height;
    double curve_part = height*(2*ramp)/M_PI;
    return quad_part+curve_part;
}

Camera
GetCamera(int frame, int nframes)
{
    double t = SineParameterize(frame, nframes, nframes/10);
    Camera c;
    c.near = 5;
    c.far = 200;
    c.angle = M_PI/6;
    c.position[0] = 40*sin(2*M_PI*t);
    c.position[1] = 40*cos(2*M_PI*t);
    c.position[2] = 40;
    c.focus[0] = 0;
    c.focus[1] = 0;
    c.focus[2] = 0;
    c.up[0] = 0;
    c.up[1] = 1;
    c.up[2] = 0;
    return c;
}

vtkImageData *
NewImage(int width, int height)
{
    vtkImageData *img = vtkImageData::New();
    img->SetDimensions(width, height, 1);
    img->AllocateScalars(VTK_UNSIGNED_CHAR, 3);

    return img;
}

void
WriteImage(vtkImageData *img, const char *filename)
{
   std::string full_filename = filename;
   full_filename += ".png";
   vtkPNGWriter *writer = vtkPNGWriter::New();
   writer->SetInputData(img);
   writer->SetFileName(full_filename.c_str());
   writer->Write();
   writer->Delete();
}

class Triangle
{
  public:
      double         X[3];
      double         Y[3];
      double         Z[3];
      int            V1 = 0;
      int            V2 = 1;
      int            V3 = 2;
      double         colors[3][3];
      double         normals[3][3]; //[vertexnumber][coordinate]
      double         shading[3];

  // would some methods for transforming the triangle in place be helpful?
        T_Type sortvertices(){
            //sort the vertices in order of decreasing y
            T_Type type = OTHER;
            int temp;
            if(Y[V1] < Y[V3]){
                std::swap(V1, V3);
            }
            if(Y[V2] < Y[V3]){
                std::swap(V2, V3);
            }
            if(Y[V1] < Y[V2]){
                std::swap(V1, V2);
            }
            //flat top - V1 and V2 are on the top
            if(Y[V1] == Y[V2]){
                if(X[V1] > X[V2]){
                    std::swap(V1, V2);
                }
                type = TOP;
            }
            //flat bottom - V1 and V2 are on the bottom
            else if(Y[V2] == Y[V3]){
                std::swap(V1, V3);
                if(X[V1] > X[V2]){
                   std::swap(V1, V2); 
                }
                type = BOTTOM;
            }
            return type;
        }
};

class Screen
{
  public:
      unsigned char   *buffer;
      double          *depthBuff;
      int width, height;

  // would some methods for accessing and setting pixels be helpful?

    //color pixel r, c with a given color
    bool colorpix(int r, int c, double color[3], double newZ, double curShad){
        if(r >= height || r < 0 || c >= width || c < 0){
            return false;
        }
        //check the z
        if(!(newZ > depthBuff[(r*width) + c])){
            return false;
        }
        depthBuff[(r*width) + c] = newZ;
        //draw the pixel in buffer
        int index = (r*width + c)*3;
        for(int i = 0; i < 3; i++){
            buffer[index+i] = ceil_441(std::min(color[i]*curShad, 1.0)*255);
        }
        return true;
    }

    void ClearScreen(){
        int npixels = width*height;
        for (int i = 0 ; i < npixels*3 ; i++)
            buffer[i] = 0;
        for(int i = 0; i<npixels; i++)
            depthBuff[i] = -1;
    }
};

std::vector<Triangle>
GetTriangles(void)
{
    vtkPolyDataReader *rdr = vtkPolyDataReader::New();
    rdr->SetFileName("proj1e_geometry.vtk");
    cerr << "Reading" << endl;
    rdr->Update();
    cerr << "Done reading" << endl;
    if (rdr->GetOutput()->GetNumberOfCells() == 0)
    {
        cerr << "Unable to open file!!" << endl;
        exit(EXIT_FAILURE);
    }
    vtkPolyData *pd = rdr->GetOutput();

    int numTris = pd->GetNumberOfCells();
    vtkPoints *pts = pd->GetPoints();
    vtkCellArray *cells = pd->GetPolys();
    vtkDoubleArray *var = (vtkDoubleArray *) pd->GetPointData()->GetArray("hardyglobal");
    double *color_ptr = var->GetPointer(0);
    //vtkFloatArray *var = (vtkFloatArray *) pd->GetPointData()->GetArray("hardyglobal");
    //float *color_ptr = var->GetPointer(0);
    vtkFloatArray *n = (vtkFloatArray *) pd->GetPointData()->GetNormals();
    float *normals = n->GetPointer(0);
    std::vector<Triangle> tris(numTris);
    vtkIdType npts;
    vtkIdType *ptIds;
    int idx;
    for (idx = 0, cells->InitTraversal() ; cells->GetNextCell(npts, ptIds) ; idx++)
    {
        if (npts != 3)
        {
            cerr << "Non-triangles!! ???" << endl;
            exit(EXIT_FAILURE);
        }
        double *pt = NULL;
        pt = pts->GetPoint(ptIds[0]);
        tris[idx].X[0] = pt[0];
        tris[idx].Y[0] = pt[1];
        tris[idx].Z[0] = pt[2];
#ifdef NORMALS
        tris[idx].normals[0][0] = normals[3*ptIds[0]+0];
        tris[idx].normals[0][1] = normals[3*ptIds[0]+1];
        tris[idx].normals[0][2] = normals[3*ptIds[0]+2];
#endif
        pt = pts->GetPoint(ptIds[1]);
        tris[idx].X[1] = pt[0];
        tris[idx].Y[1] = pt[1];
        tris[idx].Z[1] = pt[2];
#ifdef NORMALS
        tris[idx].normals[1][0] = normals[3*ptIds[1]+0];
        tris[idx].normals[1][1] = normals[3*ptIds[1]+1];
        tris[idx].normals[1][2] = normals[3*ptIds[1]+2];
#endif
        pt = pts->GetPoint(ptIds[2]);
        tris[idx].X[2] = pt[0];
        tris[idx].Y[2] = pt[1];
        tris[idx].Z[2] = pt[2];
#ifdef NORMALS
        tris[idx].normals[2][0] = normals[3*ptIds[2]+0];
        tris[idx].normals[2][1] = normals[3*ptIds[2]+1];
        tris[idx].normals[2][2] = normals[3*ptIds[2]+2];
#endif

        // 1->2 interpolate between light blue, dark blue
        // 2->2.5 interpolate between dark blue, cyan
        // 2.5->3 interpolate between cyan, green
        // 3->3.5 interpolate between green, yellow
        // 3.5->4 interpolate between yellow, orange
        // 4->5 interpolate between orange, brick
        // 5->6 interpolate between brick, salmon
        double mins[7] = { 1, 2, 2.5, 3, 3.5, 4, 5 };
        double maxs[7] = { 2, 2.5, 3, 3.5, 4, 5, 6 };
        unsigned char RGB[8][3] = { { 71, 71, 219 }, 
                                    { 0, 0, 91 },
                                    { 0, 255, 255 },
                                    { 0, 128, 0 },
                                    { 255, 255, 0 },
                                    { 255, 96, 0 },
                                    { 107, 0, 0 },
                                    { 224, 76, 76 } 
                                  };
        for (int j = 0 ; j < 3 ; j++)
        {
            float val = color_ptr[ptIds[j]];
            int r;
            for (r = 0 ; r < 7 ; r++)
            {
                if (mins[r] <= val && val < maxs[r])
                    break;
            }
            if (r == 7)
            {
                cerr << "Could not interpolate color for " << val << endl;
                exit(EXIT_FAILURE);
            }
            double proportion = (val-mins[r]) / (maxs[r]-mins[r]);
            tris[idx].colors[j][0] = (RGB[r][0]+proportion*(RGB[r+1][0]-RGB[r][0]))/255.0;
            tris[idx].colors[j][1] = (RGB[r][1]+proportion*(RGB[r+1][1]-RGB[r][1]))/255.0;
            tris[idx].colors[j][2] = (RGB[r][2]+proportion*(RGB[r+1][2]-RGB[r][2]))/255.0;
        }
    }

    return tris;
}

double calc_shading(double *normal, double *viewDirection){
    //normalize viewDirection
    double mag = viewDirection[0]*viewDirection[0] + viewDirection[1]*viewDirection[1] + viewDirection[2]*viewDirection[2];
    mag = sqrt(mag);
    for(int i = 0; i < 3; i++){
        viewDirection[i] = viewDirection[i] / mag;
    }
    //calc diffuse
    double diff = lp.Kd * abs(normal[0]*lp.lightDir[0] + normal[1]*lp.lightDir[1] + normal[2]*lp.lightDir[2]);
    //calc R 
    double LN = 2 * (lp.lightDir[0]*normal[0] + lp.lightDir[1]*normal[1] + lp.lightDir[2]*normal[2]);
    double R[3] = { (LN * normal[0]) - lp.lightDir[0],
                    (LN * normal[1]) - lp.lightDir[1],
                    (LN * normal[2]) - lp.lightDir[2]};
    //normalize R
    double mag1 = R[0]*R[0] + R[1]*R[1] + R[2]*R[2];
    mag1 = sqrt(mag1);
    for(int i = 0; i < 3; i++){
        R[i] = R[i] / mag1;
    }
    //calc specular
    double spec = viewDirection[0]*R[0] + viewDirection[1]*R[1] + viewDirection[2]*R[2]; 
    if(spec < 0){
        spec = 0;
    }
    spec = lp.Ks * pow(spec, lp.alpha);
    double res = lp.Ka + diff + spec;
    return res; 
}

void shade_triangle(Triangle *t, Camera *c){
    //calc the shading
    double viewDirection[3];
    viewDirection[0] = c->position[0] - t->X[0];
    viewDirection[1] = c->position[1] - t->Y[0];
    viewDirection[2] = c->position[2] - t->Z[0];
    t->shading[0] = calc_shading(t->normals[0], viewDirection); //Shading for 0
    viewDirection[0] = c->position[0] - t->X[1];
    viewDirection[1] = c->position[1] - t->Y[1];
    viewDirection[2] = c->position[2] - t->Z[1];
    t->shading[1] = calc_shading(t->normals[1], viewDirection); //Shading for 1
    viewDirection[0] = c->position[0] - t->X[2];
    viewDirection[1] = c->position[1] - t->Y[2];
    viewDirection[2] = c->position[2] - t->Z[2];
    t->shading[2] = calc_shading(t->normals[2], viewDirection); //Shading for 2
}
double interpolate(double A, double B, double X, double fA, double fB){
    //A and B are the end points and X is the midpoint
    //fA and fB are the values at A and B
    //this function should return F(X)
    /*
    if(A > B){
        std::swap(A, B);
        std::swap(fA, fB);
    }*/
    double t = (X - A) / (B - A);
    double fX = fA + t*(fB - fA);
    return fX;
}

void rasterize_triangle(Triangle t, Screen *screen, T_Type type){
        //find rowMin and rowMax
        double rowMin;
        double rowMax;
        if(type == BOTTOM){
            rowMin = ceil_441(t.Y[t.V1]);
            rowMax = floor_441(t.Y[t.V3]);
        }
        else{
            rowMin = ceil_441(t.Y[t.V3]);
            rowMax = floor_441(t.Y[t.V1]);
        }
        
        //go through the rows
        for(double r = rowMin; r <= rowMax; r++){
            double leftEnd;
            double rightEnd;
            //check if right triangle one left side
            if(t.X[t.V1] == t.X[t.V3]){
                leftEnd = t.X[t.V1];
            }
            //if not find with math
            else{
                double dy = t.Y[t.V3] - t.Y[t.V1];
                double dx = t.X[t.V3] - t.X[t.V1];
                double m = dy/dx;
                double b = t.Y[t.V1] - (t.X[t.V1] * m);
                leftEnd = (b + r*-1)/m*-1;
            }
            //check if right triangle on right side
            if(t.X[t.V2] == t.X[t.V3]){
                rightEnd = t.X[t.V2];
            }
            //if not find with math
            else{
                double dy = t.Y[t.V2] - t.Y[t.V3];
                double dx = t.X[t.V2] - t.X[t.V3];
                double m = dy/dx;
                double b = t.Y[t.V2] - (t.X[t.V2] * m);
                rightEnd = (b + r*-1)/m*-1;
            }
            //Calculate colorL and colorR
            //colorL
            double colorL[3];
            colorL[0] = interpolate(t.Y[t.V1], t.Y[t.V3], r, t.colors[t.V1][0], t.colors[t.V3][0]);
            colorL[1] = interpolate(t.Y[t.V1], t.Y[t.V3], r, t.colors[t.V1][1], t.colors[t.V3][1]);
            colorL[2] = interpolate(t.Y[t.V1], t.Y[t.V3], r, t.colors[t.V1][2], t.colors[t.V3][2]);
            //colorR
            double colorR[3];
            colorR[0] = interpolate(t.Y[t.V2], t.Y[t.V3], r, t.colors[t.V2][0], t.colors[t.V3][0]);
            colorR[1] = interpolate(t.Y[t.V2], t.Y[t.V3], r, t.colors[t.V2][1], t.colors[t.V3][1]);
            colorR[2] = interpolate(t.Y[t.V2], t.Y[t.V3], r, t.colors[t.V2][2], t.colors[t.V3][2]);
            //leftZ
            double leftZ = interpolate(t.Y[t.V1], t.Y[t.V3], r, t.Z[t.V1], t.Z[t.V3]);
            //rightZ
            double rightZ = interpolate(t.Y[t.V2], t.Y[t.V3], r, t.Z[t.V2], t.Z[t.V3]);
            //interpolate the shading
            double leftShad = interpolate(t.Y[t.V1], t.Y[t.V3], r, t.shading[t.V1], t.shading[t.V3]);
            double rightShad = interpolate(t.Y[t.V2], t.Y[t.V3], r, t.shading[t.V2], t.shading[t.V3]);
            //paint the pixels
            for(double c = ceil_441(leftEnd); c <= floor_441(rightEnd); c++){
                //calc color using the colors at leftend and right end
                double curcolor[3];
                curcolor[0] = interpolate(leftEnd, rightEnd, c, colorL[0], colorR[0]);
                curcolor[1] = interpolate(leftEnd, rightEnd, c, colorL[1], colorR[1]);
                curcolor[2] = interpolate(leftEnd, rightEnd, c, colorL[2], colorR[2]);
                //newZ
                double newZ = interpolate(leftEnd, rightEnd, c, leftZ, rightZ);
                //final shading
                double curShad = interpolate(leftEnd, rightEnd, c, leftShad, rightShad);
                screen->colorpix(r, c, curcolor, newZ, curShad);
            }
    }
}

void rasterize_general(Triangle t, Screen *screen){
        //cout << tri_id << endl;
        //tri_id++;
        //sort the vertices
        T_Type type = t.sortvertices();
        switch(type){
            case OTHER:{
                Triangle temp1; //V1 V2 and the other
                Triangle temp2; //V3 V2 and the other
                temp1.X[0] = t.X[t.V1];
                temp1.X[1] = t.X[t.V2];
                temp1.Y[0] = t.Y[t.V1];
                temp1.Y[1] = t.Y[t.V2];
                temp1.Y[2] = t.Y[t.V2];
                temp1.Z[0] = t.Z[t.V1];
                temp1.Z[1] = t.Z[t.V2];
                temp1.shading[0] = t.shading[t.V1];
                temp1.shading[1] = t.shading[t.V2];

                temp2.X[0] = t.X[t.V3];
                temp2.X[1] = t.X[t.V2];
                temp2.Y[0] = t.Y[t.V3];
                temp2.Y[1] = t.Y[t.V2];
                temp2.Y[2] = t.Y[t.V2];
                temp2.Z[0] = t.Z[t.V3];
                temp2.Z[1] = t.Z[t.V2];
                temp2.shading[0] = t.shading[t.V3];
                temp2.shading[1] = t.shading[t.V2];

                temp1.colors[0][0] = t.colors[t.V1][0];
                temp1.colors[0][1] = t.colors[t.V1][1];
                temp1.colors[0][2] = t.colors[t.V1][2];
                temp1.colors[1][0] = t.colors[t.V2][0];
                temp1.colors[1][1] = t.colors[t.V2][1];
                temp1.colors[1][2] = t.colors[t.V2][2];

                temp2.colors[0][0] = t.colors[t.V3][0];
                temp2.colors[0][1] = t.colors[t.V3][1];
                temp2.colors[0][2] = t.colors[t.V3][2];
                temp2.colors[1][0] = t.colors[t.V2][0];
                temp2.colors[1][1] = t.colors[t.V2][1];
                temp2.colors[1][2] = t.colors[t.V2][2];
                
                //find the split point
                //if vertical side then its easy
                if(t.X[t.V1] == t.X[t.V3]){
                    temp1.X[2] = t.X[t.V1];

                    temp2.X[2] = t.X[t.V3];
                }
                //if not find with math
                else{
                    double dy = t.Y[t.V3] - t.Y[t.V1];
                    double dx = t.X[t.V3] - t.X[t.V1];
                    double m = dy/dx;
                    double b = t.Y[t.V1] - (t.X[t.V1] * m);
                    double splitX = (b + t.Y[t.V2]*-1)/m*-1;

                    temp1.X[2] = splitX;
                    temp2.X[2] = splitX;
                }
                //interpolate the color at X2
                temp1.colors[2][0] = interpolate(t.Y[t.V1], t.Y[t.V3], temp1.Y[2], t.colors[t.V1][0], t.colors[t.V3][0]);
                temp1.colors[2][1] = interpolate(t.Y[t.V1], t.Y[t.V3], temp1.Y[2], t.colors[t.V1][1], t.colors[t.V3][1]);
                temp1.colors[2][2] = interpolate(t.Y[t.V1], t.Y[t.V3], temp1.Y[2], t.colors[t.V1][2], t.colors[t.V3][2]); 
                temp2.colors[2][0] = temp1.colors[2][0];
                temp2.colors[2][1] = temp1.colors[2][1];
                temp2.colors[2][2] = temp1.colors[2][2];
                //interpolate Z[2]
                temp1.Z[2] = interpolate(t.Y[t.V1], t.Y[t.V3], temp1.Y[2], t.Z[t.V1], t.Z[t.V3]);
                temp2.Z[2] = temp1.Z[2]; //should be the same
                //interpolate the shading 
                temp1.shading[2] = interpolate(t.Y[t.V1], t.Y[t.V3], temp1.Y[2], t.shading[t.V1], t.shading[t.V3]);
                temp2.shading[2] = temp1.shading[2];
                
                //draw the triangles
                T_Type type1 = temp1.sortvertices();
                rasterize_triangle(temp1, screen, type1);
                T_Type type2 = temp2.sortvertices();
                rasterize_triangle(temp2, screen, type2);
                break;
            }
            default:{
                rasterize_triangle(t, screen, type);
                break;
            }
        }
}

void TransformTriangle(std::vector<Triangle> *triangles, Screen *screen, Camera *c){ 
    Matrix CT = c->CameraTransform();
    Matrix VT = c->ViewTransform();
    Matrix DT = c->DeviceTransform(screen->height, screen->width);
    Matrix M;
    M = M.ComposeMatrices(CT, VT);
    M = M.ComposeMatrices(M, DT);
    
    for(Triangle t: *triangles){
            //shade before transforming
            shade_triangle(&t, c);
            //transform
            double v1[4] = {t.X[0], t.Y[0], t.Z[0], 1};
            double v2[4] = {t.X[1], t.Y[1], t.Z[1], 1};
            double v3[4] = {t.X[2], t.Y[2], t.Z[2], 1};
            double v1_out[4]; 
            double v2_out[4]; 
            double v3_out[4];

            M.TransformPoint(v1, v1_out);
            M.TransformPoint(v2, v2_out);
            M.TransformPoint(v3, v3_out);
            
            if(v1_out[3] != 0){
                v1_out[0] = v1_out[0] / v1_out[3];
                v1_out[1] = v1_out[1] / v1_out[3];
                v1_out[2] = v1_out[2] / v1_out[3];
            }
            if(v2_out[3] != 0){
                v2_out[0] = v2_out[0] / v2_out[3];
                v2_out[1] = v2_out[1] / v2_out[3];
                v2_out[2] = v2_out[2] / v2_out[3];
            }
            if(v3_out[3] != 0){
                v3_out[0] = v3_out[0] / v3_out[3];
                v3_out[1] = v3_out[1] / v3_out[3];
                v3_out[2] = v3_out[2] / v3_out[3];
            }
            
            t.X[0] = v1_out[0]; t.Y[0] = v1_out[1]; t.Z[0] = v1_out[2];
            t.X[1] = v2_out[0]; t.Y[1] = v2_out[1]; t.Z[1] = v2_out[2];
            t.X[2] = v3_out[0]; t.Y[2] = v3_out[1]; t.Z[2] = v3_out[2];

            rasterize_general(t, screen);
            
            //set the triangle back
            t.X[0] = v1[0]; t.Y[0] = v1[1]; t.Z[0] = v1[2];
            t.X[1] = v2[0]; t.Y[1] = v2[1]; t.Z[1] = v2[2];
            t.X[2] = v3[0]; t.Y[2] = v3[1]; t.Z[2] = v3[2];
    }
}

int main()
{
   vtkImageData *image = NewImage(1000, 1000);
   unsigned char *buffer = 
     (unsigned char *) image->GetScalarPointer(0,0,0);
   double depthBuff[1000*1000];
   

   std::vector<Triangle> triangles = GetTriangles();
   
   Screen screen;
   screen.buffer = buffer;
   screen.depthBuff = depthBuff;
   screen.width = 1000;
   screen.height = 1000;

   

   for(int i = 0; i < 1; i++){
        screen.ClearScreen();
        Camera c = GetCamera(i, 1000);
        TransformTriangle(&triangles, &screen, &c);
        char name[20];
        sprintf(name, "frames/frame%03d", i);
        WriteImage(image, name);

    }
}

