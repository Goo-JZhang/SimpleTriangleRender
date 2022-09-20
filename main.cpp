#include <vector>
#include <cmath>
#include <iostream>
#include <limits>
#include <ctime>
#include "tgaimage.h"
#include "model.h"
#include "geometry.h"

const TGAColor background = TGAColor(0,0,0,0);
const TGAColor black = TGAColor(0  , 0  , 0  , 255);
const TGAColor white = TGAColor(255, 255, 255, 255);
const TGAColor red   = TGAColor(255, 0  , 0  , 255);
const TGAColor green = TGAColor(0  , 255, 0  , 255);
const TGAColor blue  = TGAColor(0  , 0  , 255, 255);
Model *model = NULL;
const int width  = 800;
const int height = 800;
Vec3f light_dir(0,0,-1);

void line(Vec2i p0, Vec2i p1, TGAImage &image, TGAColor color) 
{
    bool steep=false;
    if (std::abs(p0.x-p1.x)<std::abs(p0.y-p1.y))
    {
        steep = true;
        std::swap(p0.x,p0.y);
        std::swap(p1.x,p1.y);
    }
    if (p0.x>p1.x) std::swap(p0,p1);
    int dy = std::abs(p1.y-p0.y) , dx = std::abs(p1.x-p0.x), error = 0;
    int y = p0.y;
    if (steep)
    {
        for(int x = p0.x; x<=p1.x;x++)
        {
            image.set(y,x,color);
            if(error+dy>=dx)
            {
                y+=(p1.y>p0.y?1:-1);
                //image.set(y,x,color);
                error-=dx;
            }
            error +=dy;
        }
    }
    else
    {
        for(int x = p0.x; x<=p1.x;x++)
        {
            image.set(x,y,color);
            if(error+dy>=dx)
            {
                y+=(p1.y>p0.y?1:-1);
                //image.set(x,y,color);
                error-=dx;
            }
            error +=dy;
        }
    }
}

Vec3f barycentric(Vec3f A, Vec3f B, Vec3f C, Vec3f P)
{
    //P=A+u*AB+v*AC, u,v both in [0,1]
    //u*AB+v*AC+PA=0
    Vec3f s[2];
    for(int i=2;i--;)
    {
        s[i][0]=C[i]-A[i];
        s[i][1]=B[i]-A[i];
        s[i][2]=A[i]-P[i];
    }
    Vec3f u = cross(s[0],s[1]);
    if (std::abs(u[2])>1e-6) return Vec3f(1.f-(u.x+u.y)/u.z, u.y/u.z, u.x/u.z);
    else return Vec3f(-1,1,1);
}


bool inTriangle2(Vec3f A, Vec3f B, Vec3f C, Vec3f P)
{
    Vec3f bc_screen  = barycentric(A, B, C, P);
    if (bc_screen.x<0 || bc_screen.y<0 || bc_screen.z<0) return 0;
    else return 1;
}

float dzCal(Vec3f *pts)//calculate how z changes when a point in plane ABC moves a unit in y direction
{
    Vec3f bc_S = barycentric(pts[0],pts[1],pts[2],Vec3f(0, 0, 0));
    Vec3f bc_E = barycentric(pts[0],pts[1],pts[2],Vec3f(0, 1, 0));
    return (bc_E[0]-bc_S[0])*pts[0].z + (bc_E[1]-bc_S[1])*pts[1].z + (bc_E[2]-bc_S[2])*pts[2].z;
}

void triangle3(Vec3f *pts, float *zbuffer, TGAImage &image, TGAColor color)
{
    //sort pts ascending according to x
    if(pts[0].x>pts[1].x) std::swap(pts[0],pts[1]);
    if(pts[0].x>pts[2].x) std::swap(pts[0],pts[2]);
    if(pts[1].x>pts[2].x) std::swap(pts[1],pts[2]);
    //A: pts[0];  B: pts[2]; C:pts[1];
    //we can figure out the points (int,int,float) near the edge of triangle ABC
    //for any x(int), we can corresponding y in the edge AB and the edge A-C-B
    if(int(pts[2].x)==int(pts[0].x)) return;
    int width = image.get_width();
    //y's moving direction: if C is above the line AB, is 1, otherwise -1
    int dy = (pts[1].y>pts[0].y+(pts[1].x-pts[0].x)*(pts[2].y-pts[0].y)/(pts[2].x-pts[0].x)?1:-1);
    int yS, yE;//yS is near the edge AB, while yE is near the edge A-C--B
    float z, dz = dy*dzCal(pts);// if 
    for(int x = int(floor(pts[0].x)); x<=int(pts[1].x);x++)
    {   //from A to C, if int(pts[0].x)==int(pts[1].x), the loop will be skip, 
        //since A and C are in the same column, no points in the plane (x, . , .) belongs to the triangle
        //so we need not to worry about the case of C.x-A.x == 0.f
        //calculate the range of y
        //std::cout<<x<<std::endl;
        yS = pts[0].y + ((x-pts[0].x)/(pts[2].x-pts[0].x))*(pts[2].y-pts[0].y);
        Vec3f bc_screenStart = barycentric(pts[0],pts[1],pts[2],Vec3f(x,yS,0));
        yE = pts[0].y + ((x-pts[0].x)/(pts[1].x-pts[0].x))*(pts[1].y-pts[0].y);
        z = bc_screenStart[0]*pts[0].z + bc_screenStart[1]*pts[1].z + bc_screenStart[2]*pts[2].z;
        if(bc_screenStart.x<0||bc_screenStart.y<0||bc_screenStart.z<0)
        {//if (x,yAB) not in triangle, move (x, yAB,z) to (x, yAB+dy, z+dz)
            yS+=dy;
            z+=dz;
        }
        if(inTriangle2(pts[0],pts[1],pts[2],Vec3f(x,yE,0))) yE-=dy;
        for(int y = yS; dy*y<=dy*yE; y+=dy)
        {
            //std::cout<<x<<","<<y<<std::endl;
            if(zbuffer[x+y*width]<z)
            {
                zbuffer[x+y*width]=z;
                image.set(x,y,color);
            }
        }
    }
    for(int x = int(pts[1].x)+1; x<=int(pts[2].x);x++)
    {   //calculate from C to B.
        //std::cout<<x<<std::endl;
        yS = pts[0].y + ((x-pts[0].x)/(pts[2].x-pts[0].x))*(pts[2].y-pts[0].y);
        Vec3f bc_screenStart = barycentric(pts[0],pts[1],pts[2],Vec3f(x,yS,0));
        yE = pts[1].y + ((x-pts[1].x)/(pts[2].x-pts[1].x))*(pts[2].y-pts[1].y);
        z = bc_screenStart[0]*pts[0].z + bc_screenStart[1]*pts[1].z + bc_screenStart[2]*pts[2].z;
        if(bc_screenStart.x<0||bc_screenStart.y<0||bc_screenStart.z<0)
        {//if (x,yAB) not in triangle, move (x, yAB,z) to (x, yAB+dy, z+dz)
            yS+=dy;
            z+=dz;
        }
        if(inTriangle2(pts[0],pts[1],pts[2],Vec3f(x,yE,0))) yE-=dy;
        for(int y = yS; dy*y<=dy*yE; y+=dy)
        {
            //std::cout<<x<<","<<y<<std::endl;
            if(zbuffer[x+y*width]<z)
            {
                zbuffer[x+y*width]=z;
                image.set(x,y,color);
            }
        }
    }
}

void triangle3b(Vec3f *pts, float *zbuffer, TGAImage &image, TGAColor color) 
{
    Vec2f bboxmin( std::numeric_limits<float>::max(),  std::numeric_limits<float>::max());
    Vec2f bboxmax(-std::numeric_limits<float>::max(), -std::numeric_limits<float>::max());
    Vec2f clamp(image.get_width()-1, image.get_height()-1);
    for (int i=0; i<3; i++) {
        for (int j=0; j<2; j++) {
            bboxmin[j] = std::max(0.f,      std::min(bboxmin[j], pts[i][j]));
            bboxmax[j] = std::min(clamp[j], std::max(bboxmax[j], pts[i][j]));
        }
    }
    Vec3f P;
    for (P.x=bboxmin.x; P.x<=bboxmax.x; P.x++) {
        for (P.y=bboxmin.y; P.y<=bboxmax.y; P.y++) {
            Vec3f bc_screen  = barycentric(pts[0], pts[1], pts[2], P);
            if (bc_screen.x<0 || bc_screen.y<0 || bc_screen.z<0) continue;
            P.z = 0;
            for (int i=0; i<3; i++) P.z += pts[i][2]*bc_screen[i];
            if (zbuffer[int(P.x+P.y*width)]<P.z) {
                zbuffer[int(P.x+P.y*width)] = P.z;
                image.set(P.x, P.y, color);
            }
        }
    }
}

int main(int argc, char** argv) 
{
    clock_t start = clock();
    if (2==argc) {
        model = new Model(argv[1]);
    } else {
        model = new Model("obj/african_head.obj");
    }
    float *zbuffer = new float[width*height];
    for (int i=width*height; i--; zbuffer[i] = -std::numeric_limits<float>::max());
    TGAImage image(width, height, TGAImage::RGB);
    for (int i=0; i<model->nfaces(); i++) {
        //std::cout<<model->nfaces()<<","<<i<<std::endl;
        //srand(int(time(0)));
        //int i = rand()%model->nfaces();
        //int i =1827;
        std::vector<int> face = model->face(i);
        Vec3f worlds[3];
        Vec3f pts[3];
        for (int j=0; j<3; j++)
        {
            worlds[j] = model->vert(face[j]);
            pts[j] = Vec3f((worlds[j].x+1.)*width/2. , (worlds[j].y+1.)*height/2. , worlds[j].z);
            //pts[j] = Vec3f(int((worlds[j].x+1.)*width/2. + 0.5 ), int((worlds[j].y+1.)*height/2. + 0.5), worlds[j].z);
        }
        //if(i==1015) std::cout<<pts[0]<<","<<pts[1]<<","<<pts[2]<<std::endl;
        Vec3f n = cross(worlds[2]-worlds[0],worlds[1]-worlds[0]).normalize();
        float intensity = n*light_dir;
        if (intensity>0) 
        {
            //std::cout<<i<<","<<pts[0]<<","<<pts[1]<<","<<pts[2]<<std::endl;
            triangle3(pts, zbuffer, image, TGAColor(intensity*255, intensity*255, intensity*255, 255));
        }
    }
    image.flip_vertically(); // i want to have the origin at the left bottom corner of the image
    image.write_tga_file("african_head.tga");
    delete model;
    std::cout<<(clock()-start)/1000.<<std::endl;
    return 0;
}

/*
int main()
{
    float *zbuffer = new float[width*height];
    for (int i=width*height; i--; zbuffer[i] = -std::numeric_limits<float>::max());
    TGAImage image(width, height, TGAImage::RGB);
    Vec3f pts[3];
    srand((int)time(0));
    //pts[0] = Vec3f(rand()%800, rand()%800, (rand()%100-50)/float(50));
    //pts[1] = Vec3f(rand()%800, rand()%800, (rand()%100-50)/float(50));
    //pts[2] = Vec3f(rand()%800, rand()%800, (rand()%100-50)/float(50));
    pts[0] = Vec3f(400.284, 752.482, 0.320973);
    pts[1] = Vec3f(400.324, 711.676, 0.437668);
    pts[2]=  Vec3f(425.998, 710.867, 0.432594);
    //std::cout<<pts[0]<<","<<pts[1]<<","<<pts[2]<<std::endl;
    triangle3(pts,zbuffer,image,red);
    //std::cout<<image.get(400,400)<<std::endl;
    image.flip_vertically(); // i want to have the origin at the left bottom corner of the image
    image.write_tga_file("test.tga");
    //Vec3f k = barycentric(pts[0],pts[1],pts[2],Vec3f(401, 390.364, 0.604801));
    //std::cout<<k<<std::endl;
    return 0;
}
*/
