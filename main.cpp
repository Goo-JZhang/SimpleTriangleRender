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

Vec3f camera_pos(0,0,5);
Vec3f camera_dir(0,0,-1);
Vec3f camera_up(0,1,0);
Vec3f camera_right = cross(camera_dir,camera_up);
Vec2f screen_size(1.5,1.5);
float camera_plane = 4;

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
    if (std::abs(u[2])>1e-3) return Vec3f(1.f-(u.x+u.y)/u.z, u.y/u.z, u.x/u.z);
    else return Vec3f(-1,1,1);
}


bool inTriangle2(Vec3f A, Vec3f B, Vec3f C, Vec3f P)
{
    Vec3f bc_screen  = barycentric(A, B, C, P);
    if (bc_screen.x<0 || bc_screen.y<0 || bc_screen.z<0) return 0;
    else return 1;
}

Vec3f dBCal(Vec3f *pts)//calculate how z changes when a point in plane ABC moves a unit in y direction
{
    Vec3f bc_S = barycentric(pts[0],pts[1],pts[2],Vec3f(0, 0, 0));
    Vec3f bc_E = barycentric(pts[0],pts[1],pts[2],Vec3f(0, 1, 0));
    return Vec3f((bc_E[0]-bc_S[0]), (bc_E[1]-bc_S[1]), (bc_E[2]-bc_S[2]));
}

void triangle3(Vec3f *pts, Vec2f *vts, Vec3f *vns, float *zbuffer, TGAImage &image, TGAImage &texture)
{
    //sort pts ascending according to x
    if(pts[0].x>pts[1].x) {std::swap(pts[0],pts[1]); std::swap(vts[0],vts[1]); std::swap(vns[0],vns[1]);}
    if(pts[0].x>pts[2].x) {std::swap(pts[0],pts[2]); std::swap(vts[0],vts[2]); std::swap(vns[0],vns[2]);}
    if(pts[1].x>pts[2].x) {std::swap(pts[1],pts[2]); std::swap(vts[1],vts[2]); std::swap(vns[1],vns[2]);}
    //A: pts[0];  B: pts[2]; C:pts[1];
    //we can figure out the points (int,int,float) near the edge of triangle ABC
    //for any x(int), we can corresponding y in the edge AB and the edge A-C-B
    if(int(pts[2].x)==int(pts[0].x)) return;
    int txwidth = texture.get_width();
    int txheight = texture.get_height();
    //std::cout<<txwidth<<","<<txheight<<std::endl;
    //y's moving direction: if C is above the line AB, is 1, otherwise -1
    int dy = (pts[1].y>pts[0].y+(pts[1].x-pts[0].x)*(pts[2].y-pts[0].y)/(pts[2].x-pts[0].x)?1:-1);
    int yS, yE;//yS is near the edge AB, while yE is near the edge A-C--B
    float z, dz;
    Vec3f dB = dBCal(pts);
    dz = dy*(dB[0]*pts[0].z+dB[1]*pts[1].z+dB[2]*pts[2].z);
    Vec2f Vt;
    Vec2f dVt = dy*(dB[0]*vts[0] + dB[1]*vts[1] + dB[2]*vts[2]);
    Vec3f Vn;
    Vec3f dVn = dy*(dB[0]*vns[0] + dB[1]*vns[1] + dB[2]*vns[2]);
    for(int x = int(ceil(pts[0].x)); x<=int(pts[1].x);x++)
    {   //from A to C, if int(pts[0].x)==int(pts[1].x), the loop will be skip, 
        //since A and C are in the same column, no points in the plane (x, . , .) belongs to the triangle
        //so we need not to worry about the case of C.x-A.x == 0.f
        //calculate the range of y
        if(x<0) continue;
        if(x>=width) return;
        yS = pts[0].y + ((x-pts[0].x)/(pts[2].x-pts[0].x))*(pts[2].y-pts[0].y);
        Vec3f bc_screenStart = barycentric(pts[0],pts[1],pts[2],Vec3f(x,yS,0));
        yE = pts[0].y + ((x-pts[0].x)/(pts[1].x-pts[0].x))*(pts[1].y-pts[0].y);
        z = bc_screenStart[0]*pts[0].z + bc_screenStart[1]*pts[1].z + bc_screenStart[2]*pts[2].z;
        Vt = bc_screenStart[0]*vts[0] + bc_screenStart[1]*vts[1] + bc_screenStart[2]*vts[2];
        Vn = bc_screenStart[0]*vns[0] + bc_screenStart[1]*vns[1] + bc_screenStart[2]*vns[2];
        if(bc_screenStart.x<0||bc_screenStart.y<0||bc_screenStart.z<0)
        {//if (x,yAB) not in triangle, move (x, yAB,z) to (x, yAB+dy, z+dz)
            yS += dy;
            z  += dz;
            Vt += dVt;
            Vn += dVn;
        }
        if(!inTriangle2(pts[0],pts[1],pts[2],Vec3f(x,yE,0))) yE-=dy;
        for(int y = yS; dy*y<=dy*yE; y+=dy)
        {
            if(y<0){if(dy>0) continue; else break;}
            if(y>=height){if(dy<0){z+=dz; Vn+=dVn; Vt+=dVt ; continue;} else break;}
            if(zbuffer[x+y*width]>z&&z>camera_plane)
            {
                zbuffer[x+y*width]=z;
                float intensity = -(Vn*light_dir);
                //std::cout<<intensity<<std::endl;
                if(intensity>0) image.set(x,y,intensity*texture.get(Vt.x*txwidth,Vt.y*txheight));
            }
            z +=dz;
            Vn+=dVn;
            Vt+=dVt;
        }
    }
    for(int x = int(pts[1].x)+1; x<=int(pts[2].x);x++)
    {   //calculate from C to B.
        //std::cout<<x<<std::endl;
        if(x<0) continue;
        if(x>=width) return;
        yS = pts[0].y + ((x-pts[0].x)/(pts[2].x-pts[0].x))*(pts[2].y-pts[0].y);
        Vec3f bc_screenStart = barycentric(pts[0],pts[1],pts[2],Vec3f(x,yS,0));
        yE = pts[1].y + ((x-pts[1].x)/(pts[2].x-pts[1].x))*(pts[2].y-pts[1].y);
        z = bc_screenStart[0]*pts[0].z + bc_screenStart[1]*pts[1].z + bc_screenStart[2]*pts[2].z;
        Vt = bc_screenStart[0]*vts[0] + bc_screenStart[1]*vts[1] + bc_screenStart[2]*vts[2];
        Vn = bc_screenStart[0]*vns[0] + bc_screenStart[1]*vns[1] + bc_screenStart[2]*vns[2];
        if(bc_screenStart.x<0||bc_screenStart.y<0||bc_screenStart.z<0)
        {//if (x,yAB) not in triangle, move (x, yAB,z) to (x, yAB+dy, z+dz)
            yS += dy;
            z  += dz;
            Vt += dVt;
            Vn += dVn;
        }
        if(!inTriangle2(pts[0],pts[1],pts[2],Vec3f(x,yE,0))) yE-=dy;
        for(int y = yS; dy*y<=dy*yE; y+=dy)
        {
            if(y<0){if(dy>0) continue; else break;}
            if(y>=height){if(dy<0){z+=dz; Vn+=dVn; Vt+=dVt ; continue;} else break;}
            if(zbuffer[x+y*width]>z&&z>camera_plane)
            {
                zbuffer[x+y*width]=z;
                float intensity = -(Vn*light_dir);
                //std::cout<<intensity<<std::endl;
                if(intensity>0) image.set(x,y,intensity*texture.get(Vt.x*txwidth,Vt.y*txheight));
            }
            z +=dz;
            Vn+=dVn;
            Vt+=dVt;
        }
    }
}

Vec3f World2Screen(Vec3f worldc)
{
    Vec3f temp_v = worldc-camera_pos;
    float z = temp_v*camera_dir;
    float r;
    if(z<camera_plane) r = 1;
    else r = camera_plane/z;
    Vec3f screen;
    screen.z = z;
    screen.y = ((temp_v*camera_up)*r+screen_size.y/2.)*height/screen_size.y;
    screen.x = ((temp_v*camera_right)*r+screen_size.x/2.)*width/screen_size.x;
    return screen;
}

Vec3f world2screen(Vec3f worldc)
{
    return Vec3f((worldc.x+1)*width/2., (worldc.y+1)*height/2., worldc.z);
}

int main(int argc, char** argv) 
{
    //std::cout<<0.2*red<<std::endl;
    clock_t start = clock();
    if (2==argc) {
        model = new Model(argv[1]);
    } else {
        model = new Model("obj/african_head.obj");
    }
    float *zbuffer = new float[width*height];
    for (int i=width*height; i--; zbuffer[i] = std::numeric_limits<float>::max());
    TGAImage image(width, height, TGAImage::RGB);
    model->get_texture("obj/african_head_diffuse.tga");
    TGAImage texture = model->texture();
    for (int i=0; i<model->nfaces(); i++) {
        //std::cout<<model->nfaces()<<","<<i<<std::endl;
        //srand(time(0));
        //int i = 1179;
        //std::cout<<i<<std::endl;
        std::vector<int> face = model->face(i);
        Vec3f worlds[3];
        Vec3f pts[3];
        Vec2f txs[3];
        Vec3f nvs[3];
        for (int j=0; j<3; j++)//read verts
        {
            worlds[j] = model->vert(face[3*j]);
            txs[j] = model->txvert(face[3*j+1]);
            nvs[j] = model->nvert(face[3*j+2]);
            pts[j] = World2Screen(worlds[j]);
            //std::cout<<pts[j]<<","<<txs[j]<<","<<nvs[j]<<std::endl;
            //pts[j] = Vec3f(int((worlds[j].x+1.)*width/2. + 0.5 ), int((worlds[j].y+1.)*height/2. + 0.5), worlds[j].z);
        }
        //if(i==24){
            //std::cout<<worlds[0]<<" , "<<worlds[1]<<" , "<<worlds[2]<<std::endl;
            //std::cout<<pts[0]<<" , "<<pts[1]<<" , "<<pts[2]<<std::endl;
        //} 
        triangle3(pts, txs, nvs , zbuffer, image, texture);
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
    TGAImage texture;
    texture.read_tga_file("obj/african_head_diffuse.tga");
    Vec3f pts[3], vns[3];
    Vec2f vts[3];
    pts[0] = Vec3f(400.284, 783.855, 0.167019);
    pts[1] = Vec3f(400.284, 752.482, 0.320973);
    pts[2]=  Vec3f(427.583, 750.817, 0.315395);
    //std::cout<<pts[0]<<","<<pts[1]<<","<<pts[2]<<std::endl;
    triangle3(pts, vts, vns,zbuffer,image, texture);
    //std::cout<<image.get(400,400)<<std::endl;
    image.flip_vertically(); // i want to have the origin at the left bottom corner of the image
    image.write_tga_file("test.tga");
    //Vec3f k = barycentric(pts[0],pts[1],pts[2],Vec3f(401, 390.364, 0.604801));
    //std::cout<<k<<std::endl;
    return 0;
}
*/