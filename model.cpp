#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <vector>
#include "model.h"

Model::Model(const char *filename) : verts_(), txverts_(), nverts_(), faces_() {
    std::ifstream in;
    in.open (filename, std::ifstream::in);
    if (in.fail()) return;
    std::string line;
    while (!in.eof()) {
        std::getline(in, line);
        std::istringstream iss(line.c_str());
        char trash;
        if (!line.compare(0, 2, "v ")) {
            iss >> trash;
            Vec3f v;
            for (int i=0;i<3;i++) iss >> v[i];
            verts_.push_back(v);
        } else if (!line.compare(0, 2, "f ")) {
            std::vector<int> f;
            int a,b,c;
            iss >> trash;
            while (iss >> a >> trash >> b >> trash >> c) {
                a--; // in wavefront obj all indices start at 1, not zero
                f.push_back(a);
                b--;
                f.push_back(b);
                c--;
                f.push_back(c);
            }
            faces_.push_back(f);
        } else if(!line.compare(0,3, "vt ")) {
            iss>>trash>>trash;
            Vec2f vt;
            for (int i=0; i<2;i++) iss>>vt[i];
            txverts_.push_back(vt);
        } else if(!line.compare(0,3, "vn ")){
            iss>>trash>>trash;
            Vec3f vn;
            for(int i=0; i<3;i++) iss>>vn[i];
            nverts_.push_back(vn);
        }
    }
    std::cerr << "# v# " << verts_.size() << " f# "  << faces_.size() << std::endl;
}

Model::~Model() {
}

int Model::nverts() {
    return (int)verts_.size();
}

int Model::nfaces() {
    return (int)faces_.size();
}

std::vector<int> Model::face(int idx) {
    return faces_[idx];
}

Vec3f Model::vert(int i) {
    return verts_[i];
}

Vec2f Model::txvert(int i){
    return txverts_[i];
}

Vec3f Model::nvert(int i){
    return nverts_[i];
}