#include <fstream>
#include <cstdlib>
#include <cmath>
#include <vector>

using namespace std;

struct color {
    unsigned char r;
    unsigned char g;
    unsigned char b;

    color(unsigned char r, unsigned char g, unsigned char b) : r(r), g(g), b(b) {}
    color(): r(0), g(0), b(0) {}
};

struct vec3f {
    float x;
    float y;
    float z;

    vec3f(float x, float y, float z) : x(x), y(y), z(z) { }
    vec3f(): x(0), y(0), z(0) {}

    vec3f operator-(const vec3f& b) const {
        return vec3f(this->x - b.x, this->y - b.y, this->z - b.z);
    }

    float operator*(const vec3f& b) const {
        return this->x * b.x + this->y * b.y + this->z * b.z;
    }

    vec3f operator*(float a) const {
        return vec3f(a * this->x, a * this->y, a * this->z);
    }

    vec3f& operator/=(float a) {
        this->x /= a;
        this->y /= a;
        this->z /= a;
        return *this;
    }

};

struct ray {
    vec3f r0;
    vec3f s;

    ray(const vec3f &r0, const vec3f &s) : r0(r0), s(s) {}
    ray(float x0, float y0, float z0, float sx, float sy, float sz)
            : r0(x0, y0, z0), s(sx, sy, sz) {}
    ray() {}
};

struct sphere {
    vec3f c;
    float R;
    color cl;

    sphere(const vec3f &c, float R, color cl = color(255, 255, 255)) : c(c), R(R), cl(cl) { }
    sphere(float x, float y, float z, float R, color cl = color(255, 255, 255))
            : c(x, y, z), R(R), cl(cl) {}
    sphere() {}
};

const float INF = 1e10;

float sqr(float x) {
    return x * x;
}

float intersect(const ray& r, const sphere& sp) {
    vec3f dr = sp.c - r.r0;
    float D = sqr(r.s * dr) + sqr(sp.R) - dr * dr;
    if (D < 0) return 2 * INF;
    return r.s * dr - sqrt(D);
}

vec3f view(256, 256, -256);
vector<sphere> spheres;

color get_pixel(int x, int y) {
    vec3f s = vec3f(x, y, 0) - view;
    s /= sqrt(s * s);
    ray r(view, s);
    int min_id = 0;
    float min_dist = intersect(r, spheres[0]);
    for (int i = 1; i < spheres.size(); i++) {
        float dist = intersect(r, spheres[0]);
        if (dist < min_dist) {
            min_id = i;
            min_dist = dist;
        }
    }
    if (min_dist > INF) return color(0, 0, 0);
    return spheres[min_id].cl;
}

int main() {
    spheres.push_back(sphere(256, 256, 512, 256, color(255, 0, 0)));

    ofstream fout("E:\\C++ Projects\\RayTracer\\out.ppm", ios::out | ios::binary);
    fout << "P6 512 512 255 ";
    color pixel;
    for (int i = 0; i < 512; i++) {
        for (int j = 0; j < 512; j++) {
            pixel = get_pixel(i, j);
            fout << pixel.r << pixel.g << pixel.b;
        }
    }
    fout.close();
    system("convert \"E:\\C++ Projects\\RayTracer\\out.ppm\" \"E:\\C++ Projects\\RayTracer\\out.png\"");
    return 0;
}