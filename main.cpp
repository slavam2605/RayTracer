#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <ctime>
#include <memory>

using namespace std;

struct color {
    float r;
    float g;
    float b;

    color(float r, float g, float b) : r(r), g(g), b(b) { }

    color() : r(0), g(0), b(0) { }

    color operator*(const color& o) const {
        return color(r * o.r, g * o.g, b * o.b);
    }

    color& operator+=(const color& o) {
        r += o.r;
        g += o.g;
        b += o.b;
        return *this;
    }

    color& operator/=(float a) {
        r /= a;
        g /= a;
        b /= a;
        return *this;
    }

    color& normalize() {
        if (r < 0) r = 0.0f;
        if (r > 1) r = 1.0f;
        if (g < 0) g = 0.0f;
        if (g > 1) g = 1.0f;
        if (b < 0) b = 0.0f;
        if (b > 1) b = 1.0f;
        return *this;
    }

};

const float PI = 3.14159265358979323f;

float sqr(float x) {
    return x * x;
}

struct hsv;

inline hsv hsv2hsl(const hsv& cl);

inline hsv hsl2hsv(const hsv& cl);

struct hsv {
    float h;
    float s;
    float v;


    hsv(float h, float s, float v) : h(h), s(s), v(v) { }


    hsv() { }

    hsv mix(const hsv& o) {
        hsv hsl1 = hsv2hsl(*this);
        hsv hsl2 = hsv2hsl(o);
        float x1 = hsl1.s * cos(PI * hsl1.h / 180);
        float y1 = hsl1.s * sin(PI * hsl1.h / 180);
        float z1 = hsl1.v;
        float x2 = hsl2.s * cos(PI * hsl2.h / 180);
        float y2 = hsl2.s * sin(PI * hsl2.h / 180);
        float z2 = hsl2.v;
        x1 = (x1 + x2) / 2;
        y1 = (y1 + y2) / 2;
        z1 = (z1 + z2) / 2;
        return hsl2hsv(hsv((atan2f(y1, x1) + PI) * 180 / PI, sqrt(x1 * x1 + y1 * y1), z1));
    }
};

inline hsv hsv2hsl(const hsv& cl) {
    float hue2 = (2 - cl.s) * cl.v;
    return hsv(
            cl.h,
            cl.s * cl.v / (hue2 < 1 ? hue2 : 2 - hue2),
            hue2 / 2
    );
}

inline hsv hsl2hsv(const hsv& cl) {
    float sat = cl.s * (cl.v < 0.5 ? cl.v : 1 - cl.v);
    return hsv(
            cl.h,
            2 * sat / (cl.v + sat),
            cl.v + sat
    );
}

hsv rgb2hsv(color in) {
    hsv out;
    float min, max, delta;

    min = in.r < in.g ? in.r : in.g;
    min = min < in.b ? min : in.b;

    max = in.r > in.g ? in.r : in.g;
    max = max > in.b ? max : in.b;

    out.v = max;                                // v
    delta = max - min;
    if (delta < 0.00001f) {
        out.s = 0;
        out.h = 0; // undefined, maybe nan?
        return out;
    }
    if (max > 0.0f) { // NOTE: if Max is == 0, this divide would cause a crash
        out.s = (delta / max);                  // s
    } else {
        // if max is 0, then r = g = b = 0
        // s = 0, v is undefined
        out.s = 0.0f;
        out.h = NAN;                            // its now undefined
        return out;
    }
    if (in.r >= max)                           // > is bogus, just keeps compilor happy
        out.h = (in.g - in.b) / delta;        // between yellow & magenta
    else if (in.g >= max)
        out.h = 2.0f + (in.b - in.r) / delta;  // between cyan & yellow
    else
        out.h = 4.0f + (in.r - in.g) / delta;  // between magenta & cyan

    out.h *= 60.0f;                              // degrees

    if (out.h < 0.0f)
        out.h += 360.0f;

    return out;
}


color hsv2rgb(hsv in) {
    if (in.s <= 0.0) {
        return color(in.v, in.v, in.v);
    }
    float hh = in.h / 60.0f;
    if (hh >= 6) hh = 0.0f;
    long i = hh;
    float ff = hh - i;
    float p = in.v * (1.0f - in.s);
    float q = in.v * (1.0f - in.s * ff);
    float t = in.v * (1.0f - in.s * (1.0f - ff));

    switch (i) {
        case 0:
            return color(in.v, t, p);
        case 1:
            return color(q, in.v, p);
        case 2:
            return color(p, in.v, t);
        case 3:
            return color(p, q, in.v);
        case 4:
            return color(t, p, in.v);
        case 5:
        default:
            return color(in.v, p, q);
    }
}

struct vec3f {
    float x;
    float y;
    float z;

    vec3f(float x, float y, float z) : x(x), y(y), z(z) { }

    vec3f() : x(0), y(0), z(0) { }

    vec3f operator+(const vec3f& b) const {
        return vec3f(this->x + b.x, this->y + b.y, this->z + b.z);
    }

    vec3f operator-(const vec3f& b) const {
        return vec3f(this->x - b.x, this->y - b.y, this->z - b.z);
    }

    float operator*(const vec3f& b) const {
        return this->x * b.x + this->y * b.y + this->z * b.z;
    }

    vec3f operator*(float a) const {
        return vec3f(a * this->x, a * this->y, a * this->z);
    }

    vec3f operator/(float a) const {
        return vec3f(this->x / a, this->y / a, this->z / a);
    }

    vec3f& operator/=(float a) {
        this->x /= a;
        this->y /= a;
        this->z /= a;
        return *this;
    }

    vec3f& operator*=(float a) {
        this->x *= a;
        this->y *= a;
        this->z *= a;
        return *this;
    }

};

struct ray {
    vec3f r0;
    vec3f s;

    ray(const vec3f& r0, const vec3f& s) : r0(r0), s(s) { }

    ray(float x0, float y0, float z0, float sx, float sy, float sz)
            : r0(x0, y0, z0), s(sx, sy, sz) { }

    ray() { }
};

const float INF = 1e10;

struct object {
    virtual float intersect(const ray& r) = 0;

    virtual color get_color(const ray& r, const vec3f& p, int depth) = 0;
};

vec3f view(256 / 2, 256 / 2, -256);
vector<unique_ptr<object>> scene;

vec3f rand_dir() {
    while (true) {
        int x = 2 * rand() - RAND_MAX;
        int y = 2 * rand() - RAND_MAX;
        int z = 2 * rand() - RAND_MAX;
        int len = x * x + y * y + z * z;
        if (len <= RAND_MAX * RAND_MAX) {
            float sq_len = sqrtf(len);
            return vec3f(x / sq_len, y / sq_len, z / sq_len);
        }
    }
}

const int DEPTH_LIMIT = 1;

color get_color(const ray& r, int depth);

const int NN = 500;

struct sphere : object {
    vec3f c;
    float R;
    color cl;
    bool is_light;
    bool is_mirror;

    sphere(const vec3f& c, float R, color cl = color(1, 1, 1), bool is_light = false, bool is_mirror = false)
            : c(c), R(R), cl(cl), is_light(is_light), is_mirror(is_mirror) { }

    sphere(float x, float y, float z, float R, color cl = color(1, 1, 1), bool is_light = false, bool is_mirror = false)
            : c(x, y, z), R(R), cl(cl), is_light(is_light), is_mirror(is_mirror) { }

    sphere() { }

    float intersect(const ray& r) {
        vec3f dr = c - r.r0;
        float D = sqr(r.s * dr) + sqr(R) - dr * dr;
        if (D < 1e-5) return 2 * INF;
        float result = r.s * dr - sqrt(D);
        if (result < 0) return 2 * INF; else return result;
    }

    color get_color(const ray& r, const vec3f& p, int depth) {
        if (is_light) return cl;
        if (depth > DEPTH_LIMIT) return color(0, 0, 0);
        if (is_mirror) {
            vec3f norm = (p - c) / R;
            vec3f s = r.s - norm * 2 * (r.s * norm);
            return ::get_color(ray(p, s), depth);
        }
        color light_c;
        vec3f norm = (p - c) / R;
        int n = NN;
        for (int i = 0; i < n; i++) {
            vec3f dir = rand_dir();
            if (dir * norm < 0) dir *= -1;
            light_c += ::get_color(ray(p, dir), depth + 1);
        }
        light_c /= n;
        return cl * light_c;
    }
};

struct plane : object {
    vec3f norm;
    float D;
    color cl;
    bool is_mirror;

    plane(const vec3f& norm, float D, const color& cl = color(1, 1, 1), bool is_mirror = false)
            : norm(norm), D(D), cl(cl), is_mirror(is_mirror) { }

    plane(const vec3f& norm, const vec3f& p, const color& cl = color(1, 1, 1), bool is_mirror = false)
            : norm(norm), D(norm * (-1) * p), cl(cl), is_mirror(is_mirror) { }

    float intersect(const ray& r) {
        float result = -(D + r.r0 * norm) / (r.s * norm);
        if (result < 1e-5) return 2 * INF;
        return result;
    }

    color get_color(const ray& r, const vec3f& p, int depth) {
        if (depth > DEPTH_LIMIT) return color(0, 0, 0);
        if (is_mirror) {
            vec3f s = r.s - norm * 2 * (r.s * norm);
            return ::get_color(ray(p, s), depth);
        }
        color light_c;
        int n = NN;
        for (int i = 0; i < n; i++) {
            vec3f dir = rand_dir();
            if (dir * norm < 0) dir *= -1;
            light_c += ::get_color(ray(p, dir), depth + 1);
        }
        light_c /= n;
        return cl * light_c;
    }
};

ofstream& print(ofstream& fout, const color& pixel);

color get_color(float x, float y) {
    vec3f s = vec3f(x, y, 0) - view;
    s /= sqrt(s * s);
    ray r(view, s);
    return get_color(r, 1).normalize();
}

color get_color(const ray& r, int depth) {
    int min_id = 0;
    float min_dist = scene[0]->intersect(r);
    for (int i = 1; i < scene.size(); i++) {
        float dist = scene[i]->intersect(r);
        if (dist < min_dist) {
            min_id = i;
            min_dist = dist;
        }
    }
    if (min_dist > INF) return color(0, 0, 0);
    return scene[min_id]->get_color(r, r.s * min_dist + r.r0, depth);
}

color mix2(color c1, color c2) {
    return hsv2rgb(rgb2hsv(c1).mix(rgb2hsv(c2)));
}

inline color get_pixel(int x, int y) {
    return get_color(x, y);
//    return mix2(
//            mix2(get_color(x - 0.25f, y - 0.25f),
//                 get_color(x - 0.25f, y + 0.25f)),
//            mix2(get_color(x + 0.25f, y - 0.25f),
//                 get_color(x + 0.25f, y + 0.25f)));
}

int main() {
    int w = 256;
    int h = 256;

//    scene.push_back(unique_ptr<sphere>(new sphere(w / 2 - w / 4, h / 2, w / 2, w / 2, color(1, 0, 0))));
//    scene.push_back(unique_ptr<sphere>(new sphere(w / 2 + w / 6, h / 2, w / 2, w / 2, color(0, 1, 0))));
//    scene.push_back(unique_ptr<sphere>(new sphere(w, 0, 0, w / 4, color(10, 10, 10), true)));

    int MARGIN = 4 * w / 3;

    scene.push_back(unique_ptr<sphere>(new sphere(0, 2 * h / 3, w / 2, w / 3, color(1, 0, 0), false, true)));
    scene.push_back(unique_ptr<sphere>(new sphere(w / 2 + w / 4 - w / 6, h / 4 + h / 6, w / 2, w / 6, color(0, 1, 0))));
    scene.push_back(unique_ptr<sphere>(new sphere(0, 0, w, w / 5, color(20, 20, 20), true)));
    scene.push_back(unique_ptr<sphere>(new sphere(w, 0, w, w / 5, color(20, 20, 20), true)));
    scene.push_back(unique_ptr<sphere>(new sphere(0, 0, 0, w / 5, color(20, 20, 20), true)));
    scene.push_back(unique_ptr<sphere>(new sphere(w, 0, 0, w / 5, color(20, 20, 20), true)));
    scene.push_back(unique_ptr<plane>(new plane(vec3f(-1, 0, 0), vec3f(w + MARGIN, 0, 0), color(1, 1, 0))));
    scene.push_back(unique_ptr<plane>(new plane(vec3f(1, 0, 0), vec3f(0 - MARGIN, 0, 0), color(1, 1, 0))));
    scene.push_back(unique_ptr<plane>(new plane(vec3f(0, -1, 0), vec3f(0, h + MARGIN, 0), color(1, 0, 1))));
    scene.push_back(unique_ptr<plane>(new plane(vec3f(0, 1, 0), vec3f(0, 0 - MARGIN, 0), color(1, 0, 1))));
    scene.push_back(unique_ptr<plane>(new plane(vec3f(0, 0, -1), vec3f(0, 0, w + MARGIN), color(0, 1, 1))));
    scene.push_back(unique_ptr<plane>(new plane(vec3f(0, 0, 1), vec3f(0, 0, 0 - MARGIN), color(0, 1, 1))));

//    scene.push_back(unique_ptr<sphere>(new sphere(w / 2, h / 2, w / 2, w / 2, color(1, 0, 0))));
//    scene.push_back(unique_ptr<sphere>(new sphere(w, 0, 0, w / 4, color(10, 10, 10), true)));

    view = vec3f(w / 2, h / 2, -w / 2);
    int time = clock();
//    string path = "E:\\C++ Projects\\RayTracer\\";
    string path = "C:\\Users\\slava\\ClionProjects\\RayTracer\\RayTracer\\";
    ofstream fout(path + "out.ppm", ios::out | ios::binary);
    fout << "P6 " << w << " " << h << " 255 ";
    color pixel;
    unsigned char cl[w * h * 3];
    int pointer = 0;
    int last_flushed = 0;
    for (int j = 0; j < h; j++) {
        for (int i = 0; i < w; i++) {
            pixel = get_pixel(i, j);
            cl[pointer++] = pixel.r * 255;
            cl[pointer++] = pixel.g * 255;
            cl[pointer++] = pixel.b * 255;
        }
        if (100 * j / h > last_flushed) {
            cout << "[" << 100 * j / h << "%]";
            last_flushed = 100 * j / h;
            float progress = (float) j / h;
            cout << ", estimated time: "  << float(clock() - time) / CLOCKS_PER_SEC / progress * (1 - progress) << " seconds" << endl;
        }
    }
    fout.write((char*) cl, w * h * 3);
    fout.close();
    cout << float(clock() - time) / CLOCKS_PER_SEC << " seconds" << endl;
    system(("convert \"" + path + "out.ppm\" \"" + path + "out.png\"").c_str());
    return 0;
}