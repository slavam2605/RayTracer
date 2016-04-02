#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <ctime>
#include <memory>
#include <thread>
#include <iomanip>

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

    color operator*(float a) const {
        return color(r * a, g * a, b * a);
    }

    color operator/(float a) const {
        return color(r / a, g / a, b / a);
    }

    color operator+(const color& o) const {
        return color(r + o.r, g + o.g, b + o.b);
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

struct vec3f {
    float x;
    float y;
    float z;

    vec3f(float x, float y, float z) : x(x), y(y), z(z) { }

    vec3f() : x(0.0f), y(0.0f), z(0.0f) { }

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

    ray(const vec3f& r0, const vec3f& s) : r0(r0), s(s / sqrtf(s.x * s.x + s.y * s.y + s.z * s.z)) { }

    ray() { }
};

const float INF = 1e10;
const int NN = 500;

//unsigned int rand_state;
//#define RAND_MAX 4294967295
//#define mrand() (rand_state = 214013 * rand_state + 2531011)

vec3f rand_dir() {
    while (true) {
        float x = 2.0f * rand() / RAND_MAX - 1.0f;
        float y = 2.0f * rand() / RAND_MAX - 1.0f;
        float z = 2.0f * rand() / RAND_MAX - 1.0f;
        float len = x * x +  y * y + z * z;
        if (len > 1) continue;
        return vec3f(x, y, z);
    }
}

const int DEPTH_LIMIT = 10;
const int BUNCH_LIMIT = 1;

color get_color(const ray& r, float cn, int bunch_depth, int depth);

enum class surface_type {SIMPLE, LIGHT, MIRROR, TRANSPARENT, GLOSSY};

struct object {
    float on;
    bool is_light = false;
    bool is_mirror = false;
    bool is_transparent = false;
    bool is_glossy = false;

    object(surface_type type, float on = 1) {
        this->on = on;
        switch (type) {
            case surface_type::LIGHT:
                is_light = true;
                break;
            case surface_type::MIRROR:
                is_mirror = true;
                break;
            case surface_type::TRANSPARENT:
                is_transparent = true;
                break;
            case surface_type::GLOSSY:
                is_glossy = true;
                break;
            case surface_type::SIMPLE:
                break;
        }
    }

    virtual float intersect(const ray& r) const = 0;

    virtual vec3f get_norm(const ray& r, const vec3f& p) const = 0;

    virtual bool is_inside(const ray& r, const vec3f& p) const = 0;

    virtual color get_color(const vec3f& p) const = 0;

    color diffuse(const vec3f& p, float cn, int bunch_depth, int depth, const vec3f& norm) const {
        color light_c;
        int n = NN;
        for (int i = 0; i < n; i++) {
            vec3f dir = rand_dir();
            if (dir * norm < 0) dir *= -1;
            light_c += ::get_color(ray(p + norm, dir), cn, bunch_depth + 1, depth + 1);
        }
        light_c /= n;
        return light_c;
    }

    color reflect(const ray& r, const vec3f& p, float cn, int bunch_depth, int depth, const vec3f& norm) const {
        vec3f s = r.s - norm * 2 * (r.s * norm);
        return ::get_color(ray(p + norm, s), cn, bunch_depth, depth + 1);
    }

    color glossy(const ray& r, const vec3f& p, float cn, int bunch_depth, int depth, const vec3f& norm) const {
        float cosa = -(r.s * norm);
        float sina = sqrtf(fabs(1 - cosa * cosa));
        float sinb = cn / on * sina;
        float cosb = sqrtf(fabs(1 - sinb * sinb));
        float R = sqr((sina * cosb - cosa * sinb) / (sina * cosb + cosa * sinb));

        vec3f ss = r.s - norm * 2 * (r.s * norm);
        color c_refl;
        if (R > 1e-3) {
            c_refl = ::get_color(ray(p + norm, ss), cn, bunch_depth, depth + 1);
        }
        color c_difr = get_color(p) * diffuse(p, cn, bunch_depth, depth, norm);
        return c_refl * R + c_difr * (1 - R);
    }

    color full_refract(const ray& r, const vec3f& p, float cn, float n2, int bunch_depth, int depth, const vec3f& norm) const {
        vec3f s = r.s * cn;
        float D = (n2 * n2 - cn * cn) / (s * norm) / (s * norm) + 1;

        float cosa = -(r.s * norm);
        float sina = sqrtf(fabs(1 - cosa * cosa));
        float sinb = cn / n2 * sina;
        float cosb = sqrtf(fabs(1 - sinb * sinb));
        float R = sqr((sina * cosb - cosa * sinb) / (sina * cosb + cosa * sinb));

        color c_refl;
        if (D < 1e-5f || fabs(R) > 1e-3f) {
            c_refl = reflect(r, p, cn, bunch_depth, depth, norm);
        }
        if (D < 1e-5f) return c_refl;
        vec3f ns = s + norm * (s * norm) * (sqrtf(D) - 1);

        color c_refr = ::get_color(ray(p - norm, ns), n2, bunch_depth, depth + 1);
        return c_refl * R + c_refr * (1 - R);
    }

    color get_color(const ray& r, const vec3f& p, float cn, int bunch_depth, int depth) const {
        if (depth > DEPTH_LIMIT) return color(0, 0, 0);
        if (is_light) return get_color(p);
        vec3f norm = get_norm(r, p);
        if (is_mirror) {
            return reflect(r, p, cn, bunch_depth, depth, norm);
        }
        if (is_transparent) {
            float n2 = is_inside(r, p) ? 1 : on;
            return full_refract(r, p, cn, n2, bunch_depth, depth, norm);
        }
        if (bunch_depth > BUNCH_LIMIT) return color(0, 0, 0);
        if (is_glossy) {
            return glossy(r, p, cn, bunch_depth, depth, norm);
        }
        return get_color(p) * diffuse(p, cn, bunch_depth, depth, norm);
    }

};

vec3f view(256 / 2, 256 / 2, -256);
vector<unique_ptr<object>> scene;

struct sphere : object {
    vec3f c;
    float R;
    color cl;

    sphere(const vec3f& c, float R, color cl = color(1, 1, 1), surface_type type = surface_type::SIMPLE, float on = 1)
            : object(type, on), c(c), R(R), cl(cl) { }

    sphere(float x, float y, float z, float R, color cl = color(1, 1, 1), surface_type type = surface_type::SIMPLE, float on = 1)
            : sphere(vec3f(x, y, z), R, cl, type, on) { }

    float intersect(const ray& r) const {
        vec3f dr = c - r.r0;
        float D = sqr(r.s * dr) + sqr(R) - dr * dr;
        if (D < 1e-5) return 2 * INF;
        float result = r.s * dr - sqrtf(D);
        if (result < 0) {
            result = r.s * dr + sqrtf(D);
            if (result < 0) return 2 * INF; else return result;
        } else return result;
    }

    vec3f get_norm(const ray& r, const vec3f& p) const {
        vec3f norm = (p - c) / R;
        if (r.s * norm > 0) norm *= -1;
        return norm;
    }

    bool is_inside(const ray& r, const vec3f& p) const {
        return r.s * (p - c) > 0;
    }

    color get_color(const vec3f& p) const {
        return cl;
    }

};

struct plane : object {
    vec3f norm;
    float D;

    plane(const vec3f& norm, float D, surface_type type = surface_type::SIMPLE, float on = 1)
            : object(type, on), norm(norm), D(D) { }

    plane(const vec3f& norm, const vec3f& p, surface_type type = surface_type::SIMPLE, float on = 1)
            : plane(norm, -(norm * p), type, on) { }

    float intersect(const ray& r) const {
        float result = -(D + r.r0 * norm) / (r.s * norm);
        if (result < 1e-5) return 2 * INF;
        return result;
    }

    vec3f get_norm(const ray& r, const vec3f& p) const {
        return norm;
    }

    bool is_inside(const ray& r, const vec3f& p) const {
        return r.s * norm > 0;
    }

    color get_color(const vec3f& p) const {
        return ((int) floor(p.x / 100) + (int) floor(p.z / 100)) & 1 ? color(1, 0, 0) : color(1, 1, 1);
    }

};

const int w = 512;
const int h = 512;

color get_color(float x, float y) {
    int N = 0;
    float range = 10.0f;
    color cl;
    for (int i = -N; i <= N; i++) {
        for (int j = -N; j <= N; j++) {
            vec3f vv = view + vec3f(i * range / (N ? N : 1), j * range / (N ? N : 1), 0);
            vec3f s = vec3f(2 * x - w / 2, 2 * y - h / 2, 256) - vv;
            s /= sqrtf(s * s);
            ray r(vv, s);
            cl += get_color(r, 1, 1, 1).normalize();
        }
    }
    return cl / (2 * N + 1) / (2 * N + 1);
}

color get_color(const ray& r, float cn, int bunch_depth, int depth) {
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
    return scene[min_id]->get_color(r, r.s * min_dist + r.r0, cn, bunch_depth, depth);
}

color mix2(color c1, color c2) {
    return (c1 + c2) / 2;
}

inline color get_pixel(int x, int y) {
    return get_color(x, y);
//    return mix2(
//            mix2(get_color(x - 0.25f, y - 0.25f),
//                 get_color(x - 0.25f, y + 0.25f)),
//            mix2(get_color(x + 0.25f, y - 0.25f),
//                 get_color(x + 0.25f, y + 0.25f)));
}

unsigned char cl[w * h * 3];

int progress[8];

void print(int n) {
    for (int i = 0; i < n; i++) {
        if (progress[i] < 0)
            cout << "      ";
        else
            cout << " [" << (progress[i] < 10 ? "0" : "") << progress[i] << "%]";
    }
    cout << endl;
}

void perform(int k, int n) {
    color pixel;
    int offset;
    int start = w * k / n;
    int end = w * (k + 1) / n;
    for (int i = start; i < end; i++) {
        for (int j = 0; j < h; j++) {
            pixel = get_pixel(i, j);
            offset = 3 * (i + j * w);
            cl[offset] = pixel.r * 255;
            cl[offset + 1] = pixel.g * 255;
            cl[offset + 2] = pixel.b * 255;
        }
//        if (i % 10 == 0) {
            progress[k] = 100 * (i + 1 - start) / (end - start);
            print(n);
//        }
    }
    progress[k] = -1;
}

int main() {
//    SetPriorityClass(GetCurrentProcess(), IDLE_PRIORITY_CLASS);

    int MARGIN = 4 * w / 3;

    float f = 1;
//    scene.push_back(unique_ptr<sphere>(new sphere(w / 6, 2 * h / 3, w / 2, w / 3, color(1, 0, 0), false, true)));
//    scene.push_back(unique_ptr<sphere>(new sphere(5 * w / 6, 2 * h / 3, w / 2, w / 3, color(1, 0, 0), false, true)));
//    scene.push_back(unique_ptr<sphere>(new sphere(w / 2, h / 2, w / 2, w / 2, color(1, 1, 1), false, false, true, sqrtf(2))));          // 0
//    scene.push_back(unique_ptr<sphere>(new sphere(w / 2, h + f * MARGIN, w + f * MARGIN, w / 4, color(30, 30, 30), true)));             // 1
//    scene.push_back(unique_ptr<sphere>(new sphere(0 - f * MARGIN, -f * MARGIN, w + f * MARGIN, w / 5, color(30, 30, 30), true)));       // 2
//    scene.push_back(unique_ptr<sphere>(new sphere(w + f * MARGIN, -f * MARGIN, w + f * MARGIN, w / 5, color(30, 30, 30), true)));       // 3
//    scene.push_back(unique_ptr<sphere>(new sphere(0 - f * MARGIN, -f * MARGIN, 0 - f * MARGIN, w / 5, color(30, 30, 30), true)));       // 4
//    scene.push_back(unique_ptr<sphere>(new sphere(w + f * MARGIN, -f * MARGIN, 0 - f * MARGIN, w / 5, color(30, 30, 30), true)));       // 5



    scene.push_back(unique_ptr<sphere>(new sphere(w / 2, -h, w / 2, w / 5, color(100, 100, 100), surface_type::LIGHT)));
    scene.push_back(unique_ptr<sphere>(new sphere(w / 2, 3 * h / 4, w / 2, w / 4, color(1, 0, 0), surface_type::GLOSSY, 1.2f)));
    scene.push_back(unique_ptr<sphere>(new sphere(w, 3 * h / 4, w / 2, w / 4, color(1, 0, 0), surface_type::TRANSPARENT , 1.2f)));
    scene.push_back(unique_ptr<sphere>(new sphere(0, 3 * h / 4, w / 2, w / 4, color(1, 0, 0), surface_type::MIRROR)));
    scene.push_back(unique_ptr<sphere>(new sphere(2 * w / 3, 5 * h / 6, w / 6, w / 6, color(0, 1, 0), surface_type::GLOSSY, 1.2f)));
    scene.push_back(unique_ptr<sphere>(new sphere(0, 0, 4 * w, h, color(0, 0, 1), surface_type::GLOSSY, 1.2f)));



//    scene.push_back(unique_ptr<sphere>(new sphere(w / 2, 3 * h / 4, w / 2, w / 4, color(1, 0, 0), false, false, false, 1.1f, true)));          // 0

//    scene.push_back(unique_ptr<sphere>(new sphere(w / 4, 4 * h / 5, w, w / 5, color(0, 1, 0))));
//    scene.push_back(unique_ptr<sphere>(new sphere((256 + 170) / 512.0f * w, 4 * h / 5, 100.0f / 512 * w, w / 5, color(0, 0, 1), false, false, false, 1.15f, true)));
    //scene.push_back(unique_ptr<sphere>(new sphere(3 * w / 2, h / 2, w / 2, w / 2, color(0, 0, 1))));

//    scene.push_back(unique_ptr<plane>(new plane(vec3f(-1, 0, 0), vec3f(w + MARGIN, 0, 0), color(1, 1, 0))));                            // 6
//    scene.push_back(unique_ptr<plane>(new plane(vec3f(1, 0, 0), vec3f(0 - MARGIN, 0, 0), color(1, 1, 0))));                             // 7
    scene.push_back(unique_ptr<plane>(
            new plane(vec3f(0, -1, 0), vec3f(0, h, 0), surface_type::SIMPLE)));                            // 8
//    scene.push_back(unique_ptr<plane>(new plane(vec3f(0, 1, 0), vec3f(0, 0 - MARGIN, 0), color(1, 0, 1))));                             // 9
//    scene.push_back(unique_ptr<plane>(new plane(vec3f(0, 0, -1), vec3f(0, 0, w + MARGIN), color(0, 1, 1))));                            // 10
//    scene.push_back(unique_ptr<plane>(new plane(vec3f(0, 0, 1), vec3f(0, 0, 0 - MARGIN), color(0, 1, 1))));                             // 11


    view = vec3f(w / 2, h / 2, -w / 2);

//    get_pixel(w / 2, h / 2);
//    return 0;

    int start_time = time(0);
    string path = "E:\\C++ Projects\\RayTracer\\";
//    string path = "C:\\Users\\slava\\ClionProjects\\RayTracer\\RayTracer\\";

    thread t1(perform, 0, 8), t2(perform, 1, 8), t3(perform, 2, 8), t4(perform, 3, 8),
            t5(perform, 4, 8), t6(perform, 5, 8), t7(perform, 6, 8), t8(perform, 7, 8);
    t1.join();
    t2.join();
    t3.join();
    t4.join();
    t5.join();
    t6.join();
    t7.join();
    t8.join();
//    thread t1(perform, 0, 4), t2(perform, 1, 4), t3(perform, 2, 4), t4(perform, 3, 4);
//    t1.join();
//    t2.join();
//    t3.join();
//    t4.join();
//    perform(0, 1);

    ofstream fout(path + "out.ppm", ios::out | ios::binary);
    fout << "P6 " << w << " " << h << " 255 ";
    fout.write((char*) cl, w * h * 3);
    fout.close();
    cout << float(time(0) - start_time) << " seconds" << endl;
    system(("convert \"" + path + "out.ppm\" \"" + path + "out.png\"").c_str());
    return 0;
}