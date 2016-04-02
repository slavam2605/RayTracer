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
        return hsl2hsv(hsv((atan2f(y1, x1) + PI) * 180 / PI, sqrtf(x1 * x1 + y1 * y1), z1));
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

struct object {
    color cl;
    float on;

    virtual float intersect(const ray& r) = 0;

    virtual color get_color(const ray& r, const vec3f& p, float cn, int bunch_depth, int depth) = 0;

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
        color light_c;
        int n = NN;
        for (int i = 0; i < n; i++) {
            vec3f dir = rand_dir();
            if (dir * norm < 0) dir *= -1;
            light_c += ::get_color(ray(p + norm, dir), cn, bunch_depth + 1, depth + 1);
        }
        light_c /= n;
        color c_difr = cl * light_c;
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

};

vec3f view(256 / 2, 256 / 2, -256);
vector<unique_ptr<object>> scene;

struct sphere : object {
    vec3f c;
    float R;
    bool is_light;
    bool is_mirror;
    bool is_transparent;
    bool is_glossy;

    sphere(const vec3f& c, float R, color cl = color(1, 1, 1), bool is_light = false, bool is_mirror = false,
           bool is_transparent = false, float on = 1, bool is_glossy = false)
            : c(c), R(R), is_light(is_light), is_mirror(is_mirror), is_transparent(is_transparent), is_glossy(is_glossy) {
        this->on = on;
        this->cl = cl;
    }

    sphere(float x, float y, float z, float R, color cl = color(1, 1, 1), bool is_light = false, bool is_mirror = false,
           bool is_transparent = false, float on = 1, bool is_glossy = false)
            : c(x, y, z), R(R), is_light(is_light), is_mirror(is_mirror), is_transparent(is_transparent), is_glossy(is_glossy) {
        this->on = on;
        this->cl = cl;
    }

    sphere() { }

    float intersect(const ray& r) {
        vec3f dr = c - r.r0;
        float D = sqr(r.s * dr) + sqr(R) - dr * dr;
        if (D < 1e-5) return 2 * INF;
        float result = r.s * dr - sqrtf(D);
        if (result < 0) {
            result = r.s * dr + sqrtf(D);
            if (result < 0) return 2 * INF; else return result;
        } else return result;
    }

    vec3f get_norm(const ray& r, const vec3f& p) {
        vec3f norm = (p - c) / R;
        if (r.s * norm > 0) norm *= -1;
        return norm;
    }

    bool is_inside(const ray& r, const vec3f& p) {
        return r.s * (p - c) > 0;
    }

    color get_color(const ray& r, const vec3f& p, float cn, int bunch_depth, int depth) {
        if (depth > DEPTH_LIMIT) return color(0, 0, 0);
        if (is_light) return cl;
        vec3f norm = get_norm(r, p);
        if (is_mirror) {
            return reflect(r, p, cn, bunch_depth, depth, norm);
        }
        if (is_transparent) {
            float n2 = is_inside(r, p) ? on : 1;
            return full_refract(r, p, cn, n2, bunch_depth, depth, norm);
        }
        if (bunch_depth > BUNCH_LIMIT) return color(0, 0, 0);
        if (is_glossy) {
            return glossy(r, p, cn, bunch_depth, depth, norm);
        }
        return diffuse(p, cn, bunch_depth, depth, norm);
    }
};

struct plane : object {
    vec3f norm;
    float D;
    bool is_mirror;

    plane(const vec3f& norm, float D, const color& cl = color(1, 1, 1), bool is_mirror = false)
            : norm(norm), D(D), is_mirror(is_mirror) {
        this->cl = cl;
    }

    plane(const vec3f& norm, const vec3f& p, const color& cl = color(1, 1, 1), bool is_mirror = false)
            : norm(norm), D(norm * (-1) * p), is_mirror(is_mirror) {
        this->cl = cl;
    }

    float intersect(const ray& r) {
        float result = -(D + r.r0 * norm) / (r.s * norm);
        if (result < 1e-5) return 2 * INF;
        return result;
    }

    color get_color(const ray& r, const vec3f& p, float cn, int bunch_depth, int depth) {
        if (depth > DEPTH_LIMIT) return color(0, 0, 0);
        if (is_mirror) {
            return reflect(r, p, cn, bunch_depth, depth, norm);
        }
        if (bunch_depth > BUNCH_LIMIT) return color(0, 0, 0);
        color cll;
        if (norm.y != 0) {
            if (((int) floor(p.x / 100) + (int) floor(p.z / 100)) & 1)
                cll = color(1, 0, 0);
            else
                cll = color(1, 1, 1);
        }
        color light_c;
        int n = NN;
        for (int i = 0; i < n; i++) {
            vec3f dir = rand_dir();
            if (dir * norm < 0) dir *= -1;
            light_c += ::get_color(ray(p + norm, dir), cn, bunch_depth + 1, depth + 1);
        }
        light_c /= n;
        return cll * light_c;
        //if (((int) floor(p.x / 100) + (int) floor(p.z / 100)) & 1) return color(1, 0, 0); else return color(1, 1, 1);
        //return cl * diffuse(p, cn, bunch_depth, depth, norm);
    }
};

const int w = 512;
const int h = 512;

color get_color(float x, float y) {
    int N = 0;
    float range = 2.5f;
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
//    return hsv2rgb(rgb2hsv(c1).mix(rgb2hsv(c2)));
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
        if (i % 10 == 0) {
            progress[k] = 100 * (i + 1 - start) / (end - start);
            print(n);
        }
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


    scene.push_back(unique_ptr<sphere>(new sphere(w / 2, -h, w / 2, w / 5, color(100, 100, 100), true)));          // 1

/*

    scene.push_back(unique_ptr<sphere>(new sphere(w / 2, 3 * h / 4, w / 2, w / 4, color(0.5f, 0.5f, 0), false, false, false, 1.2f, true)));          // 0
    scene.push_back(unique_ptr<sphere>(new sphere(w, 3 * h / 4, w / 2, w / 4, color(1, 0, 0), false, false, true, 1.2f)));          // 0
    scene.push_back(unique_ptr<sphere>(new sphere(0, 3 * h / 4, w / 2, w / 4, color(1, 0, 0), false, true)));          // 0

    scene.push_back(unique_ptr<sphere>(new sphere(2 * w / 3, 5 * h / 6, w / 6, w / 6, color(0, 1, 0), false, false, false, 1.2f, true)));          // 0

    scene.push_back(unique_ptr<sphere>(new sphere(0, 0, 4 * w, h, color(0, 0, 1), false, false, false, 1.2f, true)));          // 0
*/


//    scene.push_back(unique_ptr<sphere>(new sphere(w / 2, 3 * h / 4, w / 2, w / 4, color(1, 0, 0), false, false, false, 1.1f, true)));          // 0

//    scene.push_back(unique_ptr<sphere>(new sphere(w / 4, 4 * h / 5, w, w / 5, color(0, 1, 0))));
//    scene.push_back(unique_ptr<sphere>(new sphere((256 + 170) / 512.0f * w, 4 * h / 5, 100.0f / 512 * w, w / 5, color(0, 0, 1), false, false, false, 1.15f, true)));
    //scene.push_back(unique_ptr<sphere>(new sphere(3 * w / 2, h / 2, w / 2, w / 2, color(0, 0, 1))));

//    scene.push_back(unique_ptr<plane>(new plane(vec3f(-1, 0, 0), vec3f(w + MARGIN, 0, 0), color(1, 1, 0))));                            // 6
//    scene.push_back(unique_ptr<plane>(new plane(vec3f(1, 0, 0), vec3f(0 - MARGIN, 0, 0), color(1, 1, 0))));                             // 7
    scene.push_back(unique_ptr<plane>(
            new plane(vec3f(0, -1, 0), vec3f(0, h, 0), color(1, 0, 1))));                            // 8
//    scene.push_back(unique_ptr<plane>(new plane(vec3f(0, 1, 0), vec3f(0, 0 - MARGIN, 0), color(1, 0, 1))));                             // 9
//    scene.push_back(unique_ptr<plane>(new plane(vec3f(0, 0, -1), vec3f(0, 0, w + MARGIN), color(0, 1, 1))));                            // 10
//    scene.push_back(unique_ptr<plane>(new plane(vec3f(0, 0, 1), vec3f(0, 0, 0 - MARGIN), color(0, 1, 1))));                             // 11


    view = vec3f(w / 2, h / 2, -w / 2);

//    get_pixel(w / 2, h / 2);
//    return 0;

    int start_time = time(0);
//    string path = "E:\\C++ Projects\\RayTracer\\";
    string path = "C:\\Users\\slava\\ClionProjects\\RayTracer\\RayTracer\\";

//    thread t1(perform, 0, 8), t2(perform, 1, 8), t3(perform, 2, 8), t4(perform, 3, 8),
//            t5(perform, 4, 8), t6(perform, 5, 8), t7(perform, 6, 8), t8(perform, 7, 8);
//    t1.join();
//    t2.join();
//    t3.join();
//    t4.join();
//    t5.join();
//    t6.join();
//    t7.join();
//    t8.join();
    thread t1(perform, 0, 4), t2(perform, 1, 4), t3(perform, 2, 4), t4(perform, 3, 4);
    t1.join();
    t2.join();
    t3.join();
    t4.join();
//    perform(0, 1);

    ofstream fout(path + "out.ppm", ios::out | ios::binary);
    fout << "P6 " << w << " " << h << " 255 ";
    fout.write((char*) cl, w * h * 3);
    fout.close();
    cout << float(time(0) - start_time) << " seconds" << endl;
    system(("convert \"" + path + "out.ppm\" \"" + path + "out.png\"").c_str());
    return 0;
}