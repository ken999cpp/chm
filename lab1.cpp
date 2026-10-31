#include <cstdio>
#include <cmath>
#include <vector>
#ifdef _WIN32
#include <windows.h>
#endif
using namespace std;

constexpr double EPS = 1e-3;
constexpr double TOL = EPS / 2;
constexpr int MAX_IT = 1000;

typedef double (*Fn)(double);
struct Eq {
    const char* name;
    Fn f; // функция
    Fn df; // ее производная
    Fn d2f; // ее вторая производная
}; // уравнение f(x) = 0

struct Seg { 
    double a;
    double b;
}; // отрезок [a; b]

double horner(const vector<double> &p, double x, vector<double> &q) {
    q.assign(p.size(), 0);
    q[0] = p[0];
    for (size_t i = 1; i < p.size(); i++) {
        q[i] = q[i - 1] * x + p[i];
    }
    return q.back();
}

void printPoly(const vector<double> &p) {
    int n = (int)p.size() - 1;
    for (int i = 0; i <= n; i++) {
        printf("%+g*x^%d ", p[i], n - i);
    }
    printf("\n");
}

int upperBound(const vector<double> &p) {
    vector<double> q;
    for (int c = 0; c <= 1000; c++) {
        horner(p, c, q);
        bool ok = true;
        for (double v : q) {
            if (v < 0) {
                ok = false;
            }
        }
        printf("   c = %-3d строка Горнера:", c);
        for (double v : q) {
            printf(" %9.4f", v);
        }
        printf("  %s\n", ok ? "<- все >= 0" : "");
        if (ok) {
            return c;
        }
    }
    return 1000;
}

Seg polyBounds(const vector<double> &p) {
    int n = (int)p.size() - 1;
    printf("  Верхняя граница: P(x) = "); printPoly(p);
    int up = upperBound(p);
    vector<double> q(p.size());
    for (int i = 0; i <= n; i++) {
        q[i] = p[i] * ((i % 2) ? -1 : 1) + 0.0;
    }
    printf("  Нижняя граница: Q(x) = (-1)^n P(-x) = "); printPoly(q);
    int lo = upperBound(q);
    printf("  => все действительные корни лежат в [%d; %d]\n", -lo, up);
    return {(double)-lo, (double)up};
}

vector<Seg> separate(const Eq &e, double A, double B, double h) {
    vector<Seg> r;
    int n = (int)lround((B - A) / h);
    for (int i = 0; i < n; i++) {
        double a = A + i * h, b = A + (i + 1) * h;
        double fa = e.f(a), fb = e.f(b);
        if (fa * fb < 0) r.push_back({a, b});
    }
    printf("  Табулирование на [%g; %g] с шагом %g: найдено отрезков с корнем: %d\n", A, B, h, (int)r.size());
    for (auto& s : r) {
        printf("    [%8.4f; %8.4f]  f(a) = %+.5f  f(b) = %+.5f\n", s.a, s.b, e.f(s.a), e.f(s.b));
    }
    return r;
}

Seg narrow(const Eq &e, Seg s) {
    for (int i = 0; i < 60 && (e.df(s.a) * e.df(s.b) <= 0 || e.d2f(s.a) * e.d2f(s.b) <= 0); i++) {
        double c = (s.a + s.b) / 2;
        if (e.f(s.a) * e.f(c) <= 0) {
            s.b = c;
        } else {
            s.a = c;
        }
    }
    return s;
}

void bounds1(const Eq &e, Seg s, double &m1, double &M1) {
    double x = fabs(e.df(s.a)), y = fabs(e.df(s.b));
    m1 = fmin(x, y); M1 = fmax(x, y);
}

double r3(double x) { return round(x * 1000) / 1000 + 0.0; }

void check(const Eq &e, double x) {
    double l = e.f(x - EPS), r = e.f(x + EPS);
    printf("  Проверка: f(%.6f) = %+.2e;  f(x-eps) = %+.2e, f(x+eps) = %+.2e  => %s\n",
           x, e.f(x), l, r, l * r <= 0 ? "корень на отрезке [x-eps; x+eps], верно" : "ОШИБКА");
}

double iteration(const Eq &e, Seg s) {
    s = narrow(e, s);
    double m1, M1; bounds1(e, s, m1, M1);
    double K = (e.df(s.a) > 0 ? 1 : -1) * M1, q = 1 - m1 / M1;
    printf("  Метод итерации на [%.5f; %.5f]: phi(x) = x - f(x)/(%.5f), q = max|phi'| = %.5f < 1\n", s.a, s.b, K, q);
    double x = s.a, stop = q > 0 ? TOL * (1 - q) / q : INFINITY;
    printf("    n  x_n          f(x_n)\n    0  %.6f  %+.6f\n", x, e.f(x));
    for (int n = 1; n <= MAX_IT; n++) {
        double nx = x - e.f(x) / K;
        printf("  %3d  %.6f  %+.6f\n", n, nx, e.f(nx));
        bool done = fabs(nx - x) <= stop;
        x = nx;
        if (done) {
            break;
        }
    }
    printf("  Остановка при |x_n - x_(n-1)| <= (eps/2)(1-q)/q = %.2e\n  Корень: x = %.6f ~ %.3f\n", stop, x, r3(x));
    check(e, x);
    return x;
}

double chord(const Eq &e, Seg s) {
    s = narrow(e, s);
    double m1, M1; bounds1(e, s, m1, M1);
    double c = e.f(s.a) * e.d2f(s.a) > 0 ? s.a : s.b, x = c == s.a ? s.b : s.a;
    printf("  Метод хорд на [%.5f; %.5f]: неподвижный конец c = %.5f, m1 = %.5f\n", s.a, s.b, c, m1);
    printf("    n  x_n          f(x_n)      |f|/m1\n    0  %.6f  %+.6f  %.2e\n", x, e.f(x), fabs(e.f(x)) / m1);
    for (int n = 1; n <= MAX_IT && fabs(e.f(x)) / m1 > TOL; n++) {
        x = x - e.f(x) * (c - x) / (e.f(c) - e.f(x));
        printf("  %3d  %.6f  %+.6f  %.2e\n", n, x, e.f(x), fabs(e.f(x)) / m1);
    }
    printf("  Корень: x = %.6f ~ %.3f\n", x, r3(x));
    check(e, x);
    return x;
}

double newton(const Eq &e, Seg s) {
    s = narrow(e, s);
    double m1, M1; bounds1(e, s, m1, M1);
    double x = e.f(s.a) * e.d2f(s.a) > 0 ? s.a : s.b;
    printf("  Метод касательных на [%.5f; %.5f]: x0 = %.5f (f*f'' > 0), m1 = %.5f\n", s.a, s.b, x, m1);
    printf("    n  x_n          f(x_n)      |f|/m1\n    0  %.6f  %+.6f  %.2e\n", x, e.f(x), fabs(e.f(x)) / m1);
    for (int n = 1; n <= MAX_IT && fabs(e.f(x)) / m1 > TOL; n++) {
        x = x - e.f(x) / e.df(x);
        printf("  %3d  %.6f  %+.6f  %.2e\n", n, x, e.f(x), fabs(e.f(x)) / m1);
    }
    printf("  Корень: x = %.6f ~ %.3f\n", x, r3(x));
    check(e, x);
    return x;
}

double combined(const Eq &e, Seg s) {
    s = narrow(e, s);
    double t = e.f(s.a) * e.d2f(s.a) > 0 ? s.a : s.b, h = t == s.a ? s.b : s.a;
    printf("  Комбинированный метод на [%.5f; %.5f]\n", s.a, s.b);
    printf("    n  касательная  хорда        |разность|\n    0  %.6f    %.6f    %.2e\n", t, h, fabs(t - h));
    for (int n = 1; n <= MAX_IT && fabs(t - h) > 2 * TOL; n++) {
        double nh = h - e.f(h) * (t - h) / (e.f(t) - e.f(h));
        t = t - e.f(t) / e.df(t);
        h = nh;
        printf("  %3d  %.6f    %.6f    %.2e\n", n, t, h, fabs(t - h));
    }
    double x = (t + h) / 2;
    printf("  Корень: x = (t + h)/2 = %.6f ~ %.3f\n", x, r3(x));
    check(e, x);
    return x;
}

double ft(double x) {
    return x * x - 2;
}

double dft(double x) {
    return 2 * x;
}

double d2ft(double) {
    return 2;
}

bool selfTest() {
    printf("тест: f(x) = x^2 - 2 = 0, точный корень sqrt(2) = %.6f =====\n", sqrt(2.0));
    Eq e = {"x^2 - 2 = 0", ft, dft, d2ft};
    Seg s = {1, 2};
    double r[4] = {iteration(e, s), chord(e, s), newton(e, s), combined(e, s)};
    bool ok = true;
    for (double v : r) {
        ok = ok && fabs(v - sqrt(2.0)) <= EPS;
    }
    printf("тест %s\n\n", ok ? "пройден: все 4 метода дали sqrt(2) с точностью 0.001" : "не пройден");
    return ok;
}

double f21(double x) {
    return x * x * x - 0.1 * x * x + 0.4 * x + 2;
}

double df21(double x) {
    return 3 * x * x - 0.2 * x + 0.4;
}

double d2f21(double x) {
    return 6 * x - 0.2; 
}

double f22(double x) {
    return (x - 1) * (x - 1) - 0.5 * exp(x); 
}

double df22(double x) {
    return 2 * (x - 1) - 0.5 * exp(x); 
}

double d2f22(double x) {
    return 2 - 0.5 * exp(x); 
}

double f3a(double x) {
    return x + log10(x) - 0.5; 
}

double df3a(double x) {
    return 1 + 1 / (x * log(10.0)); 
}

double d2f3a(double x) {
    return -1 / (x * x * log(10.0)); 
}

double f3b(double x) {
    return x * x * x + 3 * x + 1; 
}

double df3b(double x) {
    return 3 * x * x + 3; 
}

double d2f3b(double x) {
    return 6 * x; 
}

double f3c(double x) {
    return x * x * x - 3 * x * x + 2.5; 
}

double df3c(double x) {
    return 3 * x * x - 6 * x; 
}

double d2f3c(double x) {
    return 6 * x - 6; 
}

int main() {
#ifdef _WIN32
    SetConsoleOutputCP(CP_UTF8);
#endif
    if (!selfTest()) {
        return 1;
    }

    Eq e21 = {"x^3 - 0.1x^2 + 0.4x + 2 = 0", f21, df21, d2f21};
    Eq e22 = {"(x-1)^2 = 0.5 e^x", f22, df22, d2f22};
    Eq e3a = {"x + lg(x) = 0.5", f3a, df3a, d2f3a};
    Eq e3b = {"x^3 + 3x + 1 = 0", f3b, df3b, d2f3b};
    Eq e3c = {"x^3 - 3x^2 + 2.5 = 0", f3c, df3c, d2f3c};

    printf("№ 2.1. %s =====\n", e21.name);
    Seg b = polyBounds({1, -0.1, 0.4, 2});
    for (auto& s : separate(e21, b.a, b.b, 0.1)) {
        iteration(e21, s);
    }

    printf("\n№ 2.2. %s =====\n", e22.name);
    for (auto& s : separate(e22, 0, 3, 0.1)) {
        iteration(e22, s);
    }

    printf("\n№ 3а. Метод хорд. %s =====\n", e3a.name);
    for (auto& s : separate(e3a, 0.05, 2.05, 0.1)) {
        chord(e3a, s);
    }

    printf("\n№ 3б. Метод касательных. %s =====\n", e3b.name);
    Seg b3 = polyBounds({1, 0, 3, 1});
    for (auto& s : separate(e3b, b3.a, b3.b, 0.1)) {
        newton(e3b, s);
    }

    printf("\n№ 3в. Комбинированный метод. %s =====\n", e3c.name);
    Seg b4 = polyBounds({1, -3, 0, 2.5});
    for (auto& s : separate(e3c, b4.a, b4.b, 0.1)) {
        combined(e3c, s);
    }
    return 0;
}
