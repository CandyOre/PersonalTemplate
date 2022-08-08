#include <cmath>
#include <climits>
#include <algorithm>
#include <vector>
using namespace std;

using ld = long double;
const double eps = 1e-10;

using T = ld;

inline int sgn(T x) {return abs(x) <= eps ? 0 : (x > 0 ? 1 : -1);}

struct Point {
    T x, y;
    Point(){};
    Point(T xx, T yy) { x = xx, y = yy; }

    // bool order;
    // Point(T xx, T yy) { x = xx, y = yy, order = (yy > 0 || (yy == 0 && xx < 0)); }

    // long double arg;
    // Point(T xx, T yy) { x = xx, y = yy, arg = atan2l(yy * 1.0L, xx * 1.0L); }

    Point operator + (const Point& v) const {
        return Point(x + v.x, y + v.y);
    }
    Point operator - (const Point& v) const {
        return Point(x - v.x, y - v.y);
    }
    Point operator * (const T& k) const {
        return Point(x * k, y * k);
    }
    T operator * (const Point& v) const {
        return x * v.x + y * v.y;
    }
    T operator ^ (const Point& v) const {
        return x * v.y - y * v.x;
    }
    bool operator < (const Point& v) const {
        return sgn(x - v.x) ? sgn(x - v.x) == -1 : sgn(y - v.y) == -1;
    }
    bool operator == (const Point& v) const {
        return !sgn(x - v.x) && !sgn(y - v.y);
    }

    T len2 () {
        return x * x + y * y;
    }
    ld len () {
        return sqrt(x * x + y * y);
    }
    Point unit () {
        return *this * (1.0 / this->len());
    }
    T dis2 (const Point& v) {
        return (v - *this).len2();
    }
    ld dis (const Point& v) {
        return (v - *this).len();
    }

    T proj2(const Point& v) {
        return *this * v;
    }
    ld angle(Point& v) {
        return acosl(clamp((long double)(*this * v) / (this->len() * v.len()), -1.0L, 1.0L));
    }
};

struct Segment {
    Point u, v;
    Segment(){};
    Segment(Point uu, Point vv) { u = uu, v = vv; } 

    // careful with double precisions
    bool have(const Point& p) {
        return sgn(p - u ^ p - v) == 0 && sgn((p - u) * (p - v)) <= 0;
    }

    // p in the left side of segment
    int toLeft(const Point& p) {
        return sgn(v - u ^ p - u);
    }

    bool parallel(const Segment& s) {
        return !sgn(u - v ^ s.u - s.v);
    }

    bool isInter(const Segment& s) {
        return this->toLeft(s.u) * this->toLeft(s.v) <= 0;
    }
    // double / Not parallel
    Point inter(const Segment& s) {
        return u + (v - u) * ((s.v - s.u ^ u - s.u) / (v - u ^ s.v - s.u));
    }

    // directed area * 2
    T area2To(Point& p) {
        return p - u ^ p - v;
    }
    ld len() {
        return v.dis(u);
    }
    ld disTo(Point& p) {
        return abs(area2To(p)) / len();
    }
    ld minDisTo(Point& p) {
        ld dis = min(u.dis(p), v.dis(p));
        Point dir = v - u;
        dir = {dir.y, -dir.x};
        Segment vertical(p, p + dir);
        if(vertical.isInter(*this)) {
            dis = min(dis, disTo(p));
        }
        return dis;
    }
};

struct Polygon {
    // counter clockwise
    vector<Point> p;

    Polygon(int n) { p.resize(n); }

    inline size_t nxt(const size_t i) const {return i==p.size()-1?0:i+1;}
    inline size_t pre(const size_t i) const {return i==0?p.size()-1:i-1;}

    // on the line: INT_MIN
    // inside: >= 1 (winding number)
    // outside: 0
    // O(n)
    int have(const Point& a) {
        int cnt = 0;
        for(int i = 0; i < p.size(); i++) {
            Segment s(p[i], p[(i + 1) % p.size()]);
            if(s.have(a)) return INT_MIN;
            if(sgn(s.u.y - s.v.y) == 0) continue;
            if(sgn(s.u.y - s.v.y) < 0 && s.toLeft(a) < 0) continue;
            if(sgn(s.u.y - s.v.y) > 0 && s.toLeft(a) > 0) continue;
            if(sgn(s.u.y - a.y) < 0 && sgn(s.v.y - a.y) >= 0) cnt++;
            if(sgn(s.u.y - a.y) >= 0 && sgn(s.v.y - a.y) < 0) cnt--; 
        }
        return cnt;
    }

    T area2() {
        T sum = 0;
        for(int i = 0; i < p.size(); i++) {
            sum += p[i] ^ p[nxt(i)];
        }
        return sum;
    }

    // circumference
    ld circ() {
        ld sum = 0;
        for(int i = 0; i < p.size(); i++) {
            sum += p[i].dis(p[nxt(i)]);
        }
        return sum;
    }
};

struct argcmp {
    // -x counter clockwise
    int order(const Point& a) const {
        if(sgn(a.y) < 0) return 1;
        if(sgn(a.y) > 0) return 4;
        if(sgn(a.x) < 0) return 5;
        if(sgn(a.x) > 0) return 3;
        return 2;
    }
    bool operator () (const Point& a, const Point& b) const {
        int oa = order(a), ob = order(b);
        if(oa == ob) {
            const T t = a ^ b;
            if(sgn(t) == 0) return sgn(a * a - b * b) < 0;
            return sgn(t) > 0;
        }
        return oa < ob;
    }
} cmp;

struct Convex: Polygon {
    // O(nlogn)
    Convex(vector<Point>& p): Polygon(p.size()) {
        vector<Point> st;
        sort(p.begin(), p.end());

        auto check = [&](Point& u) {
            Point back1 = st.back(), back2 = *prev(st.end(), 2);
            return Segment(back2, back1).toLeft(u) <= 0;
        };
        for(Point& u: p) {
            while(st.size() > 1 && check(u)) st.pop_back();
            st.push_back(u);
        }
        size_t cur = st.size();
        for(int i = p.size() - 2; i > -1; i--) {
            while(st.size() > cur && check(p[i])) st.pop_back();
            st.push_back(p[i]);
        }
        st.pop_back();
        this->p = st;
    }

    // O(nlogn)
    Convex operator + (const Convex& c) const {
        vector<Segment> e1(p.size()), e2(c.p.size()), e(p.size() + c.p.size());
        for(int i = 0; i < p.size(); i++) e1[i] = {p[i], p[nxt(i)]};
        for(int i = 0; i < c.p.size(); i++) e2[i] = {c.p[i], c.p[c.nxt(i)]};

        auto scmp = [](Segment& a, Segment& b){return cmp(a.v - a.u, b.v - b.u);};
        rotate(e1.begin(), min_element(e1.begin(), e1.end(), scmp), e1.end());
        rotate(e2.begin(), min_element(e2.begin(), e2.end(), scmp), e2.end());
        merge(e1.begin(), e1.end(), e2.begin(), e2.end(), e.begin(), scmp);

        vector<Point> res;
        res.reserve(p.size() + c.p.size());
        auto check = [&](Point& u) {
            Point back1 = res.back(), back2 = *prev(res.end(), 2);
            return !Segment(back2, back1).toLeft(u) && sgn((back1 - back2) * (u - back1)) >= 0;
        };
        auto u = e1[0].u + e2[0].u;
        for(const auto& s : e) {
            res.push_back(u);
            u = u + s.v - s.u;
            while(res.size() > 1 && check(u)) res.pop_back();
        }
        return {res};
    }

    // on the line: -1
    // inside: 1
    // outside: 0
    // O(logn)
    int have(const Point& a) {
        if(p.size() == 1) return a == p[0] ? -1 : 0;
        if(p.size() == 2) return Segment(p[0], p[1]).have(a) ? -1 : 0;
        if(a == p[0]) return -1;
        if(Segment(p[0], p[1]).toLeft(a) == -1 || Segment(p[0], p.back()).toLeft(a) == 1) return 0;

        auto cmp = [&](const Point& u, const Point& v) {return Segment(p[0], u).toLeft(v) == 1;};
        size_t i = lower_bound(p.begin() + 1, p.end(), a, cmp) - p.begin();
        if(i == 1) return Segment(p[0], p[i]).have(a) ? -1 : 0;
        if(i == p.size() - 1 && Segment(p[0], p[i]).have(a)) return -1;
        if(Segment(p[i-1], p[i]).have(a)) return -1;
        return Segment(p[i-1], p[i]).toLeft(a) > 0;
    }

    // O(n)
    template<typename F>
    void rotcaliper(F& func) {
        for(int i = 0, j = 1; i < p.size(); i++) {
            auto ii = this->nxt(i);
            func(p[i], p[ii], p[j]);
            Segment seg(p[i], p[ii]);
            while(seg.area2To(p[j]) <= seg.area2To(p[nxt(j)])) {
                j = nxt(j);
                func(p[i], p[ii], p[j]);
            }
        }
    }

    // O(n)
    T diameter2() {
        if (p.size() == 1) return 0;
        if (p.size() == 2) return p[0].dis2(p[1]);
        T ans = 0;
        auto func = [&](Point& u, Point& v, Point& w) {
            ans = max({ans, w.dis2(u), w.dis2(v)});
        };
        rotcaliper(func);
        return ans;
    }

    // O(nlogn)
    T minDisTo(Convex& c) {
        Convex m = *this + c;
        ld ans = LLONG_MAX;
        Point o(0, 0);
        for(int i = 0; i < m.p.size(); i++) {
            ans = min(ans, Segment(m.p[i], m.p[m.nxt(i)]).minDisTo(o));
        }
        return ans;
    }
};