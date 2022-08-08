using T = ll;

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

    T len2 () {
        return x * x + y * y;
    }
    double len () {
        return sqrt(x * x + y * y);
    }
    Point unit () {
        return *this * (1.0 / this->len());
    }
    T dis (const Point& v) {
        return (v - *this).len();
    }

    T proj2(const Point& v) {
        return *this * v;
    }
    double angle(Point& v) {
        return acosl(clamp((long double)(*this * v) / (this->len() * v.len()), -1.0L, 1.0L));
    }
};

struct Segment {
    Point u, v;
    Segment(Point uu, Point vv) { u = uu, v = vv; } 

    // careful with double precisions
    bool have(const Point& p) {
        return sgn(p - u ^ p - v) == 0 && sgn((p - u) * (p - v)) <= 0;
    }

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
};

struct Polygon {
    vector<Point> p;

    Polygon(int n) { p.resize(n); }

    inline size_t nxt(const size_t i) const {return i==p.size()-1?0:i+1;}
    inline size_t pre(const size_t i) const {return i==0?p.size()-1:i-1;}

    // on the line: INT_MIN
    // inside: 
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

[](pair<Point, int> aa, pair<Point, int> bb) {
    const auto& a = aa.first, b = bb.first;
    if(a.order == b.order) {
        const T t = a ^ b;
        return sgn(t) > 0;
    }
    return a.order < b.order;
}