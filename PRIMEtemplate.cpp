struct MT {
    ll x;
    MT() : x(0) {}
    MT(ll y): x(reduce(y)) {}

    static ll reduce(ll y) {return (y % mod + mod) % mod;}

    MT pow(ll p) const {
        if(p < 0) return pow(-p).inv();
        MT res = 1, a = *this;
        while(p) {
            if(p & 1) res *= a;
            p >>= 1, a *= a;
        }
        return res;
    }
    MT inv() const {return pow(mod - 2);}

    MT operator + (const MT& y) {ll t = x + y.x; return t < mod ? t : t - mod;}
    MT operator - (const MT& y) {ll t = x - y.x; return t >= 0 ? t : t + mod;}
    MT operator * (const MT& y) {return reduce(x * y.x);}
    MT operator / (const MT& y) {return y.inv()**this;}

    MT operator += (const MT& y) {x += y.x; return x = x < mod ? x : x - mod;}
    MT operator -= (const MT& y) {x -= y.x; return x = x >= 0 ? x : x + mod;}
    MT operator *= (const MT& y) {return x = reduce(x * y.x);}
    MT operator /= (const MT& y) {return *this *= y.inv();}

    bool operator == (const MT& y) const {return x == y.x;}

	ll operator () () {return reduce(x);}

    friend istream& operator >> (istream& i, MT& y) {ll x; i >> x; y = MT(x); return i;}
    friend ostream& operator << (ostream& o, const MT& y) {return o << reduce(y.x);}
};

template <typename T> struct Mat {
    T* elem;
    int n, m;
    Mat(int n, int m, bool diag = false) : n(n), m(m) {
        elem = (T*)malloc(n * m * sizeof(T));
        memset(elem, 0, n * m * sizeof(T));
        // only for n = m
        if(diag) rep(i, n) (*this)[i][i] = 1;
    }

    Mat(const Mat<T>& a) : n(a.n), m(a.m) {
        elem = (T*)malloc(n * m * sizeof(T));
        memcpy(elem, a.elem, n * m * sizeof(T));
    }

    // use m[i][j] to access m.elem[i][j]
    T* operator [] (int x) const {
        return &elem[x * m];
    }

    // below only for n = m
    Mat pow(ll p) const {
        if(p < 0) return pow(-p).inv();
        Mat ans(n, n, true), a(*this);
        while(p) {
            if(p & 1) ans = ans * a;
            p >>= 1, a = a * a;
        }
        return ans;
    }

    Mat inv() const {
		bool t;
		return inv(t);
	}

    Mat inv(bool& singular) const {
        singular = false;
        Mat x(*this), y(n, n, true);
        rep(i, n) {
            if(x[i][i] == 0) {
                singular = true;
                return y;
            } 
            T a = x[i][i].inv(); x[i][i] = 1;
            for(int j = i + 1; j < n; j++) x[i][j] *= a;
            for(int j = 0; j <= i; j++) y[i][j] *= a;
            for(int k = i + 1; k < n; k++) {
                T w = x[k][i];
                for(int j = i; j < n; j++) x[k][j] -= x[i][j] * w;
                for(int j = 0; j <= i; j++) y[k][j] -= y[i][j] * w;
            }
        }
        for(int i = n - 1; i > 0; i--) {
            for(int k = i - 1; k > -1; k--) {
                T w = x[k][i];
                for(int j = 0; j < n; j++) y[k][j] -= y[i][j] * w;
            }
        }
        return y;
    }
    
    // this->m = y.n
    Mat operator * (const Mat& y) const {
        Mat z(n, y.m);
        rep(i, n) rep(j, y.m) rep(k, m) {
            z[i][j] += (*this)[i][k] * y[k][j];
        }
        return z;
    }

    friend istream& operator >> (istream& is, Mat& y) {
        rep(i, y.n) rep(j, y.m) is >> y[i][j];
        return is;
    }

    friend ostream& operator << (ostream& o, const Mat& y) {
        rep(i, y.n) {rep(j, y.m) o << y[i][j] << " "; if(i < y.n - 1) cout << endl;}
        return o;
    }
};