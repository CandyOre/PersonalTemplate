#include <cstdio>
#include <cmath>
#include <cstring>
#include <climits>
#include <iostream>
#include <iomanip>
#include <chrono>
#include <numeric>
#include <functional>
#include <algorithm>
#include <vector>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <map>
#include <set>
#include <queue>
#include <stack>
using namespace std;

// Define abbr. here.
#define endl "\n"
#define IL inline
#define fi first
#define se second
#define mk make_pair
#define pb push_back
#define SZ(x) (int)(x).size()
#define ALL(x) (x).begin(), (x).end()
#define dbg1(x) cout << #x << " = " << x << ", "
#define dbg2(x) cout << #x << " = " << x << endl
#define rep(i,n) for (int i=0;i<(int)(n);++i)
#define rep1(i,n) for (int i=1;i<=(int)(n);++i)
using ll = long long;
using ld = long double;
using ull = unsigned long long;
using pii = pair<int,int>;
using umii = unordered_map<int,int>;
using v1i = vector<int>;
using v1l = vector<ll>;
using v1b = vector<bool>;

// Define const. here.
// #define INT_MAX 0x7fffFFFF
// #define LLONG_MAX 0x7fffFFFFffffFFFF
const long double pi = acosl(-1);
const ll mod = 998244353;
const double eps = 1e-10;

#ifdef ONLINE_JUDGE
const bool OJ = 1;
#else
const bool OJ = 0;
#endif


// Define variables and functions here
int n;


void init(){

}

void solve(int case_no){
	
}


const bool MulCase = 0;

signed main(){
    std::ios::sync_with_stdio(false);
    std::cin.tie(0);std::cout.tie(0);
    
    cout << fixed << setprecision(8);

    if(!OJ){
        freopen("try.in", "r", stdin);
        freopen("try.out", "w", stdout);
    }

    int T = 1;
    if(MulCase) cin >> T;
    init();
    rep(i,T) solve(i);
    
    return 0;
}