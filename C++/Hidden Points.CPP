#include<bits/stdc++.h>
using namespace std;

struct custHash {
    static uint64_t splitmix64(uint64_t x) {  
        x += 0x9e3779b97f4a7c15;
        x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9;
        x = (x ^ (x >> 27)) * 0x94d049bb133111eb;
        return x ^ (x >> 31);
    }

    size_t operator()(uint64_t x) const {
        static const uint64_t FIXED_RANDOM = chrono::steady_clock::now().time_since_epoch().count();
        return splitmix64(x + FIXED_RANDOM);
    }
};

vector<string> vec_splitter(string s) {
	s += ',';
	vector<string> res;
	while(!s.empty()) {
		res.push_back(s.substr(0, s.find(',')));
		s = s.substr(s.find(',') + 1);
	}
	return res;
}
void debug_out(
vector<string> __attribute__ ((unused)) args,
__attribute__ ((unused)) int idx, 
__attribute__ ((unused)) int LINE_NUM) { cerr << endl; } 
template <typename Head, typename... Tail>
void debug_out(vector<string> args, int idx, int LINE_NUM, Head H, Tail... T) {
	if(idx > 0) cerr << ", "; else cerr << "Line(" << LINE_NUM << ") ";
	stringstream ss; ss << H;
	cerr << args[idx] << " = " << ss.str();
	debug_out(args, idx + 1, LINE_NUM, T...);
}
#ifdef LOCAL
#define debug(...) debug_out(vec_splitter(#__VA_ARGS__), 0, __LINE__, __VA_ARGS__)
#else
#define debug(...) 42
#endif

#define FOR(i, a, b) for (int i = (a), _##i = (b); i < _##i; ++i)
#define FORD(i, a, b) for (int i = (a), _##i = (b); i >= _##i; --i)
#define REP(i, a) for (int i = 0, _##i = (a); i < _##i; ++i)
#define REPD(i,n) for(int i = (n)-1; i >= 0; --i)
#define all(x) (x).begin(), (x).end()
#define pb emplace_back

typedef long long LL;
#define double long double
#define sqrt sqrtl
const double EPS  = (double) 1e-9;


inline int cmp(double a, double b) {
    return (a < b - EPS) ? -1 : ((a > b + EPS) ? 1 : 0);
}

struct Point {
    double x, y;
    Point() { x = y = 0.0; }
    Point(double x, double y) : x(x), y(y) {}

    Point operator + (const Point& a) const { return Point(x+a.x, y+a.y); }
    Point operator - (const Point& a) const { return Point(x-a.x, y-a.y); }
    Point operator * (double k) const { return Point(x*k, y*k); }
    Point operator / (double k) const { return Point(x/k, y/k); }

    double operator * (const Point& a) const { return x*a.x + y*a.y; } // dot product
    double operator % (const Point& a) const { return x*a.y - y*a.x; } // cross product

    int cmp(Point q) const { if (int t = ::cmp(x,q.x)) return t; return ::cmp(y,q.y); }

    #define Comp(x) bool operator x (Point q) const { return cmp(q) x 0; }
    Comp(>) Comp(<) Comp(==) Comp(>=) Comp(<=) Comp(!=)
    #undef Comp

    Point conj() { return Point(x, -y); }
    double norm() { return x*x + y*y; }

    // Note: There are 2 ways for implementing len():
    // 1. sqrt(norm()) --> fast, but inaccurate (produce some values that are of order X^2)
    // 2. hypot(x, y) --> slow, but much more accurate
    // double len() { return sqrt(norm()); }
    double len() { return hypot(x, y); }
};

int n;
Point pts[100001];
Point mem[100001];
double x,y;
double d, d1, d2, d3;
double ds, ds1, ds2, ds3;
int pivot;
const double oo = (double) 1e19;

int cmpY(Point lhs, Point rhs) {
    return make_pair(lhs.y, lhs.x) < make_pair(rhs.y, rhs.x);
}

int cmpX(Point lhs, Point rhs) {
    return make_pair(lhs.x, lhs.y) < make_pair(rhs.x, rhs.y);
}

double closestpair(Point p[], int l, int r) {
    double res = oo;
    if (r - l <= 3) {
        for (int i = l; i <= r; i++) {
            for (int j = i + 1; j <= r; j++) {
                res = min(res, (p[i]-p[j]).norm() );
            }
        }
        sort(p + l, p + r + 1, cmpY);
    }
    else {
        int mid = (l + r) >> 1;
        // debug(l,r,mid);
        double xmid = p[mid].x;
        res = min(closestpair(p, l, mid), closestpair(p, mid + 1, r));
        inplace_merge(p + l, p + mid + 1, p + r + 1, cmpY);
        int cnt = 0;
        for (int i = l; i <= r; i++) {
            if ( (p[i].x - xmid) * (p[i].x - xmid) < res) mem[cnt++] = p[i];
        }
        for (int i = 0; i < cnt; i++) {
            for (int j = 1; j <= 8 && i + j < cnt; j++) {
                res = min( res, (mem[i] - mem[i + j]).norm() );
            }
        }
    }
    return res;
}

bool onLine(int pos){
    if( !cmp( d - d1 - d2, 0 )  ){ // pos nam giua
        pts[pos] = Point(d1, 0);
        return true;
    }
    if( !cmp( d1 - d2 - d, 0 ) ){  // 1 nam giua
        pts[pos] = Point(d1, 0);
        return true;
    } 
    if( !cmp( d2 - d1 - d, 0 ) ){  // 0 nam giua
        pts[pos] = Point(-d1, 0);
        return true;
    }
    return false;
}

Point getThirdVertex(){
    double x = (ds + ds1 - ds2) / d / 2.;    
    double y = sqrt(ds1 - x*x);
    debug(x,y);
    return Point(x,y);
}

void solve(){
    cin >> n;
    pivot=0;
    pts[0] = Point(0,0);
    cout << "? 1 2" << endl;
    cin >> ds; d = sqrt(ds);
    pts[1] = Point(d,0);
    FOR(i,2,n){
        cout << "? 1 " << i+1 << endl; 
        cin >> ds1; d1 = sqrt(ds1);
        cout << "? 2 " << i+1 << endl; 
        cin >> ds2; d2 = sqrt(ds2);
        if( !onLine(i) ){
            pts[i] = getThirdVertex();
            if(!pivot) {
                pivot = i;
            }else{
                cout << "? " << pivot + 1 << " " << i+1 << endl; 
                cin >> ds3;
                Point conj = pts[i].conj();
                if( abs( ds3 - (pts[i] - pts[pivot] ).norm() )  > abs( ds3 - (conj - pts[pivot] ).norm() ) ){
                    pts[i]=conj;
                }
                /*if( (pts[i] - pts[pivot] ).norm() > 10 ){
                    debug(ds3, (pts[i] - pts[pivot] ).norm());
                    debug(ds3, (pts[i].conj() - pts[pivot] ).norm());
                }*/
            }
        }
    }
    
    sort(pts, pts + n, cmpX);
    cout << "! " << setprecision(0) << fixed << closestpair(pts,0,n-1) << endl;
    /*
    double d = 1e20;
    FOR(i,0,n) FOR(j,i+1,n) d = min(d, (pts[i].x-pts[j].x)*(pts[i].x-pts[j].x) + (pts[i].y - pts[j].y)*(pts[i].y - pts[j].y ) ) ;
    cout << "! " << setprecision(0) << fixed << d << endl;
    */
}

int ntest;
int main(){
	ios::sync_with_stdio(false);		
    cin.tie(NULL);
    cin >> ntest;
    REP(i,ntest){
        solve();
    }
	return 0;
}


