#include <bits/stdc++.h>
using namespace std;

typedef long long int ll;
typedef long double ld;

#define N 100010
#define pi 3.14159265358979323846264338327950288419716939

#define mp(i, j) make_pair(i, j)
#define pb(i) push_back(i)

#define f(i, s, e, inc) for(ll i = s; i < e; i+=inc) 
#define fd(i, s, e, dec) for(ll i = s; i >= e; i-=dec)

#define p(x) cout << x << '\n'

#define MOD 1000000007
#define otherMOD 998244353

#define L 20

// Segtree
struct segtree{
    typedef ll T; //replace this

    T id = 0, t[2 * N]; 
    // id is identity: 1 for mult, 0 for add, ll max for min, ll min for max
    // Make sure N is >= siz of array
    
    T func(T a, T b) { return a+b; } // replace this
 
    void modify(ll p, T v) { //set value v at position p
        for(t[p += N] = v; p /= 2;) t[p] = func(t[2 * p], t[2 * p + 1]);
    }
 
    T query(ll l, ll r) { //query on interval [l, r)
        T resl = id, resr = id;
        for(l += N, r += N; l < r; l /= 2, r /= 2) {
            if(l & 1) resl = func(resl, t[l++]);
            if(r & 1) resr = func(t[--r], resr);
        }
        return func(resl, resr);
    }
};

class segtree{
    public:
        typedef ll T; //replace this

        T id; vector<T> t;

        segtree(){
            id = 0; //replace this
            t = vector<T>(2*N, 0);
        }

        segtree(ll siz){
            id = 0; //replace this
            t = vector<T>(2*siz, 0);
        }

        // id is identity: 1 for mult, 0 for add, ll max for min, ll min for max
        // Make sure N is >= siz of array
        
        T func(T a, T b) { return a+b; } // replace this
    
        void modify(ll p, T v) { //set value v at position p
            for(t[p += N] = v; p /= 2;) t[p] = func(t[2 * p], t[2 * p + 1]);
        }
    
        T query(ll l, ll r) { //query on interval [l, r)
            T resl = id, resr = id;
            for(l += N, r += N; l < r; l /= 2, r /= 2) {
                if(l & 1) resl = func(resl, t[l++]);
                if(r & 1) resr = func(t[--r], resr);
            }
            return func(resl, resr);
        }
};

class lztree{ 
    public:

    typedef ll T; // Type value
    typedef ll U; // Update value

    T idT = 0, t[2 * N];
    U idU = 0, d[N];
    ll x = (fill_n(d, N, idU), 0);

    T combine(T a, T b) { return a + b; }
    U combineUpdates(U b, U a) { return a + b; }
    T update(U b, T a) { return a + b; }

    /*
        Updating a combined segment == Combining updated segments
        update(x, combine(a, b)) = combine(update(x, a), update(x, b))

        Applying combined updates == Applying each update separately
        update(combineUpdates(x, y), a) = update(x, update(y, a))
    */

    void calc(ll p) { t[p] = update(d[p], combine(t[p * 2], t[p * 2 + 1])); }

    void apply(ll p, U v) {
        t[p] = update(v, t[p]);
        if(p < N) d[p] = combineUpdates(v, d[p]);
    }

    void push(ll p) {
        p += N;
        for(ll s = L - 1; s > (0); --s){
            ll i = p >> s;
            if(d[i] != idU) {
                apply(i * 2, d[i]);
                apply(i * 2 + 1, d[i]);
                d[i] = idU;
            }
        }
    }

    void modify(ll p, T v) {
        push(p);
        t[p += N] = v;
        while(p > 1) calc(p /= 2);
    }

    void modify(ll l, ll r, U v) {
        push(l), push(r - 1);
        bool cl = false, cr = false;
        for(l += N, r += N; l < r; l /= 2, r /= 2) {
            if(cl) calc(l - 1);
            if(cr) calc(r);
            if(l & 1) apply(l++, v), cl = true;
            if(r & 1) apply(--r, v), cr = true;
        }
        for(--l; r; l /= 2, r /= 2) {
            if(cl) calc(l);
            if(cr) calc(r);
        }
    }

    T query(ll l, ll r) {
        push(l), push(r - 1);
        T resl = idT, resr = idT;
        for(l += N, r += N; l < r; l /= 2, r /= 2) {
            if(l & 1) resl = combine(resl, t[l++]);
            if(r & 1) resr = combine(t[--r], resr);
        }
        return combine(resl, resr);
    }
};

//Bsearch
bool works()
{
    return false;
}

//Binary Search, goodbad is true if segment is [GOOD, BAD], false otherwise
ll bs(ll l, ll r, bool goodbad) 
{
  if(goodbad){
    // r should be maxR + 1
    while(r-l>1)
    {
      ll m = (l+r)/2;
      if(works()) l = m;
      else r = m;
    }
    return l;
    }
  else{
    // l should be minL - 1
    while(r-l>1) 
    {
      ll m = (l+r)/2;
      if(works()) r = m;
      else l = m;
    }
    return r;
  }
}
 


//Combinatorics

ll inv(ll a, ll b){return 1<a?b - inv(b%a,a)*b/a : 1;}
 
vector<ll> fact(ll n) // MOD MOD
{
  vector<ll> ans = {1, 1};
  f(i, 2, n+1, 1) ans.push_back((ans.back() * i)%MOD);
  return ans;
}
 
ll nck(vector<ll>& v, ll n, ll k) //MOD
{
  return v[n] * inv(v[n-k], MOD) % MOD * inv(v[k], MOD) % MOD;
}


//Number Theory

ll binpow(ll a, ll b) 
{ 
	ll res = 1;
	while (b > 0) {
		if (b & 1) res *= a;
		a *= a;
		res %= MOD;
		a %= MOD;
		b >>= 1;
	}
	return res;
}

ll gcd(ll a, ll b)
{ if(b==0)return a; else return gcd(b,a%b); }
 
ll lcm(ll a, ll b)
{ return a*b/gcd(a,b); }

map<ll,ll> primeFac(ll x)
{
  map<ll,ll> mp; 
  if(x==0 || x==1) return mp;
  else
  {
    while(x%2==0)
    {
      x/=2; mp[2]++;
    }
    for(ll i = 3; i <= sqrt(x); i = i+2) 
    { 
      while (x%i==0) 
      {
        mp[i]++;
        x = x/i; 
      } 
    } 
    if(x>2) mp[x]++;
    return mp;
  }
}
