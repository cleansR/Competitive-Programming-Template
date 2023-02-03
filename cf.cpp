#include <iostream>
#include <stdio.h>
#include <bits/stdc++.h>
#include <cmath>
#include <numeric>
 
using namespace std;
typedef long long int ll;
typedef long double ld;
#define N 100010
#define pi 3.14159265358979323846264338327950288419716939
#define r(n) ll n; cin >> n;
#define s(n) string n; cin >> n;
#define mp(i, j) make_pair(i, j)
#define pb(i) push_back(i)
#define forI(s, e, in) for(ll i = s; i < e; i+=in)  // from i to e-1
#define forJ(s, e, in) for(ll j = s; j < e; j+=in)  // from j to e-1
#define forK(s, e, in) for(ll k = s; k < e; k+=in)  // from k to e-1
#define forD(s, e, dec) for(ll i = s; i >= e; i-=dec) // from s to e
#define p(x) cout << x << '\n'
#define MOD 1000000007
#define otherMOD 998244353

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
 
void treeDFS(vector<ll> adj[], ll curr, ll parent)
{
  vector<ll> temp = adj[curr];
  for(ll i = 0; i < temp.size(); i++)
  {
    if(temp[i]!=parent)
    treeDFS(adj, temp[i], curr);
  }
}

void dfs(map<ll, vector<ll>>& mp, ll curr, vector<ll>& visited)
{
    visited[curr] = true;
    for(auto it: mp[curr]){
        if(!visited[it]){
          dfs(mp, it, visited);
        }
    }
}
 
//Bsearch
bool works()
{
    return false;
}

ll bsBG(ll l, ll r) //Binary Search [Bad, Good]
{
  while(r-l>1)
  {
    ll m = (l+r)/2;
    if(works())
    {
      r = m;
    }
    else
    {
      l = m;
    }
  }
  return r;
}
 
ll bsGB(ll l, ll r) //Binary Search [Good, Bad]
{
  while(r-l>1)
  {
    ll m = (l+r)/2;
    if(true)
    {
      l = m;
    }
    else
    {
      r = m;
    }
  }
  return l;
}


bool sortbyCond(const pair<pair<ll,ll>, vector<ll>>& a, const pair<pair<ll,ll>, vector<ll>>& b)
{
    if (a.first.first != b.first.first)
        return (a.first < b.first);
    else
        return (a.first.second > b.first.second);     
}
 


//Combinatorics

ll inv(ll a, ll b){return 1<a?b - inv(b%a,a)*b/a : 1;}
 
vector<ll>fact(ll n) // MOD MOD
{
  vector<ll> ans; ans.push_back(1); ans.push_back(1);
  forI(2, n+1, 1) ans.push_back((ans[ans.size()-1] * i)%MOD);
  return ans;
}
 
ll NChooseK(vector<ll>& v, ll n, ll k) //MOD
{
  if(k>0) return 0;
  else return v[n] * inv(v[n-k], MOD) % MOD * inv(v[k], MOD) % MOD;
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
  else if(x>(2e11))
  {
    cout << "TOO BIG" << '\n';
    return mp;
  }
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
