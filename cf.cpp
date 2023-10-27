#include <bits/stdc++.h>
using namespace std;

typedef long long int ll;
typedef long double ld;

#include <ext/pb_ds/assoc_container.hpp>
#include <ext/pb_ds/tree_policy.hpp>
using namespace __gnu_pbds;
  
#define ordered_set tree<int, null_type,less<int>, rb_tree_tag,tree_order_statistics_node_update>

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


long long power(long long x, long long y, long long p) {
    long long result = 1;
    for(; y; y >>= 1, x = x * x % p) {
        if(y & 1) {
            result = result * x % p;
        }
    }
    return result;
}
 
long long inverse(long long x, long long p) {
    return power(x, p - 2, p);
}
 
struct Hash {
    int length;
    const int mod1 = 1e9 + 7, mod2 = 1e9 + 9;
    const int p1 = 31, p2 = 37;
    vector<int> hash1, hash2;
    pair<int, int> hash_pair;
 

    inline static vector<int> inv_pow1, inv_pow2;
    inline static int inv_size = 1;
     
    Hash() {}
 
    Hash(const string& s) {
        length = s.size();
        hash1.resize(length);
        hash2.resize(length);
 
        int h1 = 0, h2 = 0;
        long long p_pow1 = 1, p_pow2 = 1;
        for(int i = 0; i < length; i++) {
            h1 = (h1 + (s[i] - 'a' + 1) * p_pow1) % mod1;
            h2 = (h2 + (s[i] - 'a' + 1) * p_pow2) % mod2;
            p_pow1 = (p_pow1 * p1) % mod1;
            p_pow2 = (p_pow2 * p2) % mod2;
            hash1[i] = h1;
            hash2[i] = h2;
        }
        hash_pair = make_pair(h1, h2);
 
        if(inv_size < length) {
            for(; inv_size < length; inv_size <<= 1);
             
            inv_pow1.resize(inv_size, -1);
            inv_pow2.resize(inv_size, -1);
 
            inv_pow1[inv_size - 1] = inverse(power(p1, inv_size - 1, mod1), mod1);
            inv_pow2[inv_size - 1] = inverse(power(p2, inv_size - 1, mod2), mod2);
             
            for(int i = inv_size - 2; i >= 0 && inv_pow1[i] == -1; i--) {
                inv_pow1[i] = (1LL * inv_pow1[i + 1] * p1) % mod1;
                inv_pow2[i] = (1LL * inv_pow2[i + 1] * p2) % mod2;
            }
        }
    }
 
    int size() {
        return length;
    }
 
    pair<int, int> prefix(const int index) {
        return {hash1[index], hash2[index]};
    }
 
    pair<int, int> substr(const int l, const int r) {
        if(l == 0) {
            return {hash1[r], hash2[r]};
        }
        int temp1 = hash1[r] - hash1[l - 1];
        int temp2 = hash2[r] - hash2[l - 1];
        temp1 += (temp1 < 0 ? mod1 : 0);
        temp2 += (temp2 < 0 ? mod2 : 0);
        temp1 = (temp1 * 1LL * inv_pow1[l]) % mod1;
        temp2 = (temp2 * 1LL * inv_pow2[l]) % mod2;
        return {temp1, temp2};
    }
 
    bool operator==(const Hash& other) {
        return (hash_pair == other.hash_pair);
    }
};











// trie

struct arr {
    int a[26];
    arr() {}
    void clear() { memset(a,-1,sizeof(a)); }
    int& operator[](int i) { return a[i]; }
};

// N is the max # of letters in a word
struct trie {
    int cnt, prefix_cnt[N], word_cnt[N];
    arr to[N];

    trie() { cnt = N-1; }

    void clear() { for(int i = 0; i < cnt; i++) prefix_cnt[i] = word_cnt[i] = 0, to[i].clear(); cnt = 1; }

    void add(const string& s) {
        int u = 0;
        for(const char& c: s)  {
            if(to[u][c-'a'] == -1) to[u][c-'a'] = cnt++;
            u = to[u][c-'a'];
            prefix_cnt[u]++;
        }
        word_cnt[u]++;
    }

    int prefix_count(const string& s) {
        int u = 0;
        for (const char& c: s) {
            if (to[u][c-'a'] == -1) return 0;
            u = to[u][c-'a'];
        }
        return prefix_cnt[u];
    }
} tr;
 


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

struct lztree{ 
    
    typedef ll T; // Type value
    typedef ll U; // Update value

    T idT = 0, t[2 * N];
    U idU = 0, d[N];
    ll x = (fill_n(d, N, idU), 0);

    T combine(T a, T b) { return max(a, b); }
    U combineUpdates(U b, U a) { a + b; }
    T update(U b, T a) { return a + b; }

    /*
        Updating a combined segment == Combining updated segments
        update(x, combine(a, b)) = combine(update(x, a), update(x, b))

        Applying combined updates == Applying each update separately
        update(combineUpdates(x, y), a) = update(x, update(y, a))

        Usually u change the combine method for querying
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

// DSU

ll id[N], sz[N];

void setUp(ll siz)
{
    f(i, 0, siz, 1){
        id[i] = i; sz[i] = 1;
    }
}

ll find(ll v) {
    if (v == id[v])
        return v;
    return id[v] = find(id[v]);
}

void merge(ll u, ll v) {
    ll a = find(u);
    ll b = find(v);
    if (a != b) {
        if (sz[a] < sz[b])
            swap(a, b);
        id[b] = a;
        sz[a] += sz[b];
    }
}


//Bsearch
bool works()
{
    return false;
}

// [good bad]
ll bs(ll l, ll r) 
{
    // r should be maxR + 1
    while(r-l>1)
    {
        ll m = (l+r)/2;
        if(works()) l = m;
        else r = m;
    }
    return l;
}

// [bad good]
ll bs(ll l, ll r)
{
    // l should be minL - 1
    while(r-l>1) 
    {
        ll m = (l+r)/2;
        if(works()) r = m;
        else l = m;
    }
    return r;
}


// Custom comparator
struct cmp{
    // Return true if a goes before b
    bool operator() (pair<ll, ll> a, pair<ll, ll> b) const{
        if(a.first < b.first) return true;
        else if(a.first > b.first) return false;
        else{
            return a.second > b.second;
        }
    }
};


//Combinatorics

ll inv(ll a, ll b)
{
    return 1<a?b - inv(b%a,a)*b/a : 1;
}
 
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
{ 
    if(b==0) return a; else return gcd(b,a%b); 
}
 
ll lcm(ll a, ll b)
{ 
    return a*b/gcd(a,b); 
}

void seive()
{
    int n = 1e4; 
    bool prime[n+1];
    memset(prime, true, sizeof(prime));
  
    for (int p = 2; p * p <= n; p++) {
        if (prime[p]) {

            for (int i = p * p; i <= n; i += p)
                prime[i] = false;
        }
    }
    prime[0] = false; prime[1] = false;
}

map<ll,ll> primeFac(ll x)
{
    map<ll,ll> mp; 
    if(x==0 || x==1) return mp;
    else{
        while(x%2==0){
            x/=2; mp[2]++;
        }
        for(ll i = 3; i <= sqrt(x); i = i+2) { 
            while (x%i==0) {
                mp[i]++; x = x/i; 
            } 
        } 
        if(x>2) mp[x]++;
        return mp;
    }
}


// 2D prefix sums

// The parameters are the subrectangle of the original array
ll access2dPFS(vector<vector<ll>>& pfs, ll sr, ll br, ll sc, ll bc)
{
    return pfs[br+1][bc+1] - pfs[br+1][sc] - pfs[sr][bc+1] + pfs[sr][sc];
}

vector<vector<ll>> create2dPFS(vector<vector<ll>>& v)
{
    ll n = v.size(); ll m = v[0].size();
    vector<vector<ll>> pfs(n+1, vector<ll>(m+1, 0)); pfs[1][1] = v[0][0];
    f(i, 1, n+1, 1){
        f(j, 1, m+1, 1){
            pfs[i][j] = v[i-1][j-1] + pfs[i-1][j] + pfs[i][j-1] - pfs[i-1][j-1];
        }
    }
}


// Graph

void dfs(vector<ll> adj[], ll curr, vector<bool>& visited)
{
    visited[curr] = true; 
    for(auto it: adj[curr]){
        if(!visited[it]) dfs(adj, it, visited);
    }
}

void bfs(vector<ll> adj[], ll src, vector<ll>& dist)
{
    fill(dist.begin(), dist.end(), -1);
    queue<ll> q; q.push(src); dist[src] = 0;
    while(q.size() > 0){
        ll v = q.front(); q.pop();
        for(auto it: adj[v]){
            if(dist[it] == -1){
                q.push(it); dist[it] = dist[v] + 1;
            }
        }
    }
}

/*
    Returns a pair {{end1, end2}, dist}
    Assumes adj is 1 indexed (should work if 0 indexed also)
*/
pair<pair<ll, ll>, ll> diameter(vector<ll> adj[], ll n)
{
    ll end1, end2, dist;
    ll start = 1; vector<ll> distFromStart(n+1, -1);
    bfs(adj, start, distFromStart);

    ll mx = 0; ll furthestNode;
    f(i, 0, distFromStart.size(), 1){
        if(distFromStart[i] > mx){
            mx = distFromStart[i];
            furthestNode = i;
        }
    }
    
    end1 = furthestNode;

    vector<ll> distFromFurthest(n+1, -1);
    bfs(adj, furthestNode, distFromFurthest);
    mx = 0; 
    f(i, 0, distFromFurthest.size(), 1){
        if(distFromFurthest[i] > mx){
            mx = distFromFurthest[i];
            furthestNode = i;
        }
    }

    end2 = furthestNode;

    dist = mx;

    return {{end1, end2}, dist};
}


struct LCA {
    vector<ll> height, euler, first, segtree;
    vector<bool> visited;
    ll n;

    LCA(vector<vector<ll>> &adj, ll root = 0) {
        n = adj.size();
        height.resize(n);
        first.resize(n);
        euler.reserve(n * 2);
        visited.assign(n, false);
        dfs(adj, root);
        ll m = euler.size();
        segtree.resize(m * 4);
        build(1, 0, m - 1);
    }

    void dfs(vector<vector<ll>> &adj, ll node, ll h = 0) {
        visited[node] = true;
        height[node] = h;
        first[node] = euler.size();
        euler.push_back(node);
        for (auto to : adj[node]) {
            if (!visited[to]) {
                dfs(adj, to, h + 1);
                euler.push_back(node);
            }
        }
    }

    void build(ll node, ll b, ll e) {
        if (b == e) {
            segtree[node] = euler[b];
        } else {
            ll mid = (b + e) / 2;
            build(node << 1, b, mid);
            build(node << 1 | 1, mid + 1, e);
            ll l = segtree[node << 1], r = segtree[node << 1 | 1];
            segtree[node] = (height[l] < height[r]) ? l : r;
        }
    }

    ll query(ll node, ll b, ll e, ll l, ll R) {
        if (b > R || e < l)
            return -1;
        if (b >= l && e <= R)
            return segtree[node];
        ll mid = (b + e) >> 1;

        ll left = query(node << 1, b, mid, l, R);
        ll right = query(node << 1 | 1, mid + 1, e, l, R);
        if (left == -1) return right;
        if (right == -1) return left;
        return height[left] < height[right] ? left : right;
    }

    ll lca(ll u, ll v) {
        ll left = first[u], right = first[v];
        if (left > right)
            swap(left, right);
        return query(1, 0, euler.size() - 1, left, right);
    }
};

int main() {
    // reading in entire lines
    string s;
    getline(cin, s);
    cout << s << '\n';
    cout << fixed << setprecision(15) << 15 << '\n';
}
