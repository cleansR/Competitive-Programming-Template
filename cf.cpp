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
 
void dfs(vector<ll> adj[], ll marked[], ll chosen[], ll curr, ll parent)
{
  marked[curr] = true;
  if(chosen[parent]==false)
  {
    chosen[curr] = true;
  }
  else 
  chosen[curr] = false;
  for(ll i = 0; i <adj[curr].size(); i++)
  {
    if(!marked[adj[curr][i]])
    {
      dfs(adj, marked, chosen, adj[curr][i], curr);
    }
  }
}
 
bool isPrime(ll n)
{
  for(ll i = 2; i<= (ll)sqrt(n); i++)
  {
    if(n%i==0)
    {
      return false;
    }
  }
  return true;
}
 
vector<ll> toBinary(ll num) // first index of vector is smallest bit
{
  vector<ll> ans;
  ll p = 0;
  while(pow(2, p)<num)
    p++;
  if(pow(2, p)!=num)
    p--;
  while(p>=0)
  {
    if(num>=pow(2, p))
    {
      num -= (ll)pow(2, p);
      ans.push_back(1);
    }
    else
    ans.push_back(0);
    p--;
  }
 
  ll left = 0; ll right = ans.size()-1;
  while(left<right)
  {
    ll temp = ans[left];
    ans[left] = ans[right];
    ans[right] = temp; 
    left++;
    right--;
  }
  return ans;  
}
 
bool isPalindrome(string s)
{
  ll l = 0; ll r = s.length(); r--;
  while(l<r)
  {
    if(s[l]==s[r])
    {
      l++; r--;
    }
    else return false;
  }
  return true;
}
 
bool isSubstring(string mainS, string subS)
{
  if (mainS.find(subS) != string::npos) 
  {
    return true;
    } 
  return false;
 
}
 
bool works()
{
 return false;
}
 
void readArr(ll arr[], ll size)
{
  for(ll i = 0; i < size; i++)
  {
    cin >> arr[i];
  }
}
 
vector<ll> prefixSumA(ll arr[], ll size)
{
  vector<ll> pF;
  pF.push_back(0);
  for(ll i = 1; i <=size; i++)
  {
    pF.push_back(pF[i-1]+arr[i-1]);
  }
  return pF;
}
 
vector<ll> prefixSumV(vector<ll>& v)
{
  vector<ll> pF;
  pF.push_back(0);
  for(ll i = 1; i <=v.size(); i++)
  {
    pF.push_back(pF[i-1]+v[i-1]);
  }
  return pF;
}
 
ll binToLL(ll arr[], ll size)
{
  ll ans = 0;
  for(ll i = 0; i < 32; i++)
  {
    if(arr[i]==1)
    {
      ans+= (ll) pow(2, i);
    }
  }
  return ans;
}
 
void treeDFS(vector<ll> adj[], ll curr, ll parent)
{
  vector<ll> temp = adj[curr];
  for(ll i = 0; i < temp.size(); i++)
  {
    if(temp[i]!=parent)
    treeDFS(adj, temp[i], curr);
  }
}
 
void printArray(ll arr[], ll size)
{
  forI(0, size, 1){
    cout << arr[i] << " ";
  }
  cout << '\n';
}
 
void printVector(vector<ll>& v)
{
  forI(0, v.size(), 1){
    cout << v[i] << " ";
  }
  cout << '\n';
}
 
void printSet(set<ll>& s)
{
  for(auto ptr = s.begin(); ptr!=s.end(); ptr++)
  {
    cout << (*ptr) << " ";
  }
  cout << '\n';
}
 
void printMap(map<ll,ll>& mp)
{
  for (const auto &item: mp)
  {
    cout << item.first << " " << item.second << '\n';
  }
  cout << '\n';
}
 
void tf(bool b)
{
  if(b) cout << "YES" << '\n';
  else cout << "NO" << '\n';
}
 
bool sortbyCond(const pair<pair<ll,ll>, vector<ll>>& a, const pair<pair<ll,ll>, vector<ll>>& b)
{
    if (a.first.first != b.first.first)
        return (a.first < b.first);
    else
        return (a.first.second > b.first.second);     
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
 
ll gcd(ll a, ll b)
{ if(b==0)return a; else return gcd(b,a%b); }
 
ll lcm(ll a, ll b)
{ return a*b/gcd(a,b); }
 
bool validparentheses(string s)
{
  stack<char> st; ll n = s.length();
  forI(0, n, 1)
  {
    if(s[i]=='(') st.push(s[i]);
    else
    {
      if(st.size()==0)
        return false;
      else
        st.pop();
    }
  }
    if(st.size()==0)
    return true; else return false;
}
 
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
 
void deez_nuts()
{
  
}

int main()
{
                    ios::sync_with_stdio(false);
                cin.tie(NULL);
            ll sus = 1;
        cin >> sus;
  for(ll amogus = 0; amogus < sus; amogus++)
  { 
    deez_nuts();
  }
}