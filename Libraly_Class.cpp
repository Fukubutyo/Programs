#include<iostream>
#include<algorithm>
#include<string>
#include<cstdlib>
#include<map>
#include<iomanip>
#include<sstream>
#include<vector>
#include<stack>
#include<math.h>
#include<queue>
#include<complex>
#include<random>
#include<ctime>
using namespace std;

//ユークリッドの互除法　a,bは最大公約数を求めたい２つの数
class Euclid_Gojyohou{
    public:
    long long int gcd(long long int a, long long int b) {
        long long int tmp;
        long long int r = 1;
        if (b > a) {
            tmp = a;
            a = b;
            b = tmp;
        }
        r = a % b;
        while (r != 0) {


            a = b;
            b = r;
            r = a % b;

        }
        return b;
    }
};

class Extend_Euclid_Gojyohou{
    public:
    long long int memo[2000]={},GCD;
    long long int keisuu[2000]={};
    long long int gcd(long long int a, long long int b) {
        memo[0]=a;
        keisuu[0]=1;
        memo[1]=b;
        
        long long int tmp,count=2,p,q;
        long long int r = 1;
        if (b > a) {
            tmp = a;
            a = b;
            b = tmp;
        }
        r = a % b;
        memo[2]=r;
        keisuu[1]=-(a/b);
        while (r != 0) {


            a = b;
            b = r;
            r = a % b;
                
            memo[count+1]=r;
            keisuu[count]=-(a/b);
            count++;
        }
        GCD=b;
        count-=2;
        p=1;
        q=keisuu[count];
        if(q>=0)cout<<p<<"*"<<memo[count-1]<<"+"<<q<<"*"<<memo[count]<<"="<<GCD<<endl;
        else{cout<<p<<"*"<<memo[count-1]<<q<<"*"<<memo[count]<<"="<<GCD<<endl;}
        for(int i=count-1;i>=1;i--){
            long long int tmptmp=p;
            p=q;
            q=tmptmp+q*keisuu[i];
            if(q>0)cout<<p<<"*"<<memo[i-1]<<"+"<<q<<"*"<<memo[i]<<"="<<GCD<<endl;
            else{cout<<p<<"*"<<memo[i-1]<<q<<"*"<<memo[i]<<"="<<GCD<<endl;}
        }
        //cout<<p<<"*"<<memo[0]<<"+"<<q<<"*"<<memo[1]<<"="<<GCD;
        return b;
    }
};

//aからbまでの積を求める(mod)

//階乗先に計算
class Kaizyou{
    public:
    long long int *kaizyou;
    Kaizyou(){};
    Kaizyou(long long int n,long long int mod){
        kaizyou=new long long int[n+1];
        kai(1,n,mod);
        return;
    }

    long long int kai(long long int a, long long int b, long long int mod) {
        long long int tmp = 1;
        for (long long int i = a; i <= b; i++) {
            tmp *= i;
            tmp %= mod;

            kaizyou[i] = tmp;
        }
        return tmp;
    }
};

//累乗(繰り返し2乗法)　aのb乗(mod)を求める。
class Ruizyou{
    public:
    long long int rui(long long int a, long long int b, long long mod) {
        int memo[65] = {};
        long long int A[65] = {};
        long long int tmp = 1;
        for (int i = 0; i < 65; i++) {
            memo[i] = b % 2;
            b /= 2;
        }

        A[0] = a;
        A[0] %= mod;

        for (int i = 1; i < 65; i++) {
            A[i] = A[i - 1] * A[i - 1];
            A[i] %= mod;
        }
        for (int i = 0; i < 65; i++) {
            if (memo[i] == 1) {
                tmp *= A[i];
                tmp %= mod;
            }
        }
        tmp %= mod;
        return tmp;
    }
};
//コンビネーション計算
class Combination{
    public:
        Kaizyou * Kaizyou_Instance;
        Ruizyou Ruizyou_Instance;
    Combination(long long int n,long long int mod){
        Kaizyou_Instance=new Kaizyou(n,mod);
    }
    long long int comb(long long int n, long long int r, long long int mod) {
        long long int tmp;

        tmp = (Kaizyou_Instance->kaizyou[n] * Ruizyou_Instance.rui(Kaizyou_Instance->kaizyou[r], mod - 2, mod)) % mod;
        tmp *= Ruizyou_Instance.rui(Kaizyou_Instance->kaizyou[n - r], mod - 2, mod);
        tmp %= mod;
        if (tmp < 0) { tmp = (mod - tmp) % mod; }
        return tmp;
    }
};	

//素数判定　素数ならそのまま返す　合成数なら０を返す。
class Sosu{
	public:
	long long int sosu(long long int n) {
		for (int i = 2; i*i <= n; i++) {
			if (n%i == 0) { return 0; }
		}
		return n;

	}
};

//union find木(グループわけができる)
class Union_Find{
	public:
	Union_Find(int MAX){
		par =new long long int[MAX+1];
		
		for (int i = 0; i < MAX; i++) {
			par[i] = i;
		}
	}
	long long int * par;
		//木の根を求める
	int unifind_root(int x) {
		if (par[x] == x) {
			return x;
		}
		else {
			return par[x] = unifind_root(par[x]);
		}
	}
		//同じグループか判定
	bool unifind_same(int x, int y) {
		return unifind_root(x) == unifind_root(y);
	}
		//グループ併合
	int unifind_unite(int x, int y) {
		x = unifind_root(x);
		y = unifind_root(y);
		if (x == y) {
			return 0;
		}
		par[x] = y;
	}
};

//CRT(中国式剰余定理)
class ChineseRemainderTheorem{
	public:
long long int Euclid_CRT(long long int p, long long int q, long long int *x, long long int *y) {
		long long int tmp, r, a1 = 1, b1 = 0, a2 = 0, b2 = 1, amemo, bmemo;
		r = p % q;

		while (r != 0) {

			amemo = a1 - (p / q)*a2;
			bmemo = b1 - (p / q)*b2;
			a1 = a2;
			b1 = b2;
			a2 = amemo;
			b2 = bmemo;
			p = q;
			q = r;
			r = p % q;
		}

		*x = a2;
		*y = b2;
		return q;
	}
	//pで割ってpr余りqで割ってqr余るような数
	long long int CRT(long long int p, long long int q, long long int rp, long long int rq, long long int *gcd) {
		long long int x, y, tmp, GCD, psix, ten[15] = { 1,10,100,1000,10000,100000,1000000,10000000,100000000,1000000000,10000000000,100000000000,1000000000000,10000000000000,100000000000000 };
		if (q > p) {
			tmp = p;
			p = q;
			q = tmp;
			tmp = rp;
			rp = rq;
			rq = tmp;
		}
		rp %= p;
		rq %= q;

		GCD = Euclid_CRT(p, q, &x, &y);
		*gcd = GCD;
		long long int lastp, lastr = 0;
		lastp = p * q / GCD;
		x *= ((rq - rp) / GCD);
		y *= ((rq - rp) / GCD);

		for (int i = 0; i <= 12; i += 6) {
			lastr += ((((p / ten[i]) % 1000000)*x) % lastp)*ten[i];

		}

		lastr += rp;
		lastr %= (lastp);
		if (lastr < 0) { lastr = lastp + lastr; }
		return lastr;
	}
};

//約数列挙 累乗の関数とともに用いる。
class yakusuu_rekkyo{
	public:
	long long int soinnsuu[2100000] = {}, yakusuu[2100000] = {1}, kosuu = 1;
	yakusuu_rekkyo(long long int n,long long int mod) {
		Ruizyou Ruizyou_Instance;
		long long int N = n;
		for (int i = 2; i * i <= N; i++) {
			int tmp = 1;
			while (n % i == 0) {
				soinnsuu[i]++;
				n /= i;
				tmp++;
			}
			int tmp2 = tmp;

			while (tmp!=1) {
				for (int j = 0; j < kosuu; j++) {
					yakusuu[j + (tmp - 1) * kosuu] = Ruizyou_Instance.rui(i, tmp - 1, mod)*yakusuu[j];
				}
				tmp--;
			}
			kosuu *= tmp2;

		}
		if (n != 1) {
			for (int j = 0; j < kosuu; j++) {
				yakusuu[j + kosuu] = Ruizyou_Instance.rui(n, 1, mod) * yakusuu[j];
			}
			kosuu *= 2;
		}
	}
};
//ベルマンフォード(ある一点からの最短経路)
class Bellman_Ford{
	public:
	const long long int INF = 5000000000000000000;

	struct edge {
		long long int from, to, cost;
	};
	edge ed[200101];
	long long int d[201000];
	long long int V, E;

	void bellman_Ford(int s) {
		for (int i = 0; i < V; i++) { d[i] = INF; }
		d[s] = 0;
		while (true) {
			bool update = false;
			for (int i = 0; i < E; i++) {
				edge e = ed[i];
				if (d[e.from] != INF && d[e.to] > d[e.from] + e.cost) {
					d[e.to] = d[e.from] + e.cost;
					update = true;
				}
			}
			if (!update) { break; }
		}
	}
};
//ダイクストラ(ある一点からの最短経路)
class Dijkstra{
	const long long int INF = 5000000000000000000;

	struct edge {
		long long int  to, cost;
	};
	typedef pair<long long int, long long int> P;
	vector<edge> G[200101];
	long long int d[201000];
	long long int V, E;

	void dijkstra(long long int s) {

		priority_queue<P, vector<P>, greater<P> >que;
		fill(d, d + V, INF);
		d[s] = 0;
		que.push(P(0, s));

		while (!que.empty()) {
			P p = que.top();
			que.pop();
			long long int v = p.second;
			if (d[v] < p.first)continue;
			for (int i = 0; i < G[v].size(); i++) {
				edge e = G[v][i];
				if (d[e.to] > d[v] + e.cost) {
					d[e.to] = d[v] + e.cost;
					que.push(P(d[e.to], e.to));
				}
			}
		}
	}
};
//ワーシャルフロイド
class Warshall_Floyd{
	public:
	long long int **d,nmemo;
	Warshall_Floyd(long long int n) {	
		nmemo=n;	
		long long int **d=new long long int*[n];
		for(int i=0;i<n;i++){
			d[i]=new long long int [n];
		}
		for (int k = 0; k < n; k++) {
			for (int i = 0; i < n; i++) {
				for (int j = 0; j < n; j++) {
					d[i][j] = min(d[i][j], d[i][k] + d[k][j]);
				}
			}
		}
	}
	void Do_Warshall_Floyd() {		
		for (int k = 0; k < nmemo; k++) {
			for (int i = 0; i < nmemo; i++) {
				for (int j = 0; j < nmemo; j++) {
					d[i][j] = min(d[i][j], d[i][k] + d[k][j]);
				}
			}
		}
	}
};


//最長単調増加列,0-indexed,要素数Nと配列a
class LIS{
    public:
	const long long int INF=5000000000000000000;
	long long int LIS_Do(long long int N,long long int a[]){
		long long int LIS_table[250005];
		long long int LIS_start=250000,now=LIS_start,ng=now,ok=now+1,mid;
		for(int i=0;i<=LIS_start;i++){
			LIS_table[i]=INF;
		}
		for(int i=0;i<N;i++){
			ng=now;
			ok=LIS_start+1;
			
			while(ok-ng>1){
				mid=(ng+ok)/2;
				if(LIS_table[mid]>a[i]){
					ng=mid;
				}else{ok=mid;}
			}
			
			LIS_table[ng]=a[i];
			if(ng==now){
				now--;
			}
			
		}
		for(int i=LIS_start;i<=now;i++){
			cout<<LIS_table[i]<<endl;
		}
		return LIS_start-now;
	}
    
};


//文字列strとstr.substr(i,str.length()-1)の先頭何文字が一致するかメモ
class Z_Algorithm{
	long long int * prefix_length;
	Z_Algorithm(string str){
		prefix_length=new long long int[str.length()];  
		prefix_length[0]=str.length();
		int L=0,R=0;
		for(int i=1;i<str.length();i++){
			if(R<=i){
				L=i;
				R=i;
				while(R<str.length()&&str[R]==str[R-L])R++;
				prefix_length[i]=R-L;
				
			}else{
				int k=i-L;
				if(prefix_length[k]<R-i){
					prefix_length[i]=prefix_length[k];
				}else{
					L=i;
					while(R<str.length()&&str[R]==str[R-L])R++;
					prefix_length[i]=R-L;
					
				}
			}
		}
	}
};

//セグメント木
class SegmentTree{
    public:
	long long int dataSize,identityElement;
	pair<long long int,long long int>* treeDataRange;
	long long int* treeData;

	//要素数nと単位元e
	SegmentTree(long long int n,long long int e){
		dataSize=2;
		while(dataSize<n){
			dataSize*=2;
		}
		identityElement=e;
		treeData=new long long int [2*dataSize+1];
		treeDataRange=new pair<long long int,long long int> [2*dataSize+1];
		for(int i=dataSize;i<2*dataSize;i++){
			treeDataRange[i]={i-dataSize,i-dataSize+1};
		}
		for(int i=dataSize-1;i>=1;i--){
			treeDataRange[i]={treeDataRange[2*i].first,treeDataRange[2*i+1].second};
		}
		for(int i=1;i<2*dataSize;i++){
			treeData[i]=e;
		}
	}

    long long int DataChange(long long int id,long long int data){
        treeData[id+dataSize]=data;
        DataChangeLoop((id+dataSize)/2);
		return 0;
    }
    long long int DataChangeLoop(long long int id){
        treeData[id]=Calculation(treeData[2*id],treeData[2*id+1]);
        if(id==1){return 0;}
        DataChangeLoop(id/2);
    }

	long long int Calculation(long long int x,long long int y){
		return x+y;
	}
	//[l,r)の結果を出力
	long long int OutPut(long long int l,long long int r){
		return OutPutLoop(1,l,r);
	}
	
	long long int OutPutLoop(long long int now,long long int l,long long int r){
		
		if(l>=r){return identityElement;}

        if(l==treeDataRange[now].first&&r==treeDataRange[now].second){return treeData[now];}

		return Calculation(OutPutLoop(2*now,l,min(treeDataRange[2*now].second,r)),OutPutLoop(2*now+1,max(treeDataRange[2*now+1].first,l),r));
		
	}
};



//0-indexedの幅優先探索
class BreadthFirstSearch
{
public:
	long long int **Used,**BFSResult;
	char **HW;
	long long int MaxH,MaxW;
	stack<pair<long long int,long long int> >Root;
	
	BreadthFirstSearch(long long int H,long long int W){
			MaxH=H;
			MaxW=W;
			HW= new char*[H+1];
			Used= new long long int*[H+1];
			BFSResult= new long long int*[H+1];
			for(int i=0;i<H;i++)
			{
					HW[i]= new char[W+1];
					Used[i]= new long long int[W+1];
					BFSResult[i]= new long long int[W+1];
			}
	}
	
	void DoBFS(pair<long long int,long long int>Start){
		for(int i=0;i<MaxH;i++)
		{
			for(int j=0;j<MaxW;j++)
			{
				Used[i][j]=0;
				BFSResult[i][j]=0;
			}
		}
		queue<pair<long long int,long long int> >task;
		task.push(Start);
		Used[Start.first][Start.second]=1;
		
		while(!task.empty())
		{
			long long int h=task.front().first;
			long long int w=task.front().second;
			for(int i=-1;i<=1;i++)
			{
				for(int j=-1;j<=1;j++)
				{
					if(h+i<0||w+j<0||h+i>=MaxH||w+j>=MaxW||i!=0&&j!=0)
					{
							continue;
					}
					if(Used[h+i][w+j]==0&&HW[h+i][w+j]=='.')
					{
							Used[h+i][w+j]=1;
							BFSResult[h+i][w+j]=BFSResult[h][w]+1;
							task.push({h+i,w+j});
					}
				}
			}
			task.pop();
		}
	}
	
	void RootCheck(pair<long long int,long long int>Goal){
		pair<long long int,long long int>now=Goal;
		while(!Root.empty()){
			Root.pop();
		}
		Root.push(Goal);
		while(BFSResult[now.first][now.second]!=0)
		{
			bool key=false;
			for(int i=-1;i<=1;i++)
			{
				for(int j=-1;j<=1;j++)
				{
					if(now.first+i<0||now.second+j<0||now.first+i>=MaxH||now.second+j>=MaxW||i!=0&&j!=0)
					{
						continue;
					}
					if(BFSResult[now.first+i][now.second+j]==BFSResult[now.first][now.second]-1&&HW[now.first+i][now.second+j]=='.')
					{
						Root.push({now.first+i,now.second+j});
						key=true;
						now={now.first+i,now.second+j};
						break;
					}
				}
				if(key)
				{
					break;
				}
			}
		}
	}
        
};


//INFを宣言すること。Edgeに辺の情報を入れよう。
class Graph{
//メンバ変数の説明
//Edge:一つ目の番号が辺が出る方向
//MaxV,MaxE:頂点数と辺の数
//Result系:各メソッドの答え
//InfinteUpdate:BellmanFordで各番号の頂点までの最短経路が更新されるか
//
public:
	vector<pair<long long int,long long int> > *Edge;
	vector<pair<pair<long long int,long long int>,long long int> >PrimEdge;
	long long int MaxV,MaxE;
	long long int *BellmanFordResult,*DijkstraResult,**WarshallFloydResult;
	bool *InfiniteUpdate;
	bool NegativeLoop=false;
	long long int PrimResult=0;
	bool NotConected=false;


	Graph(long long int InputV,long long int InputE){
		MaxV=InputV;
		MaxE=InputE;
		Edge=new vector<pair<long long int,long long int> >[InputV+1];

	}
	
	void BellmanFord(long long int Start) {
		BellmanFordResult=new long long int[MaxV];
		InfiniteUpdate=new bool[MaxV];
		for (int i = 0; i < MaxV; i++) { BellmanFordResult[i] = INF;InfiniteUpdate[i]=false; }
		BellmanFordResult[Start] = 0;
		long long int turn=0;
		while (true) {
			turn++;
			bool update = false;
			for (int i = 0; i < MaxV; i++) {
				for(int j=0;j<Edge[i].size();j++){
					pair<long long int,long long int> e = Edge[i][j];
					if(turn>=MaxV&&InfiniteUpdate[i])
					{
						InfiniteUpdate[e.first]=true;
					}
					if (BellmanFordResult[i] != INF && BellmanFordResult[e.first] > BellmanFordResult[i] + e.second) {
						BellmanFordResult[e.first] = BellmanFordResult[i] + e.second;
						update = true;
						if(turn>=MaxV){
							InfiniteUpdate[e.first]=true;
						}
					}
				}
			}
			if (!update) { 
				NegativeLoop=false;
				break; 
			}
			if(turn>=2*MaxV-1){
				NegativeLoop=true;
				break;
			}
		}
	}

	void Dijkstra(long long int Start) {
		DijkstraResult=new long long int[MaxV];
		priority_queue<pair<long long int,long long int>, vector<pair<long long int,long long int> >, greater<pair<long long int,long long int> > >que;
		fill(DijkstraResult, DijkstraResult + MaxV, INF);
		DijkstraResult[Start] = 0;
		que.push(pair<long long int,long long int>(0, Start));

		while (!que.empty()) {
			pair<long long int,long long int> p = que.top();
			que.pop();
			long long int v = p.second;
			if (DijkstraResult[v] < p.first)continue;
			for (int i = 0; i < Edge[v].size(); i++) {
				pair<long long int,long long int> e = Edge[v][i];
				if (DijkstraResult[e.first] > DijkstraResult[v] + e.second) {
					DijkstraResult[e.first] = DijkstraResult[v] + e.second;
					que.push(pair<long long int,long long int>(DijkstraResult[e.first], e.first));
				}
			}
		}
	}

	void WarshallFloyd() {
		WarshallFloydResult=new long long int*[MaxV];
		for(int i=0;i<MaxV;i++){
			WarshallFloydResult[i]=new long long int[MaxV];
		}
		for(int i=0;i<MaxV;i++){
			for(int j=0;j<Edge[i].size();j++){
				WarshallFloydResult[i][Edge[i][j].first]=Edge[i][j].second;
			}
		}

		for (int k = 0; k < MaxV; k++) {
			for (int i = 0; i < MaxV; i++) {
				for (int j = 0; j < MaxV; j++) {
					WarshallFloydResult[i][j] = min(WarshallFloydResult[i][j], WarshallFloydResult[i][k] + WarshallFloydResult[k][j]);
				}
			}
		}
	}
	void Prim(){
		NotConected=false;
		priority_queue<pair<long long int,pair<long long int,long long int> >,vector<pair<long long int,pair<long long int,long long int> > > ,greater<pair<long long int,pair<long long int,long long int> > > >PriorityEdge;
		bool *IncludedVertex=new bool[MaxV];
		fill(IncludedVertex,IncludedVertex+MaxV,false);
		PrimEdge.clear();
		for(int i=0;i<Edge[0].size();i++)
		{
			PriorityEdge.push({Edge[0][i].second,{0,Edge[0][i].first}});
		}
		IncludedVertex[0]=true;
		for(int i=0;i<MaxV-1;i++)
		{
			while(1)
			{
				if(PriorityEdge.empty())
				{
					NotConected=true;
					break;
				}
				pair<long long int,pair<long long int,long long int> > pe = PriorityEdge.top();
				PriorityEdge.pop();
				//cout<<pe.first<<" "<<pe.second.first<<" "<<pe.second.second<<endl;
				if(!IncludedVertex[pe.second.second])
				{
					IncludedVertex[pe.second.second]=true;
					PrimResult+=pe.first;
					PrimEdge.push_back({pe.second,pe.first});
					//cout<<"size"<<pe.second.second<<" "<<Edge[pe.second.second].size()<<endl;
					for(int j=0;j<Edge[pe.second.second].size();j++)
					{
						PriorityEdge.push({Edge[pe.second.second][j].second,{pe.second.second,Edge[pe.second.second][j].first}});
					}
					break;
				}
					
					
			}
		}
	}
};
