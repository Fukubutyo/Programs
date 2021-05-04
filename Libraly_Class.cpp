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
    long long int Memo[200]={},GCD,Now=1;
    long long int Keisuu[200]={};
    bool SwapMemo=false;
    pair<long long int,long long int> RunEEG(long long int a, long long int b) {
        if(a<b){
            SwapMemo=true;
            swap(a,b);
        }
        if(a%b==0){
            if(SwapMemo){
                return{1,0};
            }else{
                return{0,1};
            }
        }
        for(long long int i = 0; i < 100; i++){
            Memo[i]=0;
            Keisuu[i]=0;
        }
        Memo[0]=a;
        Memo[1]=b;
        Now=1;
        
        while (Memo[Now-1]%Memo[Now]!=0)
        {
            Now++;
            Memo[Now]=Memo[Now-2]%Memo[Now-1];
        }
        Now--;

        Keisuu[Now-1]=1;
        Keisuu[Now]=-Memo[Now-1]/Memo[Now];
        
        while (Now>=2)
        {
            Now--;
            Keisuu[Now]+=(-Keisuu[Now+1]*(Memo[Now-1]/Memo[Now]));
            Keisuu[Now-1]=Keisuu[Now+1];
            //cout<<Now<<endl;
        }

        if(SwapMemo){
            return{Keisuu[1],Keisuu[0]};
        }else{
            return{Keisuu[0],Keisuu[1]};
        }
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
		kaizyou[0]=1;
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
        long long int r=1;
        while(b>0){
            if(1&b){
                r=(r*a)%mod;
            }
            a=(a*a)%mod;
            b=(b>>1);
        }
        return r;
    }
};
//コンビネーション計算
class Combination{
    public:
        Kaizyou * Kaizyou_Instance;
        Ruizyou Ruizyou_Instance;
		long long int CombinationMod;
    Combination(long long int n,long long int mod){
        Kaizyou_Instance=new Kaizyou(n,mod);
		CombinationMod=mod;
    }
    long long int comb(long long int n, long long int r) {
        long long int tmp;
		if(r>n){return 0;}
        tmp = (Kaizyou_Instance->kaizyou[n] * Ruizyou_Instance.rui(Kaizyou_Instance->kaizyou[r], CombinationMod - 2, CombinationMod)) % CombinationMod;
        tmp *= Ruizyou_Instance.rui(Kaizyou_Instance->kaizyou[n - r], CombinationMod - 2, CombinationMod);
        tmp %= CombinationMod;
        if (tmp < 0) { tmp = (CombinationMod - tmp) % CombinationMod; }
        return tmp;
    }
	//start~endまでのn+iCiの和を出力
	long long int combRange1(long long int n,long long int start,long long int end){
		if(start!=0)return (comb(n+end+1,end)-comb(n+start,start-1)+3*CombinationMod)%CombinationMod;
		else return comb(n+end+1,end)%CombinationMod;
	}
	//start1~end1のiとstart2~end2のjまでのi+jCiの和を出力
	long long int combRange2(long long int start1,long long int end1,long long int start2,long long int end2){
		return (comb(end1+end2+2,end1+1)-comb(end1+start2+1,start2)-comb(end2+start1+1,start1)+comb(start1+start2,start1)+5*CombinationMod)%CombinationMod;		
	}
	
	long long int Lucas(long long int n,long long int r){
        long long int ret=1;
        while(r!=0){
            ret*=comb(n%CombinationMod,r%CombinationMod);
            n/=CombinationMod;
            r/=CombinationMod;
        }
        return ret%CombinationMod;
    }
};	



//行列累乗.実行前にgyouretuKouseiを実行すること
class gyouretuRuizyou{
    public:
    long long int ***gyouretu;
    long long int **gyouretuRes,**gyouretuTmp;
    long long int *kitei,*kiteiRes;
    long long int gyouretuMod,gyouretuSize,ruiMax;
    gyouretuRuizyou(long long int inputGyouretuSize,long long int inputRui,long long int inputMod){
        gyouretu=new long long int**[inputRui+1];
        kitei=new long long int[inputGyouretuSize];
        kiteiRes=new long long int[inputGyouretuSize];
        gyouretuRes=new long long int*[inputGyouretuSize];
        gyouretuTmp=new long long int*[inputGyouretuSize];

        for (long long int i = 0; i <= inputRui; i++)
        {
            gyouretu[i]=new long long int*[inputGyouretuSize];
            for (long long int j = 0; j < inputGyouretuSize; j++)
            {
                gyouretu[i][j]=new long long int[inputGyouretuSize];
            }
            
        }

        for (long long int i = 0; i < inputGyouretuSize; i++)
        {
            gyouretuRes[i]=new long long int[inputGyouretuSize];
            gyouretuTmp[i]=new long long int[inputGyouretuSize];
        }
        ruiMax=inputRui;
        gyouretuMod=inputMod;
        gyouretuSize=inputGyouretuSize;
        return ;
    }

    long long int inputGyouretuData(long long int inputI,long long int inputJ,long long int inputData){
        gyouretu[0][inputI][inputJ]=inputData;
        return 0;
    }

    long long int inputKiteiData(long long int inputI,long long int inputData){
        kitei[inputI]=inputData;
        return 0;
    }

    long long int gyouretuKousei(){
        for (long long int i = 1; i <= ruiMax; i++)
        {
            for (long long int j = 0; j < gyouretuSize; j++)
            {
                for (long long int k = 0; k < gyouretuSize; k++)
                {
                    gyouretu[i][j][k]=0;
                    for (long long int p = 0; p < gyouretuSize; p++)
                    {
                        gyouretu[i][j][k]+=gyouretu[i-1][j][p]*gyouretu[i-1][p][k];
                        gyouretu[i][j][k]%=gyouretuMod;
                    }
                    
                }
                
            }
            
        }
        return 0;
    }
    long long int doGyouretuRuizyou(long long int inputRui){
        
        for (long long int i = 0; i < gyouretuSize; i++)
        {
            
            for (long long int j = 0; j < gyouretuSize; j++)
            {
                
                if(i==j){
                    gyouretuRes[i][j]=1;
                }else{
                    gyouretuRes[i][j]=0;
                }
                
            }
            
        }
        
        for (long long int i = 0; i <= ruiMax; i++)
        {
            if(inputRui%2==1){
                for (long long int j = 0; j < gyouretuSize; j++)
                {
                
                   for (long long int k = 0; k < gyouretuSize; k++)
                   {
                       gyouretuTmp[j][k]=0;
                       for (long long int p = 0; p < gyouretuSize; p++)
                       {
                           gyouretuTmp[j][k]+=gyouretuRes[j][p]*gyouretu[i][p][k];
                           gyouretuTmp[j][k]%=gyouretuMod;
                       }
                   }
                   
                }

                for (long long int j = 0; j < gyouretuSize; j++)
                {
                
                   for (long long int k = 0; k < gyouretuSize; k++)
                   {
                       gyouretuRes[j][k]=gyouretuTmp[j][k];
                   }
                   
                }
                
            }

            inputRui/=2;
        }

        
        for (long long int i = 0; i < gyouretuSize; i++)
        {
            kiteiRes[i]=0;
            for (long long int j = 0; j < gyouretuSize; j++)
            {
                kiteiRes[i]+=gyouretuRes[i][j]*kitei[j];
                kiteiRes[i]%=mod;
            }
            
        }
        
        return 0;
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
		size =new long long int[MAX+1];
		for (int i = 0; i < MAX; i++) {
			par[i] = i;
            size[i]=1;
		}
	}
	long long int * par;
    long long int * size;
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
	void unifind_unite(int x, int y) {
		x = unifind_root(x);
		y = unifind_root(y);
		if (x == y) {
			return ;
		}
		par[x] = y;
        size[y] +=size[x];
        return;
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
		return 0;
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
	~BreadthFirstSearch(){
	for(int i=0;i<MaxH;i++)
	{
			free(HW[i]);
			free(Used[i]);
			free(BFSResult[i]);
	}
	free(HW);
	free(Used);
	free(BFSResult);
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
	vector<pair<pair<long long int,long long int>,long long int> >*FlowEdge;
	
    long long int *InEdge;
	long long int MaxV,MaxE;
	long long int *BellmanFordResult,*DijkstraResult,**WarshallFloydResult,*BFSResult;
	vector<long long int>TopologicalSortdVertex;
	bool *InfiniteUpdate;
	bool NegativeLoop=false;
	long long int PrimResult=0;
	bool NotConected=false;
    bool *FlowUsed;

	Graph(long long int InputV,long long int InputE){
		MaxV=InputV;
		MaxE=InputE;
		Edge=new vector<pair<long long int,long long int> >[InputV+1];

	}
	
	void BFS(long long int Start){
		BFSResult=new long long int[MaxV];
		for (long long int i = 0; i < MaxV; i++)
		{
			BFSResult[i]=INF;
		}
		BFSResult[Start]=0;
		queue<long long int>que;
		que.push(Start);
		while(!que.empty()){
			for (long long int i = 0; i < Edge[que.front()].size(); i++)
			{
				if(BFSResult[Edge[que.front()][i].first]==INF){
					BFSResult[Edge[que.front()][i].first]=BFSResult[que.front()]+1;
					que.push(Edge[que.front()][i].first);
				}
			}
			
		}
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
		for (int k = 0; k < MaxV; k++) {
			for (int i = 0; i < MaxV; i++) {
				WarshallFloydResult[k][i] = INF;
				
			}
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
    //おそらく未完成なトポロジカルソート
	void TopologicalSort(){
		bool notupdate=false;
		InEdge=new long long int[MaxV];
		TopologicalSortdVertex.clear();
		queue<long long int>task,firsttask;
		for(int i=0;i<MaxV;i++){
			InEdge[i]=0;
		}
		for(int i=0;i<MaxV;i++){
			for(int j=0;j<Edge[i].size();j++){
				InEdge[Edge[i][j].first]++;
			}
		}
		for(int i=0;i<MaxV;i++){
			if(InEdge[i]==0){
				firsttask.push(i);
			}

		}
		while(!firsttask.empty()){
			long long int i=firsttask.front();
			firsttask.pop();
			TopologicalSortdVertex.push_back(i);			
			InEdge[i]=-1;
			for(int j=0;j<Edge[i].size();j++){
				InEdge[Edge[i][j].first]--;
				if(InEdge[Edge[i][j].first]==0){
					task.push(Edge[i][j].first);
				}
			}
			

		}
		while(!task.empty()){
			long long int tmp=task.front();
			TopologicalSortdVertex.push_back(tmp);
			task.pop();
			if(InEdge[tmp]==0){
				InEdge[tmp]=-1;
				for(int j=0;j<Edge[tmp].size();j++){
					InEdge[Edge[tmp][j].first]--;
					if(InEdge[Edge[tmp][j].first]==0){
						task.push(Edge[tmp][j].first);
					}
				}
			}
			
		}
	}
    //これを最初に記述
    void NetworkFlow(){
        FlowEdge=new vector<pair<pll,long long int> >[MaxV];
        FlowUsed=new bool[MaxV];
        for (long long int i = 0; i < MaxV; i++)
        {
            FlowUsed[i]=false;
        }
        
    }
    void NetworlFlowAddEdge(long long int from,long long int to,long long int cap){
        FlowEdge[from].push_back({{to,cap},FlowEdge[to].size()});
        FlowEdge[to].push_back({{from,0},FlowEdge[from].size()-1});
    }
    long long int Flowdfs(long long int v,long long int t,long long int f){
        if(v==t){return f;}
        FlowUsed[v]=true;
        for (long long int i = 0; i < FlowEdge[v].size(); i++)
        {
            pair<pll,long long int> &e=FlowEdge[v][i];
            if(!FlowUsed[e.first.first]&&e.first.second>0){
                long long int d=Flowdfs(e.first.first,t,min(f,e.first.second));
                if(d>0){
                    e.first.second-=d;
                    FlowEdge[e.first.first][e.second].first.second+=d;
                    return d;
                }
            }
        }
        return 0;
    }
    long long int MaxFlow(long long int s,long long int t){
        long long int flow=0;
        for (;;)
        {
            for (long long int i = 0; i < MaxV; i++)
            {
                FlowUsed[i]=false;
            }
            
            long long int f=Flowdfs(s,t,INF);
            
            if(f==0){return flow;}
            flow+=f;
        }
        
    }
};

class LCS{
    public:
    long long int **dp;
    long long int **longest;
    long long int result;
    string LCSstr;
    LCS(long long int one,long long int two){
        dp=new long long int*[one];
        longest=new long long int *[one];
        for(int i=0;i<one;i++){
            dp[i]=new long long int[two];
            longest[i]=new long long int[two];
        } 
    }
    void DoLCS(string str1,string str2){
        LCSstr.clear();
        for(int i=0;i<str1.length();i++){
            for(int j=0;j<str2.length();j++){
                dp[i][j]=0;
                longest[i][j]=0;
            }
            
        }

        for(int i=0;i<str1.length();i++){
            if(str1[i]==str2[0]){dp[i][0]=1;}
            else{dp[i][0]=0;}
            for(int j=1;j<str2.length();j++){
                if(i!=0)dp[i][j]=longest[i-1][j-1];
                if(str1[i]==str2[j]){dp[i][j]++;}
            }
            
            if(i!=0){
                longest[i][0]=max({longest[i-1][0],dp[i][0]});
                for(int j=1;j<str2.length();j++){
                    longest[i][j]=max({longest[i-1][j],dp[i][j],longest[i][j-1]});
                }
            }else{
                longest[0][0]=dp[i][0];
                for(int j=1;j<str2.length();j++){
                    longest[0][j]=max({dp[0][j],longest[0][j-1]});
                }
            }
        }

        pair<long long int ,long long int>rev;
        rev={str1.length(),str2.length()};
        long long int now=longest[str1.length()-1][str2.length()-1];
        bool key=false;
        //cout<<now<<endl;
        //DEBUG_OUTPUT_ARRAY2_BOX(dp,rev.first,rev.second);
        while(now!=0){
            key=false;
            for(int j=0;j<rev.first;j++){
                if(now==dp[j][rev.second-1]&&(now-1==longest[max(j-1,0)][max(rev.second-1-1,0LL)]||j==0||rev.second-1==0)){
                    
                    LCSstr.push_back(str1[j]);
                    rev={j,rev.second-1};
                    now--;
                    key=true;
                    //cout<<"REV>>"<<rev.first<<" "<<rev.second<<endl;
                    break;
                }
            }
            if(key){continue;}
            for(int j=0;j<rev.second;j++){
                if(now==dp[rev.first-1][j]&&(now-1==longest[max(rev.first-1-1,0LL)][max(j-1,0)]||rev.first-1==0||j==0)){
                    now--;
                    LCSstr.push_back(str2[j]);
                    rev={rev.first-1,j};
                    key=true;
                    
                    //cout<<"REV>>"<<rev.first<<" "<<rev.second<<endl;
                    break;
                }
            }
            if(!key){
                rev={rev.first-1,rev.second-1};
            }        
        }
        reverse(LCSstr.begin(),LCSstr.end());
        result=LCSstr.length();
    }
    ~LCS(){}
};

//自分以下の数で自分と互いに素なものの個数を出力O(√N)
long long int EularTotientFunction(long long int Input){
    long long int InputCopy=Input;
    for(long long int i = 2; i*i <= Input; i++){
        if(Input%i==0){
            while(Input%i==0){
                Input/=i;
            }
            InputCopy=InputCopy*(i-1)/i;
        }
    }
    if(Input!=1){
        InputCopy=InputCopy*(Input-1)/Input;
    }
    return InputCopy;
}