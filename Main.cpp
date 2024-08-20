//  Gabriela Aparecida Silva Cetano - 202110681
//  João Augusto Dias Neto - 202110228
//  Mathias Silva Sousa - 201920352

#include <iostream>
#include <cstring>
#include <vector>
#include <string>
#include <sstream>
#include <algorithm>
#include <queue>
#include <set>

using namespace std;

#define DIRECIONADO 1
#define NAO_DIRECIONADO 0

typedef vector<vector<pair<int, pair<int, int>>>> grafo;


grafo lista_adj, grafo_associado, lista_adj_rev;

vector<vector<int>> capacidade;

vector< pair<long long, pair<int, int>> > edgeList;

vector<int> funcoes, visitado;

int n_arestas, n_vertices;

bool b_direcionado;

string direcionado;

void leitura_direcionado() {
    // Itera sobre o número de arestas do grafo
    for (int i = 0; i < n_arestas; ++i) {
        int id, a, b, p;
        cin >> id >> a >> b >> p;  // Lê os valores de entrada: id da aresta, vértice de origem, vértice de destino e peso

        // Adiciona a aresta na lista de adjacências do grafo direcionado
        lista_adj[a].push_back({b, {p, id}});

        // Adiciona a aresta no grafo associado, que é não direcionado (armazena em ambos os sentidos)
        grafo_associado[a].push_back({b, {p, id}});
        grafo_associado[b].push_back({a, {p, id}});

        // Adiciona a aresta na lista de adjacências do grafo reverso (para algoritmos que precisem do grafo transposto)
        lista_adj_rev[b].push_back({a, {p, id}});

        // Define a capacidade da aresta no grafo direcionado (para uso em fluxos, por exemplo)
        capacidade[a][b] = p;
    } 
}


void leitura_nao_direcionado() {
    for (int i = 0; i < n_arestas; ++i) {
        int id, a, b, p;
        cin >> id >> a >> b >> p;
        
        // Adiciona a aresta nas duas direções na lista de adjacência, pois o grafo não é direcionado.
        lista_adj[a].push_back({b, {p, id}});
        lista_adj[b].push_back({a, {p, id}});
        
        // Armazena a aresta no grafo associado, em ambos os sentidos.
        grafo_associado[a].push_back({b, {p, id}});
        grafo_associado[b].push_back({a, {p, id}});
        
        // Adiciona a aresta na lista de arestas (edgeList) com seus vértices e peso.
        edgeList.push_back(make_pair(p, pair<int, int>(a, b)));
    } 
}


void dfs(int u, bool print) {
    // Marca o vértice 'u' como visitado.
    visitado[u] = 1;
    
    // Itera sobre todos os vértices adjacentes ao vértice 'u' no grafo associado.
    for (auto &v: grafo_associado[u]) {
        // Se o vértice 'v' ainda não foi visitado, prossegue com a busca em profundidade.
        if (!visitado[v.first]) {
            // Se a flag 'print' estiver ativada, imprime o vértice 'v'.
            if (print) {
                cout << v.first << " ";
            }
            // Chama recursivamente a função DFS para o vértice 'v'.
            dfs(v.first, print);
        }
    }
}



bool conexo() {
    visitado.resize(n_vertices, 0);
    dfs(0, false);
    for (int i = 1; i < n_vertices; ++i) {
        if (!visitado[i]) return false;
    }
    return true;
}

// para ser bipartido, o grafo deve ser colorivel com 2 cores
// ou seja, se um vertice V possui cor 1, todos os seus vizinhos devem ter cor 0
// para isso, basta fazer um BFS verificando se os vizinhos possuem cor diferente
// retorna 1 caso seja bipartido, 0 caso contrario
bool bipartido() {
    vector<int> side(n_vertices, -1);
    bool eh_bipartido = true;
    queue<int> q;
    for (int src = 0; src < n_vertices; ++src) {
        //se nao foi pintado
        if (side[src] == -1) {
            //entra na fila
            q.push(src);
            side[src] = 0;
            //BFS
            while (!q.empty()) {
                int v = q.front();
                q.pop();
                for (auto p : lista_adj[v]) {
                    int u = p.first;
                    //se ainda nao foi pintado
                    if (side[u] == -1) {
                        //vizinho recebe cor oposta
                        side[u] = side[v] ^ 1;
                        q.push(u);
                    } else {
                        //se a cor for igual, eh_bipartido vai ser false pra sempre
                        eh_bipartido &= side[u] != side[v];
                    }
                }
            }
        }
    }
    return eh_bipartido;
}

//para grafo nao direcionado, basta checar se todos os vertices possuem grau par
// para direcionado, grau de entrada (indeg) deve ser igual ao grau de saida (outdeg)
// retorna 1 caso seja euleriano, 0 caso contrario
bool euleriano() {
    if (!conexo()) return false;
    if (!b_direcionado) {
        for (int i = 0; i < n_vertices; ++i) {
            if (lista_adj[i].size() % 2 != 0) return false;
        }
        return true;
    }
    vector<int> indeg(n_vertices, 0), outdeg(n_vertices, 0);
    
    for (int i = 0; i < n_vertices; ++i) {
        
        outdeg[i] = lista_adj[i].size();
        for (int j = 0; j < lista_adj[i].size(); ++j) {
            int v = lista_adj[i][j].first;
            indeg[v]++;
        }
        
    }
    for (int i = 0; i < n_vertices; ++i) {
        if (indeg[i] != outdeg[i]) return false;
    }
    return true;
}

int conta_ciclo;

#define NAO_VISITADO -1
#define EXPLORADO -2
#define VISITADO -3

vector<int> vis, pai;

void dfs_cycle(int u, int chamou) { 
	vis[u] = EXPLORADO;
	for (int j = 0, v; j < (int)lista_adj[u].size(); j++) {
		v = lista_adj[u][j].first;
        if (v == chamou) continue; // para evitar loop infinito
		if (vis[v] == NAO_VISITADO) { 
			pai[v] = u; 
			dfs_cycle(v, u);
		}
		else if (vis[v] == EXPLORADO) {
			conta_ciclo++; // aresta de retorno, ou seja, ciclo
		}
	}
	vis[u] = VISITADO; 
}

//retorna verdadeiro se existe pelo menos um ciclo
bool ciclo() {
    vis.resize(n_vertices, NAO_VISITADO);
    pai.resize(n_vertices);
    for (int i = 0; i < n_vertices; ++i) {
        if (vis[i] == NAO_VISITADO) {
            dfs_cycle(i, i);
        }
    }
    return conta_ciclo > 0; 
}
// retorna a quantidade de componentes conexas 
int componentes_conexas() {
    int qtd = 0;
    for (int i = 0; i < n_vertices; ++i) visitado[i] = 0;
    for (int i = 0; i < n_vertices; ++i) {
        if (!visitado[i]) {
            qtd++;
            dfs(i, false);
        }
    }
    return qtd;
} 


int parent = 0, numSCC;
vector<int>comp, ts;

void revdfs(int u) {
    vis[u] = true;
    for(int i = 0, v; i < (int)lista_adj_rev[u].size(); i++) {
    	v = lista_adj_rev[u][i].first;
        if(!vis[v]) revdfs(v);
    }
    ts.push_back(u);
}

void dfs(int u) {
    vis[u] = true; comp[u] = parent;
    for(int i = 0, v; i < (int)lista_adj[u].size(); i++) {
    	v = lista_adj[u][i].first;
        if(!vis[v]) dfs(v);
    }
}

void kosaraju() {
    // Redimensiona os vetores 'comp' e 'vis' para o número de vértices.
    comp.resize(n_vertices);
    vis.resize(n_vertices);
    
    // Inicializa o vetor 'vis' com 0 (não visitado).
    fill(vis.begin(), vis.end(), 0);
    
    // Realiza a primeira etapa de Kosaraju: DFS no grafo original para preencher a ordem de finalização dos vértices.
    for (int i = 0; i < n_vertices; i++) {
        if (!vis[i]) {
            revdfs(i);
        }
    }
    
    // Reinicializa o vetor 'vis' para a segunda etapa.
    fill(vis.begin(), vis.end(), 0);
    
    // Inicializa o contador de componentes fortemente conexos (SCC).
    numSCC = 0;
    
    // Realiza a segunda etapa de Kosaraju: DFS no grafo transposto (invertido) na ordem de finalização dos vértices.
    for (int i = n_vertices - 1; i >= 0; i--) {
        if (!vis[ts[i]]) {
            parent = ts[i];
            dfs(ts[i]);
            numSCC++; // Incrementa o contador de SCCs.
        }
    }
}


int counter, rootChildren, root, temArticulacao;
vector<int> num, low, parente, articulationVertex;
int pontes;

void tarjan(int u) {
    // Inicializa 'low[u]' e 'num[u]' com o valor do contador e incrementa o contador.
    low[u] = num[u] = counter++;
    
    // Itera sobre todos os vértices adjacentes a 'u'.
    for (int j = 0, v; j < (int)lista_adj[u].size(); j++) {
        v = lista_adj[u][j].first;
        
        // Se 'v' ainda não foi visitado.
        if (num[v] == NAO_VISITADO) {
            // Define 'u' como pai de 'v'.
            parente[v] = u;
            if (u == root) rootChildren++; // Conta filhos do root.
            
            // Chama a função recursivamente para 'v'.
            tarjan(v);
            
            // Se o menor número de 'v' não é menor que o número de 'u', 'u' é um ponto de articulação.
            if (low[v] >= num[u]) articulationVertex[u] = true;
            
            // Se o menor número de 'v' é maior que o número de 'u', incrementa o contador de pontes.
            if (low[v] > num[u]) pontes++;
            
            // Atualiza 'low[u]' com o menor valor entre 'low[u]' e 'low[v]'.
            low[u] = min(low[u], low[v]);
        }
        // Se 'v' já foi visitado e não é o pai de 'u'.
        else if (v != parente[u]) {
            // Atualiza 'low[u]' com o menor valor entre 'low[u]' e 'num[v]'.
            low[u] = min(low[u], num[v]);
        }
    }
}


int qtd_articulacoes;

int articulacoes(int print) {
    counter = 0;
    qtd_articulacoes = 0;
    num.resize(n_vertices, NAO_VISITADO);
    low.resize(n_vertices, 0);
    parente.resize(n_vertices, 0);
    articulationVertex.resize(n_vertices, 0);
    for (int i = 0; i < n_vertices; i++)
		if (num[i] == NAO_VISITADO) {
			root = i; rootChildren = 0; tarjan(i);
			articulationVertex[root] = (rootChildren >= 1);
		}
    if (print) {
        for (int i = 0; i < n_vertices; ++i) if (articulationVertex[i]) qtd_articulacoes++;
        if (qtd_articulacoes) {
            for (int i = 0; i < n_vertices; ++i) {
                if (articulationVertex[i]) {
                    cout << i << ' ';
                }
            }
            cout << endl;
        }
        else cout << -1 << endl;

    }
    return pontes;
}

void dfs_tree(int u) {
    visitado[u] = 1;
    for (auto &v: lista_adj[u]) {
        if (!visitado[v.first]) {
            cout << v.second.second << " ";
            dfs_tree(v.first) ;
        }
    }
}

void arvore_dfs() {
    fill(visitado.begin(), visitado.end(), 0);
    dfs_tree(0);
    cout << endl;
}

void bfs_tree() {
    vector<int> visitado(n_vertices, -1);
    queue<int> q;
    int src = 0;
    q.push(src);
    visitado[src] = 0;
    //BFS
    while (!q.empty()) {
        int v = q.front();
        q.pop();
        for (auto p : lista_adj[v]) {
            int u = p.first;
            if (visitado[u] == -1) {
                visitado[u] = visitado[v] + 1;
                cout << p.second.second << " ";
                q.push(u);
            } 
        }
    }
    cout << endl;
    
}

class UnionFind {
private:
	vector<int> parent, rank;
public:
	UnionFind(int N) {
		rank.assign(N+9, 0);
		parent.assign(N+9, 0);
		for (int i = 0; i < N; i++) parent[i] = i;
	}
	int find(int i) {
		while(i != parent[i]) i = parent[i];
		return i;
	}
	bool isSameSet(int i, int j) {
		return find(i) == find(j);
	}
	void unionSet (int i, int j) {
		if (isSameSet(i, j)) return;
		int x = find(i), y = find(j);
		if (rank[x] > rank[y]) parent[y] = x;
		else {
			parent[x] = y;
			if (rank[x] == rank[y]) rank[y]++;
		}
	}
};

typedef pair<int, int> ii;
typedef long long ll;


ll kruskal() {
	ll cost = 0;
	UnionFind UF(n_vertices);
	pair<int, ii> edge;
	sort(edgeList.begin(), edgeList.end());
	for (int i = 0; i < n_arestas; i++) {
		edge = edgeList[i];
		if (!UF.isSameSet(edge.second.first, edge.second.second)) { 
			cost += edge.first;
			UF.unionSet(edge.second.first, edge.second.second);
		}
	}
	return cost;
}

vector<int> toposort;	//Ordem reversa!

void topo(int u) {
	vis[u] = true;
	for (int j = 0, v; j < (int)lista_adj[u].size(); j++) {
		v = lista_adj[u][j].first;
		if (!vis[v]) topo(v);
	}
	toposort.push_back(u);
    cout << u << " "; 
}

void topologicalsort() {
    vis.assign(n_vertices, 0);
    for (int i = 0; i < n_vertices; ++i) {
        if (!vis[i]) {
            topo(i);
        }
    }
    cout << endl;
}

#define INF 0x3f3f3f3f

int dijkstra(int s, int t) {
    vector<int> dist(n_vertices, INF);
	set<ii> pq;
	dist[s] = 0;
	pq.insert(ii(0, s));
	while(!pq.empty()) {
		int u = pq.begin()->second;
		pq.erase(pq.begin());
		for(int i=0; i<(int)lista_adj[u].size(); i++) {
			int v = lista_adj[u][i].first;
			int w = lista_adj[u][i].second.first;
			if (dist[v] > dist[u] + w) {
				pq.erase(ii(dist[v], v));
				dist[v] = dist[u] + w;
				pq.insert(ii(dist[v], v));
			}
		}
	}
	return dist[t];
}


// vector<vector<int>> capacidade;

int bfs(int s, int t, vector<int>& parent) {
    fill(parent.begin(), parent.end(), -1);
    parent[s] = -2;
    queue<pair<int, int>> q;
    q.push({s, INF});

    while (!q.empty()) {
        int atual = q.front().first;
        int fluxo = q.front().second;
        q.pop();

        for (auto vertice : lista_adj[atual]) {
            int prox = vertice.first;
            if (parent[prox] == -1 && capacidade[atual][prox]) {
                parent[prox] = atual;
                int novo_fluxo = min(fluxo, capacidade[atual][prox]);
                if (prox == t)
                    return novo_fluxo;
                q.push({prox, novo_fluxo});
            }
        }
    }

    return 0;
}

int maxflow(int s, int t) {
    int fluxo = 0;
    vector<int> parent(n_vertices);
    int novo_fluxo;

    while (novo_fluxo = bfs(s, t, parent)) {
        fluxo += novo_fluxo;
        int atual = t;
        while (atual != s) {
            int anterior = parent[atual];
            capacidade[anterior][atual] -= novo_fluxo;
            capacidade[atual][anterior] += novo_fluxo;
            atual = anterior;
        }
    }

    return fluxo;
}

void fecho() {
    visitado.assign(n_vertices, 0);
    dfs(0, true);
    cout << endl;
}

int main () {

    //leitura das funcoes a serem testadas
    string the_string;
    getline(cin, the_string);
    istringstream iss(the_string);
    for (int funcao; iss >> funcao; )
    {
        funcoes.push_back(funcao);
    }
    sort(funcoes.begin(), funcoes.end());

    cin >> n_vertices >> n_arestas;
    cin >> direcionado;
    capacidade.resize(n_vertices);
    for (int i = 0; i < n_vertices; ++i) {
        capacidade[i].assign(n_vertices, 0);
    }
    lista_adj.resize(n_vertices);
    grafo_associado.resize(n_vertices);
    visitado.resize(n_vertices);
    lista_adj_rev.resize(n_vertices);

    //resetando o vetor de vertices visitados
    for (int i = 0; i < n_vertices; ++i) visitado[i] = 0;
    

    //controle para saber se eh direcionado ou nao
    if (direcionado.compare("direcionado") == 0) {
        b_direcionado = DIRECIONADO;
        leitura_direcionado();
    }
    else {
        b_direcionado = NAO_DIRECIONADO;
        leitura_nao_direcionado();
    }

    //ordenando os vertices vizinhos de cada vertice
    for (int i = 0; i < n_vertices; ++i) sort(lista_adj[i].begin(), lista_adj[i].end());

    for (int i = 0; i < funcoes.size(); ++i) {
        // cout << "f: " <<  funcoes[i] << endl; //apagar futuramente
        //funcoes indexadas em 0, assim como nos casos
        switch (funcoes[i])
        {
            case 0:
                cout << conexo() << endl;
                break;
            case 1:
                cout << bipartido() << endl;
                break;
            case 2:
                cout << euleriano() << endl;
                break;
            case 3:
                cout << ciclo() << endl;
                break;
            case 4:
                if (b_direcionado) {
                    cout << -1 << endl;
                    break;
                }
                cout << componentes_conexas() << endl;
                break;
            case 5:
                if (!b_direcionado) {
                    cout << -1 << endl;
                    break; 
                }
                kosaraju();
                cout << numSCC << endl;
                break;
            case 6:
                if (b_direcionado) {
                    cout << -1 << endl;
                    break;
                }
                articulacoes(true);
                break;
            case 7:
                if (b_direcionado) {
                    cout << -1 << endl;
                    break;
                }
                cout << articulacoes(false) << endl;
                break;
            case 8:
                arvore_dfs();
                break;
            case 9:
                bfs_tree();
                break;
            case 10:
                if (b_direcionado || !conexo()) {
                    cout << -1 << endl;
                    break;
                }
                cout << kruskal () << endl;
                break;
            case 11:
                if (!b_direcionado || ciclo()) {
                    cout << -1 << endl;
                    break;
                }
                topologicalsort();
            case 12:
                if (b_direcionado) {
                    cout << -1 << endl;
                    break;
                }
                cout << dijkstra(0, n_vertices-1) << endl;
                break;
            case 13:
                if (!b_direcionado) {
                    cout << -1 << endl;
                    break;
                }
                cout << maxflow(0, n_vertices-1) << endl;
                break;
            case 14:
                if (!b_direcionado) {
                    cout << -1 << endl;
                    break;
                }
                fecho();
                break;
            default:
                break;
        }
    }
    return 0;
}
