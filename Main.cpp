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


//Função para ver se o grafo é conexo (todos os vértices estão conectados)
bool conexo() {
    visitado.resize(n_vertices, 0); // inicializando o vetor visitado, ele é redimensionado para ter o mesmo tamanho do número de vértices (n_vertices) do grafo
    dfs(0, false); //inicia uma busca em profundidade a partir do vértice '0'
    for (int i = 1; i < n_vertices; ++i) {
        if (!visitado[i]) return false; // após realizar a DFS, a função verifica se todos os vértices foram visitados. Se encontrar algum vértice que não foi visitado, retorna false, indicando que o grafo não é conexo.
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
bool euleriano(){
    if (!conexo()) return false; // verifica se o grafo é conexo, se não for conexo, não tem como ser eureliano
    if (!b_direcionado) /* Se o grafo for não direcionado */{
        for (int i = 0; i < n_vertices; ++i) {
            if (lista_adj[i].size() % 2 != 0) return false; //Verifica o grau de cada vértice, para um grafo ser eureliano, todos os vértices devem ter grau par. Se qualquer vértice tiver um grau ímpar, retorna false
        }
        return true;
    }
    //Se o grafo for direcionado
    vector<int> indeg(n_vertices, 0), outdeg(n_vertices, 0); ////Inicializa os vetores para graus de Entrada e Saída de cada vértice
    
    for (int i = 0; i < n_vertices; ++i) {
        
        outdeg[i] = lista_adj[i].size(); ////outdeg[i] é  o tamanho da lista de adjacências de i, o que representa o grau de saída do vértice i
        for (int j = 0; j < lista_adj[i].size(); ++j) {
            int v = lista_adj[i][j].first;
            indeg[v]++; //Percorre a lista de adjacências de cada vértice i. Para cada vértice adjacente v, o grau de entrada de v (indeg[v]) é incrementado.
        }
        
    }
    for (int i = 0; i < n_vertices; ++i) /*Verificação dos graus de entrada e saída*/ {
        if (indeg[i] != outdeg[i]) return false; //Para um grafo direcionado ser euleriano, o grau de entrada deve ser igual ao grau de saída para todos os vértices. 
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
    // Inicializa os contadores e variáveis necessárias
    counter = 0;
    qtd_articulacoes = 0;
    
    // Redimensiona os vetores para o número de vértices do grafo
    num.resize(n_vertices, NAO_VISITADO);       // Vetor para armazenar o número de visitação de cada vértice (todos inicialmente não visitados)
    low.resize(n_vertices, 0);                  // Vetor para armazenar o menor número de visita alcançável a partir de um vértice
    parente.resize(n_vertices, 0);              // Vetor para armazenar o pai de cada vértice na árvore de DFS
    articulationVertex.resize(n_vertices, 0);   // Vetor para marcar quais vértices são pontos de articulação (articulation points)
    
    // Itera sobre todos os vértices do grafo
    for (int i = 0; i < n_vertices; i++) {
        if (num[i] == NAO_VISITADO) {            // Se o vértice ainda não foi visitado
            root = i;                            // Define o vértice atual como a raiz da árvore DFS
            rootChildren = 0;                    // Inicializa o contador de filhos da raiz
            tarjan(i);                           // Chama a função Tarjan para processar o vértice
            
            // Verifica se a raiz é um ponto de articulação
            articulationVertex[root] = (rootChildren >= 1);
        }
    }

    // Se a função foi chamada com print != 0, imprime os pontos de articulação
    if (print) {
        // Conta quantos pontos de articulação existem
        for (int i = 0; i < n_vertices; ++i) {
            if (articulationVertex[i]) {
                qtd_articulacoes++;
            }
        }

        // Se houver pelo menos um ponto de articulação, imprime-os
        if (qtd_articulacoes) {
            for (int i = 0; i < n_vertices; ++i) {
                if (articulationVertex[i]) {
                    cout << i << ' ';  // Imprime o índice do ponto de articulação
                }
            }
            cout << endl;
        } else {
            cout << -1 << endl;  // Se não houver pontos de articulação, imprime -1
        }
    }

    // Retorna a quantidade de pontes no grafo (ou outro valor dependendo da implementação)
    return pontes;
}


void dfs_tree(int u) {
    visitado[u] = 1;  // Marca o vértice atual 'u' como visitado

    // Itera sobre a lista de adjacência do vértice 'u'
    for (auto &v: lista_adj[u]) {
        if (!visitado[v.first]) {  // Se o vértice 'v.first' não foi visitado
            cout << v.second.second << " ";  // Imprime o segundo valor associado à aresta (presumivelmente uma etiqueta ou peso)
            dfs_tree(v.first);  // Chama recursivamente a função DFS para explorar o vértice 'v.first'
        }
    }
}


void arvore_dfs() {
    // Inicializa o vetor 'visitado' com zeros, marcando todos os vértices como não visitados
    fill(visitado.begin(), visitado.end(), 0);
    
    // Inicia a busca em profundidade (DFS) a partir do vértice 0
    dfs_tree(0);
    
    // Imprime uma nova linha após a execução da DFS
    cout << endl;
}


void bfs_tree() {
    // Inicializa o vetor 'visitado' com -1, marcando todos os vértices como não visitados
    vector<int> visitado(n_vertices, -1);
    
    // Cria uma fila para o BFS e insere o vértice de origem (src = 0)
    queue<int> q;
    int src = 0;
    q.push(src);
    visitado[src] = 0;  // Marca o vértice de origem como visitado
    
    // Início do loop BFS (Busca em Largura)
    while (!q.empty()) {
        int v = q.front();  // Obtém o vértice na frente da fila
        q.pop();  // Remove o vértice da fila
        
        // Itera sobre os vizinhos do vértice 'v'
        for (auto p : lista_adj[v]) {
            int u = p.first;  // Obtém o vértice adjacente
            
            if (visitado[u] == -1) {  // Se o vértice adjacente não foi visitado
                visitado[u] = visitado[v] + 1;  // Marca o vértice 'u' como visitado e define a distância a partir da origem
                cout << p.second.second << " ";  // Imprime o segundo valor associado à aresta (presumivelmente uma etiqueta ou peso)
                q.push(u);  // Adiciona o vértice 'u' na fila para continuar a BFS
            } 
        }
    }
    
    // Imprime uma nova linha após a execução do BFS
    cout << endl;
}


class UnionFind {
private:
    vector<int> parent, rank;  // Vetores para armazenar o pai de cada elemento e a "altura" das árvores

public:
    // Construtor que inicializa o Union-Find para 'N' elementos
    UnionFind(int N) {
        rank.assign(N + 9, 0);   // Inicializa o vetor de rank com zeros, com tamanho N+9
        parent.assign(N + 9, 0); // Inicializa o vetor de pais com zeros, com tamanho N+9
        for (int i = 0; i < N; i++) parent[i] = i;  // Cada elemento é seu próprio pai (representante do conjunto)
    }

    // Função para encontrar o pai de um elemento 'i'
    int find(int i) {
        while(i != parent[i]) i = parent[i];  // Segue o caminho até encontrar o representante
        return i;  // Retorna o representante do conjunto ao qual 'i' pertence
    }

    // Verifica se os elementos 'i' e 'j' pertencem ao mesmo conjunto
    bool isSameSet(int i, int j) {
        return find(i) == find(j);  // Retorna true se ambos tiverem o mesmo representante
    }

    // Função para unir os conjuntos que contêm os elementos 'i' e 'j'
    void unionSet(int i, int j) {
        if (isSameSet(i, j)) return;  // Se já pertencem ao mesmo conjunto, não faz nada
        int x = find(i), y = find(j); // Encontra os representantes dos conjuntos de 'i' e 'j'

        // Realiza a união com base no rank (técnica de union by rank)
        if (rank[x] > rank[y]) parent[y] = x;  // Se o rank de 'x' é maior, 'x' se torna o pai de 'y'
        else {
            parent[x] = y;  // Caso contrário, 'y' se torna o pai de 'x'
            if (rank[x] == rank[y]) rank[y]++;  // Se os ranks eram iguais, incrementa o rank de 'y'
        }
    }
};


typedef pair<int, int> ii;
typedef long long ll;


ll kruskal() {
    ll cost = 0;  // Inicializa o custo total da árvore geradora mínima (MST) como 0
    UnionFind UF(n_vertices);  // Cria uma estrutura de Union-Find para gerenciar os conjuntos de vértices
    pair<int, ii> edge;  // Define uma variável para armazenar uma aresta (peso, vértice 1, vértice 2)

    // Ordena a lista de arestas em ordem crescente com base no peso
    sort(edgeList.begin(), edgeList.end());

    // Itera sobre todas as arestas do grafo
    for (int i = 0; i < n_arestas; i++) {
        edge = edgeList[i];  // Obtém a i-ésima aresta da lista de arestas

        // Verifica se os dois vértices da aresta estão em conjuntos diferentes
        if (!UF.isSameSet(edge.second.first, edge.second.second)) { 
            cost += edge.first;  // Adiciona o peso da aresta ao custo total
            UF.unionSet(edge.second.first, edge.second.second);  // Une os conjuntos dos dois vértices
        }
    }

    return cost;  // Retorna o custo total da árvore geradora mínima
}


vector<int> toposort;  // Vetor que armazenará a ordem topológica dos vértices (em ordem reversa)

// Função recursiva para realizar a ordenação topológica
void topo(int u) {
    vis[u] = true;  // Marca o vértice 'u' como visitado

    // Itera sobre todos os vértices adjacentes a 'u'
    for (int j = 0, v; j < (int)lista_adj[u].size(); j++) {
        v = lista_adj[u][j].first;  // Obtém o vértice adjacente 'v'
        
        // Se o vértice 'v' ainda não foi visitado, chama a função topo para ele
        if (!vis[v]) topo(v);
    }

    // Após visitar todos os vértices adjacentes, adiciona 'u' ao vetor de ordenação topológica
    toposort.push_back(u);

    // Opcional: Imprime o vértice 'u' (na ordem de processamento, que será a ordem inversa da ordem topológica)
    cout << u << " "; 
}


void topologicalsort() {
    // Inicializa o vetor 'vis' com zeros, marcando todos os vértices como não visitados
    vis.assign(n_vertices, 0);
    
    // Percorre todos os vértices do grafo
    for (int i = 0; i < n_vertices; ++i) {
        // Se o vértice 'i' não foi visitado, chama a função 'topo' para ele
        if (!vis[i]) {
            topo(i);
        }
    }
    
    // Imprime uma nova linha após a execução da ordenação topológica
    cout << endl;
}


#define INF 0x3f3f3f3f

int dijkstra(int s, int t) {
    // Inicializa o vetor de distâncias com infinito para todos os vértices
    vector<int> dist(n_vertices, INF);

    // Cria uma fila de prioridade para gerenciar os vértices a serem processados
    set<ii> pq;

    // Define a distância do vértice de origem 's' como 0
    dist[s] = 0;

    // Adiciona o vértice de origem à fila de prioridade com distância 0
    pq.insert(ii(0, s));

    // Executa o algoritmo de Dijkstra
    while(!pq.empty()) {
        // Obtém o vértice com a menor distância da fila de prioridade
        int u = pq.begin()->second;
        pq.erase(pq.begin());

        // Itera sobre todos os vértices adjacentes a 'u'
        for(int i = 0; i < (int)lista_adj[u].size(); i++) {
            int v = lista_adj[u][i].first;  // Vértice adjacente
            int w = lista_adj[u][i].second.first;  // Peso da aresta para o vértice 'v'

            // Se encontrar uma distância menor para 'v', atualiza a distância e ajusta a fila de prioridade
            if (dist[v] > dist[u] + w) {
                // Remove a entrada antiga para o vértice 'v' da fila de prioridade
                pq.erase(ii(dist[v], v));
                // Atualiza a distância para o vértice 'v'
                dist[v] = dist[u] + w;
                // Adiciona o vértice 'v' atualizado à fila de prioridade
                pq.insert(ii(dist[v], v));
            }
        }
    }

    // Retorna a menor distância do vértice de origem 's' ao vértice destino 't'
    return dist[t];
}


int bfs(int s, int t, vector<int>& parent) {
    // Inicializa o vetor 'parent' com -1, indicando que nenhum vértice foi visitado
    fill(parent.begin(), parent.end(), -1);

    // Marca o vértice de origem 's' com -2 para indicar que é o ponto de partida
    parent[s] = -2;

    // Cria uma fila para a busca em largura (BFS) e insere o vértice de origem 's' com fluxo infinito
    queue<pair<int, int>> q;
    q.push({s, INF});

    // Executa o BFS
    while (!q.empty()) {
        // Obtém o vértice atual e o fluxo atual da frente da fila
        int atual = q.front().first;
        int fluxo = q.front().second;
        q.pop();

        // Itera sobre todos os vértices adjacentes ao vértice atual
        for (auto vertice : lista_adj[atual]) {
            int prox = vertice.first;  // Próximo vértice
            // Verifica se o próximo vértice ainda não foi visitado e se há capacidade restante
            if (parent[prox] == -1 && capacidade[atual][prox]) {
                // Atualiza o pai do próximo vértice e calcula o novo fluxo
                parent[prox] = atual;
                int novo_fluxo = min(fluxo, capacidade[atual][prox]);
                // Se o próximo vértice é o destino 't', retorna o novo fluxo
                if (prox == t)
                    return novo_fluxo;
                // Caso contrário, insere o próximo vértice na fila com o novo fluxo
                q.push({prox, novo_fluxo});
            }
        }
    }

    // Se o destino 't' não for alcançado, retorna 0
    return 0;
}


int maxflow(int s, int t) {
    int fluxo = 0;  // Inicializa o fluxo total como 0
    vector<int> parent(n_vertices);  // Vetor para armazenar o caminho encontrado pelo BFS
    int novo_fluxo;  // Variável para armazenar o fluxo encontrado na iteração do BFS

    // Enquanto houver um caminho de fluxo disponível do vértice 's' ao vértice 't'
    while (novo_fluxo = bfs(s, t, parent)) {
        fluxo += novo_fluxo;  // Adiciona o fluxo encontrado ao fluxo total

        int atual = t;  // Começa do vértice de destino 't'
        // Atualiza as capacidades das arestas no caminho encontrado
        while (atual != s) {
            int anterior = parent[atual];  // Obtém o vértice anterior no caminho
            // Reduz a capacidade da aresta no sentido de 'anterior' para 'atual'
            capacidade[anterior][atual] -= novo_fluxo;
            // Aumenta a capacidade da aresta no sentido oposto ('atual' para 'anterior')
            capacidade[atual][anterior] += novo_fluxo;
            atual = anterior;  // Move para o vértice anterior no caminho
        }
    }

    return fluxo;  // Retorna o fluxo máximo encontrado
}


void fecho() {
    // Inicializa o vetor 'visitado' com 0, indicando que nenhum vértice foi visitado
    visitado.assign(n_vertices, 0);
    
    // Executa a busca em profundidade (DFS) a partir do vértice 0
    // O segundo parâmetro 'true' pode indicar que a DFS deve imprimir os vértices visitados ou realizar outra ação específica
    dfs(0, true);
    
    // Imprime uma nova linha após a execução da DFS
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
    if (direcionado.compare("direcionado") == 0)
    {
        // Se o tipo de grafo é "direcionado"
        b_direcionado = DIRECIONADO; // Define a variável 'b_direcionado' para o valor correspondente a grafo direcionado
        leitura_direcionado();       // Chama a função para ler os dados de um grafo direcionado
    }
    else
    {
        // Se o tipo de grafo não é "direcionado" (assumindo que é "não direcionado")
        b_direcionado = NAO_DIRECIONADO; // Define a variável 'b_direcionado' para o valor correspondente a grafo não direcionado
        leitura_nao_direcionado();       // Chama a função para ler os dados de um grafo não direcionado
    }

    // Ordena os vértices vizinhos de cada vértice
    for (int i = 0; i < n_vertices; ++i)
    {
        sort(lista_adj[i].begin(), lista_adj[i].end()); // Ordena a lista de adjacências para o vértice 'i'
    }

    // Processa as funções especificadas
    for (int i = 0; i < funcoes.size(); ++i)
    {
        // cout << "f: " <<  funcoes[i] << endl; // Comentado para depuração futura, exibe a função atual sendo processada

        // As funções são indexadas a partir de 0
        switch (funcoes[i])
        {
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
