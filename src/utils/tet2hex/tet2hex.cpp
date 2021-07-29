#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <cstring>

enum elem_type_t {
    ELEMENT_THREE_NODE_TRIANGLE = 2,
    ELEMENT_FOUR_NODE_QUADRANGLE = 3,
    ELEMENT_FOUR_NODE_TETRAHEDRON = 4,
    ELEMENT_EIGHT_NODE_HEXAHEDRON = 5
};

inline void sort(int *arr, int n) {
    for (int i = 0; i < n; i++) {
        for (int j = 1; j < n-i; j++) {
            if (arr[j-1] > arr[j]) {
                int tmp = arr[j-1];
                arr[j-1] = arr[j];
                arr[j] = tmp;
            }
        }
    }
}


struct patch_t {
    int id;
    int dim;
    std::string name;
};


struct node_t {
    double x, y, z;
};

struct tface_t {
    int nodes[3];
    int tag;
    int c;
};

struct tcell_t {
    int nodes[4];
    int tag;
    int c;
};

struct tedge_t {
    int nodes[2];
    int c;
};

struct qface_t {
    int nodes[4];
    int tag;
};

struct qcell_t {
    int nodes[8];
    int tag;
};



typedef std::map<int, patch_t>  patches_t;
typedef std::vector<node_t>     nodes_t;
typedef std::vector<tface_t>    tfaces_t;
typedef std::vector<tcell_t>    tcells_t;
typedef std::vector<tedge_t>    tedges_t;
typedef std::vector<qface_t>    qfaces_t;
typedef std::vector<qcell_t>    qcells_t;

typedef std::tuple<int, int, int> tri_t;
typedef std::tuple<int, int> edg_t;

bool operator<(const tri_t t1,  const tri_t t2) {
    return std::tie(std::get<0>(t1), std::get<1>(t1), std::get<2>(t1)) < std::tie(std::get<0>(t2), std::get<1>(t2), std::get<2>(t2));
}

std::map<tri_t, int> tri2index;
std::map<edg_t, int> edg2index;


tri_t make_tri(int *nodes) {
    sort(nodes, 3);
    return std::make_tuple(nodes[0], nodes[1], nodes[2]);
}

tri_t make_tri(int n1, int n2, int n3) {
    int nodes[] = {n1, n2, n3};
    return make_tri(nodes);
}

edg_t make_edg(int *nodes) {
    sort(nodes, 2);
    return std::make_tuple(nodes[0], nodes[1]);
}

edg_t make_edg(int n1, int n2) {
    int nodes[] = {n1, n2};
    return make_edg(nodes);
}

inline void get_line_upper(std::ifstream &fin, std::string &line) {
    std::getline(fin, line);
    for (auto &c: line) {
        c = toupper(c);
    }
}


patches_t patches;
nodes_t nodes;
tcells_t tcells;
tfaces_t tfaces;
tedges_t tedges;

qcells_t qcells;
qfaces_t qfaces;


tedge_t& get_edge(int n1, int n2) {
    edg_t edg = make_edg(n1, n2);
    if (edg2index.find(edg) == edg2index.end()) {
        tedge_t edge{};
        edge.nodes[0] = std::get<0>(edg);
        edge.nodes[1] = std::get<1>(edg);
        node_t &node1 = nodes[n1];
        node_t &node2 = nodes[n2];
        node_t center{};
        center.x = (node1.x+node2.x)/2.;
        center.y = (node1.y+node2.y)/2.;
        center.z = (node1.z+node2.z)/2.;
        nodes.push_back(center);
        edge.c = nodes.size()-1;
        tedges.push_back(edge);
        edg2index[edg] = tedges.size() - 1;
    }

    return tedges[edg2index[edg]];
}

tface_t& get_face(int n1, int n2, int n3) {
    tri_t tri = make_tri(n1, n2, n3);
    if (tri2index.find(tri) == tri2index.end()) {
        tface_t f{};
        f.nodes[0] = std::get<0>(tri);
        f.nodes[1] = std::get<1>(tri);
        f.nodes[2] = std::get<2>(tri);
        node_t &node0 = nodes[f.nodes[0]];
        node_t &node1 = nodes[f.nodes[1]];
        node_t &node2 = nodes[f.nodes[2]];
        node_t center{};
        center.x = (node0.x+node1.x+node2.x)/3.;
        center.y = (node0.y+node1.y+node2.y)/3.;
        center.z = (node0.z+node1.z+node2.z)/3.;
        f.c = nodes.size();
        f.tag = -1;
        tfaces.push_back(f);
        nodes.push_back(center);
        tri2index[tri] = tfaces.size() - 1;
    }
    return tfaces[tri2index[make_tri(n1, n2, n3)]];
}

int main(int argc, char* argv[]) {
    if (argc != 2) {
        std::cout << "Usage: msh-tet2hex <file_name>" << std::endl;
        return 1;
    }

    std::string sfin(argv[1]);
    std::ifstream fin(sfin);

    std::string line;

    while (!fin.eof()) {
        get_line_upper(fin, line);
        if (line[0] == '$') {
            if (line.find("$PHYSICALNAMES") != std::string::npos) {
                int n, dim, id;
                std::string name;
                fin >> n;
                for (int i = 0; i < n; i++) {
                    fin >> dim >> id >> name;
                    patch_t p;
                    p.id = id;
                    p.dim = dim;
                    p.name = name;
                    patches[id] = p;
                }
            }
            else if (line.find("$NODES") != std::string::npos) {
                int n, id;
                node_t node;
                fin >> n;
                for (int i = 0; i < n; i++) {
                    fin >> id >> node.x >> node.y >> node.z;
                    nodes.push_back(node);
                }

            }
            else if (line.find("$ELEMENTS") != std::string::npos) {
                /*
                $Elements
                number-of-elements
                elm-number elm-type number-of-tags < tag > ... node-number-list
                ...
                $EndElements
                */
                int n, id, type, nt, tag, tmp;
                fin >> n;
                for (int i = 0; i < n; i++) {
                    fin >> id >> type >> nt >> tag;
                    for (int j = 1; j < nt; j++) {
                        fin >> tmp;
                    }
                    if (type == ELEMENT_THREE_NODE_TRIANGLE) {
                        tface_t f{};
                        f.tag = tag;
                        fin >> f.nodes[0] >> f.nodes[1] >> f.nodes[2];
                        f.nodes[0]--;
                        f.nodes[1]--;
                        f.nodes[2]--;
                        sort(f.nodes, 3);
                        node_t &n0 = nodes[f.nodes[0]];
                        node_t &n1 = nodes[f.nodes[1]];
                        node_t &n2 = nodes[f.nodes[2]];
                        node_t c;
                        c.x = (n0.x+n1.x+n2.x)/3.;
                        c.y = (n0.y+n1.y+n2.y)/3.;
                        c.z = (n0.z+n1.z+n2.z)/3.;
                        f.c = nodes.size();
                        nodes.push_back(c);
                        tfaces.push_back(f);
                        tri2index[make_tri(f.nodes)] = tfaces.size() - 1;
                    }
                    else if (type == ELEMENT_FOUR_NODE_TETRAHEDRON) {
                        tcell_t c{};
                        c.tag = tag;
                        fin >> c.nodes[0] >> c.nodes[1] >> c.nodes[2] >> c.nodes[3];
                        c.nodes[0]--; c.nodes[1]--; c.nodes[2]--; c.nodes[3]--;
                        node_t &n0 = nodes[c.nodes[0]];
                        node_t &n1 = nodes[c.nodes[1]];
                        node_t &n2 = nodes[c.nodes[2]];
                        node_t &n3 = nodes[c.nodes[3]];
                        node_t cen;
                        cen.x = (n0.x+n1.x+n2.x+n3.x)/4.;
                        cen.y = (n0.y+n1.y+n2.y+n3.y)/4.;
                        cen.z = (n0.z+n1.z+n2.z+n3.z)/4.;
                        c.c = nodes.size();
                        nodes.push_back(cen);
                        tcells.push_back(c);
                    }
                    else {
                        std::cerr << "Wrong element type: " << type << ". Id: " << id << std::endl;
                        return 1;
                    }
                }
            }
        }
    }


    for (auto face: tfaces) {
        if (face.tag == -1) continue;
        int n[3];
        for (int iv = 0; iv < 3; iv++) {
            for (int j = 0; j < 3; j++) {
                n[j] = face.nodes[(iv + j) % 3];
            }
            tedge_t &e1 = get_edge(n[0], n[1]);
            tedge_t &e2 = get_edge(n[0], n[2]);
            qface_t qface{};
            qface.nodes[0] = n[0];
            qface.nodes[1] = e1.c;
            qface.nodes[2] = face.c;
            qface.nodes[3] = e2.c;
            qface.tag = face.tag;
            qfaces.push_back(qface);
        }
    }


    for (auto cell: tcells) {
        int n[4];
        for (int iv = 0; iv < 4; iv++) {
            for (int j = 0; j < 4; j++) {
                n[j] = cell.nodes[(iv + j) % 4];
            }
            tface_t &f_abc = get_face(n[0], n[1], n[3]);
            tface_t &f_acd = get_face(n[0], n[1], n[2]);
            tface_t &f_abd = get_face(n[0], n[2], n[3]);

            tedge_t &e_ab = get_edge(n[0], n[3]);
            tedge_t &e_ac = get_edge(n[0], n[1]);
            tedge_t &e_ad = get_edge(n[0], n[2]);

            int qn[8] = {
                    e_ad.c, f_abd.c, cell.c, f_acd.c,
                    n[0],  e_ab.c,  f_abc.c, e_ac.c
            };

            qcell_t qcell{};
            memcpy(qcell.nodes, qn, sizeof(int)*8);
            qcell.tag = cell.tag;
            qcells.push_back(qcell);
        }
    }


    std::string sfout(sfin);
    int pos = sfout.find(".msh");
    sfout.erase(pos, 4);
    sfout += ".quad.msh";
    std::ofstream fout(sfout);
    fout.precision(16);
    fout << "$MeshFormat" << std::endl;
    fout << "2.2 0 8" << std::endl;
    fout << "$EndMeshFormat" << std::endl;
    fout << "$PhysicalNames" << std::endl;
    fout << patches.size() << std::endl;
    for (auto p: patches) {
        fout << p.second.dim << ' ' << p.second.id << ' ' << p.second.name << std::endl;
    }
    fout << "$EndPhysicalNames" << std::endl;
    fout << "$Nodes" << std::endl;
    fout << nodes.size() << std::endl;
    int id = 0;
    for (auto n: nodes) {
        fout << ++id << ' ' << n.x << ' ' << n.y << ' ' << n.z << std::endl;
    }
    fout << "$EndNodes" << std::endl;
    fout << "$Elements" << std::endl;
    fout << qcells.size()+qfaces.size() << std::endl;
    id = 0;
    for (auto e: qfaces) {
        fout << ++id << ' ' << ELEMENT_FOUR_NODE_QUADRANGLE << ' ' << 1 << ' ' << e.tag << ' '
             << e.nodes[0]+1 << ' ' << e.nodes[1]+1 << ' ' << e.nodes[2]+1 << ' ' << e.nodes[3]+1 << std::endl;
    }
    for (auto e: qcells) {
        fout << ++id << ' ' << ELEMENT_EIGHT_NODE_HEXAHEDRON << ' ' << 1 << ' ' << e.tag << ' '
                << e.nodes[0]+1 << ' ' << e.nodes[1]+1 << ' ' << e.nodes[2]+1 << ' ' << e.nodes[3]+1 << ' '
                << e.nodes[4]+1 << ' ' << e.nodes[5]+1 << ' ' << e.nodes[6]+1 << ' ' << e.nodes[7]+1 << std::endl;
    }
    fout << "$EndElements" << std::endl;

    fout.close();

    return 0;
}
