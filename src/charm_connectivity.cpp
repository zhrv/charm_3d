//
// Created by zhrv on 19.10.17.
//

#include "charm_globals.h"
#include <map>
#include <set>

using std::map;
using std::set;
using std::array;
using std::make_pair;

typedef p4est_connectivity_t* (*charm_conn_reader_t)(charm_ctx_t*);

const int charm_face_corners[6][4] =
        {{ 0, 2, 4, 6 },
         { 1, 3, 5, 7 },
         { 0, 1, 4, 5 },
         { 2, 3, 6, 7 },
         { 0, 1, 2, 3 },
         { 4, 5, 6, 7 }};

typedef struct charm_fhash_elem {
    p4est_topidx_t* key;
    sc_array_t*     val;
} charm_hash_elem_t;

typedef struct charm_fhash {
    p4est_topidx_t      size;
    charm_hash_elem_t*  el;
} charm_fhash_t;

typedef struct charm_face_info {
    p4est_topidx_t  tree;
    int             type;
} charm_face_info_t;

typedef struct {
    char  name[128];
    int   id;
    int   dim;
} charm_patch_t;

typedef array<p4est_topidx_t, CHARM_FACES>      charm_farr_t;
typedef set<p4est_topidx_t>                     charm_vset_t;
typedef map<charm_vset_t, p4est_topidx_t>       charm_cmap_t;
typedef map<charm_vset_t, charm_farr_t>         charm_fmap_t;

static sc_array_t *patches;





int charm_face_cmp(const void *_key1, const void *_key2)
{
    p4est_topidx_t key1[4];
    p4est_topidx_t key2[4];
    if (_key1 == NULL || _key2 == NULL) {
        return 0;
    }
    memcpy(key1, _key1, 4*sizeof(p4est_topidx_t));
    memcpy(key2, _key2, 4*sizeof(p4est_topidx_t));
    p4est_topidx_bsort(key1, 4);
    p4est_topidx_bsort(key2, 4);
    if (    key1[0] == key2[0] &&
            key1[1] == key2[1] &&
            key1[2] == key2[2] &&
            key1[3] == key2[3] ) {
        return 1;
    }
    return 0;
}

unsigned charm_face_hash(const void *v)
{
    const p4est_topidx_t *fi = (p4est_topidx_t *) v;

#ifdef P4_TO_P8
    return p4est_topidx_hash4 (fi);
#else
    return p4est_topidx_hash2 (fi);
#endif
}

charm_fhash_t* charm_fhash_new(p4est_topidx_t size)
{
    int i;
    charm_fhash_t * fh = CHARM_ALLOC(charm_fhash_t, 1);
    fh->size = size;
    fh->el = CHARM_ALLOC(charm_hash_elem_t, size);
    for (i = 0; i < size; i++) {
        fh->el[i].key = NULL;
        fh->el[i].val = sc_array_new(sizeof(charm_face_info_t));
    }
    return fh;
}

void charm_fhash_insert(charm_fhash_t* fh, const p4est_topidx_t* key, charm_face_info_t* val)
{
    charm_face_info_t * pval;
    unsigned idx = charm_face_hash(key) % fh->size;
    while(fh->el[idx].key != NULL) {
        if (charm_face_cmp(fh->el[idx].key, key)) {
            pval = (charm_face_info_t*)sc_array_push(fh->el[idx].val);
            memcpy(pval, val, sizeof(charm_face_info_t));
            return;
        }
        idx = (idx + 1) % fh->size;
    }
    fh->el[idx].key = CHARM_ALLOC(p4est_topidx_t, 4);
    memcpy(fh->el[idx].key, key, 4*sizeof(p4est_topidx_t));
    pval = (charm_face_info_t*)sc_array_push(fh->el[idx].val);
    memcpy(pval, val, sizeof(charm_face_info_t));
}


sc_array_t* charm_fhash_lookup(charm_fhash_t* fh, p4est_topidx_t* key)
{
    unsigned c = 0;
    unsigned idx = charm_face_hash(key) % fh->size;
    while(!charm_face_cmp(fh->el[idx].key, key)) {
        idx = (idx + 1) % fh->size;
        c++;
        if (c > fh->size) {
            return 0;
        }
    }
    return fh->el[idx].val;
}

/*
 * Read a line from a file. Obtained from:
 * http://stackoverflow.com/questions/314401/
 * how-to-read-a-line-from-the-console-in-c/314422#314422
 *
 * Using this avoids a dependence on IEEE Std 1003.1-2008 (``POSIX.1'') for the
 * getline function.
 */
static char        *
charm_connectivity_getline_upper (FILE * stream)
{
    char               *line = CHARM_ALLOC (char, 1024), *linep = line;
    size_t              lenmax = 1024, len = lenmax;
    int                 c;

    if (line == NULL)
        return NULL;

    for (;;) {
        c = fgetc (stream);
        c = toupper (c);
        if (c == EOF && linep == line) {
            CHARM_FREE (linep);
            return NULL;
        }

        if (--len == 0) {
            char               *linen;

            len = lenmax;
            lenmax *= 2;

            linen = CHARM_REALLOC (linep, char, lenmax);
            if (linen == NULL) {
                CHARM_FREE (linep);
                return NULL;
            }

            line = linen + (line - linep);
            linep = linen;
        }
        if ((*line++ = c) == '\n')
            break;
    }
    *line = '\0';
    return linep;
}


//int charm_conn_faces_is_equal(p4est_topidx_t  *face_vert_1, p4est_topidx_t  *face_vert_2)
//{
//    return charm_face_cmp(face_vert_1, face_vert_2);
//}


charm_bnd_t * charm_conn_bnd_find_by_id(charm_ctx_t* ctx, int id)
{
    size_t i;
    int patch_found = 0;
    charm_patch_t *p;
    charm_bnd_t   *bnd;
    char *name = NULL;
    // find by name
    for (i = 0; i < patches->elem_count; i++) {
        p = (charm_patch_t*)sc_array_index(patches, i);
        if (p->id == id) {
            patch_found = 1;
            name = p->name;
        }
    }



    CHARM_ASSERT(patch_found);

    for (i = 0; i < ctx->bnd->elem_count; i++) {
        bnd = (charm_bnd_t*)sc_array_index(ctx->bnd, i);
        if (strcmp(name, bnd->name) == 0) {
            return bnd;
        }
    }
    CHARM_LERRORF("Patch (id=%d, name='%s') is nod defined in task.yaml\n", id, name);
    return NULL;
}


charm_reg_t * charm_conn_reg_find_by_id(charm_ctx_t* ctx, int id)
{
    size_t i;
    int patch_found = 0;
    charm_patch_t *p;
    charm_reg_t   *reg;
    char *name = NULL;
    // find by name
    for (i = 0; i < patches->elem_count; i++) {
        p = (charm_patch_t*)sc_array_index(patches, i);
        if (p->id == id) {
            patch_found = 1;
            name = p->name;
        }
    }



    CHARM_ASSERT(patch_found);

    for (i = 0; i < ctx->reg->elem_count; i++) {
        reg = (charm_reg_t*)sc_array_index(ctx->reg, i);
        if (strcmp(name, reg->name) == 0) {
            return reg;
        }
    }
    CHARM_LERRORF("Region (id=%d, name='%s') is nod defined in task.yaml\n", id, name);
    return NULL;
}


void charm_conn_find_tree_by_face(charm_ctx_t *ctx, p4est_connectivity_t  *conn, charm_fhash_t* fh, p4est_topidx_t  *face_vert, int8_t face_type)
{
    int i;
    charm_tree_attr_t *attr;

    sc_array_t* arr = charm_fhash_lookup(fh, face_vert);
    CHARM_ASSERT(arr != NULL);
    for (i = 0; i < arr->elem_count; i++) {
        charm_face_info_t * fi = (charm_face_info_t *)sc_array_index(arr, i);
        attr = (charm_tree_attr_t *)&(conn->tree_to_attr[sizeof(charm_tree_attr_t)*fi->tree]);
        attr->bnd[fi->type] = charm_conn_bnd_find_by_id(ctx, face_type);
    }
}


int8_t charm_conn_parse_face(char *line, p4est_topidx_t *fv)
{
    int retval;
    int       node, type;

    retval = sscanf (line, "%d %d", &node, &type);
    if (retval != 2) {
        CHARM_LERROR ("Premature end of file");
        return 0;
    }

    if (type == 3) {
        int tagc, tag2;
        retval = sscanf (line, "%d %d %d %d %d %d %d %d %d",
                         &node, &type, &tagc, &fv[4], &tag2,
                         &(fv[0]), &(fv[1]), &(fv[3]), &(fv[2]) );
        if (retval != 9) {
            CHARM_LERROR ("Premature end of file");
            return 0;
        }
        fv[0]--;
        fv[1]--;
        fv[2]--;
        fv[3]--;
        return 1;
    }

    return 0;
}


int8_t charm_conn_parse_cell(char *line, p4est_topidx_t *fv)
{
    int i, retval;
    int       node, type;

    retval = sscanf (line, "%d %d", &node, &type);
    if (retval != 2) {
        CHARM_LERROR ("Premature end of file");
        return 0;
    }

    if (type == 5) {
        int tagc, tag2;
        retval = sscanf (line, "%d %d %d %d %d %d %d %d %d %d %d %d %d",
                         &node, &type, &tagc, &fv[8], &tag2,
                         &(fv[0]), &(fv[1]), &(fv[3]), &(fv[2]),
                         &(fv[4]), &(fv[5]), &(fv[7]), &(fv[6]));
        if (retval != 13) {
            CHARM_LERROR ("Premature end of file");
            return 0;
        }
        for (i = 0; i < 8; i++) {
            fv[i]--;
        }
        return 1;
    }

    return 0;
}


p4est_connectivity_t * charm_conn_reader_msh (charm_ctx_t* ctx)
{
    int                 retval;
    p4est_topidx_t      num_vertices = 0, num_patches = 0, num_trees = 0, num_trees_real = 0, tree;
    int                 i;
    int                 ele_count;
    int8_t              face;

    char                *line;

    long long int       node, type;
    charm_real_t              x, y, z;

    p4est_connectivity_t *conn = NULL;
    FILE                 *fid  = NULL;
    fpos_t                fpos;

    charm_fhash_t *     fh;
    charm_face_info_t   fi;
    p4est_topidx_t      fkey[4];

    charm_patch_t      *patch;
    charm_int_t         mpi_rank;

    charm_cmap_t        cmap;
    charm_fmap_t        fmap;
    charm_int_t         mpi_size;

    char *filename;

    CHARM_ASSERT(ctx->msh->type != CHARM_MESH_UNKNOWN);
    CHARM_ASSERT(ctx->msh->filename);

    filename = ctx->msh->filename;

    CHARM_GLOBAL_ESSENTIALF ("Reading connectivity from '%s'.\n", filename);

    fid = fopen (filename, "rb");
    if (fid == NULL) {
        CHARM_LERRORF ("Failed to open %s\n", filename);
        goto dead;
    }

    charm_real_t          *vert;
    p4est_topidx_t  *trees;

    while(1) {
        line = charm_connectivity_getline_upper (fid);

        if (line == NULL) {
            break;
        }

        if (line[0] == '$') {
            if (strstr(line, "$PHYSICALNAMES")) {
                CHARM_FREE(line);
                line = charm_connectivity_getline_upper (fid);
                sscanf(line, "%d", &num_patches);
                patches = sc_array_new(sizeof(charm_patch_t));
                for (i = 0; i < num_patches; i++) {
                    CHARM_FREE(line);
                    line = charm_connectivity_getline_upper (fid);
                    patch = (charm_patch_t*)sc_array_push(patches);
                    sscanf(line, "%d %d \"%[^\"]", &(patch->dim), &(patch->id), patch->name);
                }
            }
            else if (strstr(line, "$NODES")) {
                CHARM_FREE(line);
                line = charm_connectivity_getline_upper (fid);
                sscanf(line, "%d", &num_vertices);
                vert = CHARM_ALLOC(charm_real_t, num_vertices*3);
                for (i = 0; i < num_vertices; i++) {
                    CHARM_FREE(line);
                    line = charm_connectivity_getline_upper (fid);
                    retval = sscanf (line, "%lld %lf %lf %lf", &node, &x, &y, &z);
                    if (retval != 4) {
                        CHARM_LERROR ("Premature end of file");
                        CHARM_FREE (line);
                        CHARM_FREE (vert);
                        return NULL;
                    }
                    vert[3 * i + 0] = x;
                    vert[3 * i + 1] = y;
                    vert[3 * i + 2] = z;
                }
            }
            else if(strstr(line, "$ELEMENTS")) {
                fgetpos (fid, &fpos);
                CHARM_FREE(line);
                line = charm_connectivity_getline_upper (fid);
                sscanf(line, "%d", &num_trees);
                trees = CHARM_ALLOC(p4est_topidx_t, num_trees*8);
                fh = charm_fhash_new(num_trees*6);
                for (i = 0; i < num_trees; i++) {
                    CHARM_FREE(line);
                    line = charm_connectivity_getline_upper (fid);
                    retval = sscanf (line, "%lld %lld", &node, &type);
                    if (retval != 2) {
                        CHARM_LERROR ("Premature end of file");
                        CHARM_FREE (line);
                        CHARM_FREE (vert);
                        CHARM_FREE (trees);
                        return NULL;
                    }

                    if (type == 5) {
                        int tagc, tag1, tag2, j, k;
                        retval = sscanf (line, "%lld %lld %d %d %d %d %d %d %d %d %d %d %d", &node, &type, &tagc, &tag1, &tag2,
                                         &trees[num_trees_real*CHARM_CHILDREN+0], &trees[num_trees_real*CHARM_CHILDREN+1], &trees[num_trees_real*CHARM_CHILDREN+3], &trees[num_trees_real*CHARM_CHILDREN+2],
                                         &trees[num_trees_real*CHARM_CHILDREN+4], &trees[num_trees_real*CHARM_CHILDREN+5], &trees[num_trees_real*CHARM_CHILDREN+7], &trees[num_trees_real*CHARM_CHILDREN+6] );
                        if (retval != 13) {
                            CHARM_LERROR ("Premature end of file");
                            CHARM_FREE (line);
                            CHARM_FREE (vert);
                            CHARM_FREE (trees);
                            return NULL;
                        }
                        for (j = 0; j < CHARM_CHILDREN; j++) {
                            --trees[num_trees_real*CHARM_CHILDREN+j];
                        }
                        for (j = 0; j < CHARM_FACES; j++) {
                            fi.type = j;
                            fi.tree = num_trees_real;
                            for (k = 0; k < CHARM_HALF; k++) {
                                fkey[k] = trees[num_trees_real*CHARM_CHILDREN+charm_face_corners[j][k]];
                            }
                            charm_fhash_insert(fh, fkey, &fi);
                        }
                        ++num_trees_real;
                    }
                }
            }
        }

        CHARM_FREE(line);
    }

    conn = p4est_connectivity_new(num_vertices, num_trees_real, 0, 0, 0, 0);
    memcpy(conn->vertices, vert, conn->num_vertices*3*sizeof(charm_real_t));
    memcpy(conn->tree_to_vertex, trees, conn->num_trees*CHARM_CHILDREN*sizeof(p4est_topidx_t));

    CHARM_FREE(vert);
    CHARM_FREE(trees);

    /*
     * Fill tree_to_tree and tree_to_face to make sure we have a valid
     * connectivity.
     */
    for (tree = 0; tree < conn->num_trees; ++tree) {
        for (face = 0; face < CHARM_FACES; ++face) {
            conn->tree_to_tree[CHARM_FACES * tree + face] = tree;
            conn->tree_to_face[CHARM_FACES * tree + face] = face;
        }
    }
    CHARM_ASSERT (p4est_connectivity_is_valid (conn));


    /* Compute real tree_to_* fields and complete (edge and) corner fields. */
    p4est_connectivity_complete (conn);

#ifdef P4EST_WITH_METIS
    sc_MPI_Comm_size(sc_MPI_COMM_WORLD, &mpi_size);
    if (mpi_size > 1)
        p4est_connectivity_reorder (sc_MPI_COMM_WORLD, 0, conn, P4EST_CONNECT_FACE);
#endif /* P4EST_WITH_METIS */



    p4est_connectivity_set_attr(conn, sizeof(charm_tree_attr_t));
    charm_tree_attr_t * attr;
    for (tree = 0, attr = (charm_tree_attr_t*)(conn->tree_to_attr); tree < conn->num_trees; ++tree, ++attr) {
        for (face = 0; face < CHARM_FACES; ++face) {
            attr->bnd[face] = NULL;
        }
        attr->reg = NULL;
        set<charm_int_t> sv;
        for (int vert = 0; vert < CHARM_CHILDREN; vert++) {
            sv.insert(conn->tree_to_vertex[tree*CHARM_CHILDREN+vert]);
        }
        cmap.insert(make_pair(sv, tree));

        charm_farr_t farr;
        farr.fill(-1);

        for (int j = 0; j < CHARM_FACES; j++) {
            sv.clear();
            for (int k = 0; k < CHARM_HALF; k++) {
                sv.insert(conn->tree_to_vertex[tree*CHARM_CHILDREN+charm_face_corners[j][k]]);
            }
            if (fmap.find(sv) == fmap.end()) {
                fmap[sv] = farr;
            }
            fmap[sv][j] = tree;
        }
    }

    fsetpos (fid, &fpos);
    line = charm_connectivity_getline_upper(fid);
    ele_count = 0;
    sscanf(line, "%d", &ele_count);
    CHARM_FREE(line);
    p4est_topidx_t fv[32];
    attr = (charm_tree_attr_t*)(conn->tree_to_attr);
    for (; ele_count > 0; ele_count--) {
        line = charm_connectivity_getline_upper(fid);
        if (charm_conn_parse_face(line, fv)) {
            //charm_conn_find_tree_by_face(ctx, conn, fh, fv, fv[4]);
            set<charm_int_t> sv;
            for (int v = 0; v < CHARM_HALF; v++) {
                sv.insert(fv[v]);
            }
            charm_farr_t &fa = fmap[sv];
            for (int f = 0; f < CHARM_FACES; f++) {
                if (fa[f] > -1) {
                    attr[fa[f]].bnd[f] = charm_conn_bnd_find_by_id(ctx, fv[4]);
                }
            }
        }
        else if (charm_conn_parse_cell(line, fv)) { // @todo предполагаем, что последовательность считывания совпадет с предыдущей последовательностью
            //charm_conn_find_tree_by_cell(ctx, conn, ch, fv, fv[4]);
            set<charm_int_t> sv;
            for (int v = 0; v < CHARM_CHILDREN; v++) {
                sv.insert(fv[v]);
            }
            i = cmap[sv];
            attr[i].reg = charm_conn_reg_find_by_id(ctx, fv[8]);
            //attr++;
        }
        CHARM_FREE(line);
    }

    sc_array_destroy(patches);
    retval = fclose (fid);
    fid = NULL;
    if (retval) {
        CHARM_LERRORF ("Failed to close %s\n", filename);
        goto dead;
    }

    cmap.clear();
    fmap.clear();

    CHARM_GLOBAL_ESSENTIALF
    ("New connectivity with %lld trees and %lld vertices\n",
     (long long) conn->num_trees, (long long) conn->num_vertices);

    return conn;

    dead:
    /* clean up on error */
    if (fid != NULL) {
        fclose (fid);
    }
    if (conn != NULL) {
        p4est_connectivity_destroy (conn);
    }
    return NULL;
}





/**
 * INP-file reader
 */

int
charm_connectivity_read_inp_stream (FILE * stream,
                                    p4est_topidx_t * num_vertices,
                                    p4est_topidx_t * num_trees,
                                    charm_real_t *vertices,
                                    p4est_topidx_t * tree_to_vertex)
{
    int                 reading_nodes = 0, reading_elements = 0;
    int                 lines_read = 0, lines_free = 0;
    char               *line;
    p4est_topidx_t      num_nodes = 0;
    p4est_topidx_t      num_elements = 0;
    int                 fill_trees_and_vertices = (vertices != NULL &&
                                                   tree_to_vertex != NULL);

    CHARM_ASSERT ((vertices == NULL && tree_to_vertex == NULL) ||
                  (vertices != NULL && tree_to_vertex != NULL));

    for (;;) {
        line = charm_connectivity_getline_upper (stream);

        if (line == NULL) {
            break;
        }

        ++lines_read;

        /* check for control line */
        if (line[0] == '*') {
            reading_elements = reading_nodes = 0;
            if (strstr (line, "*NODE")) {
                reading_nodes = 1;
                ++lines_free;
                CHARM_FREE (line);
                continue;
            }
            else if (strstr (line, "*ELEMENT")) {
                if (
#ifdef P4_TO_P8
                        strstr (line, "TYPE=C3D8")
#else
                    strstr (line, "TYPE=C2D4") || strstr (line, "TYPE=CPS4")
             || strstr (line, "TYPE=S4")
#endif
                        ) {
                    reading_elements = 1;
                    ++lines_free;
                    CHARM_FREE (line);
                    continue;
                }
            }
        }

        if (reading_nodes) {
            if (fill_trees_and_vertices) {
                long long int       node;
                charm_real_t              x, y, z;
                int                 retval;

                retval = sscanf (line, "%lld, %lf, %lf, %lf", &node, &x, &y, &z);
                if (retval != 4) {
                    CHARM_LERROR ("Premature end of file");
                    CHARM_FREE (line);
                    return 1;
                }

                if (node > *num_vertices) {
                    CHARM_LERRORF
                    ("Encountered vertex %lld that will not fit in vertices"
                             " array of length %lld.  Are the vertices contiguously"
                             " numbered?\n", node, (long long int) *num_vertices);
                    CHARM_FREE (line);
                    return 1;
                }

                vertices[3 * (node - 1) + 0] = x;
                vertices[3 * (node - 1) + 1] = y;
                vertices[3 * (node - 1) + 2] = z;
            }

            ++num_nodes;
        }
        else if (reading_elements) {
            if (fill_trees_and_vertices) {
                long long int       element_number;
                long long int       v[CHARM_CHILDREN];
                int                 n;
                int                 retval;

                /* Note that when we read in the
                 * vertices we switch from right-hand
                 * vertex ordering to z-order
                 */
                retval = sscanf (line, "%lld, %lld, %lld, %lld, %lld"
#ifdef P4_TO_P8
                                         ", %lld, %lld, %lld, %lld"
#endif
                        , &element_number, &v[0], &v[1], &v[3], &v[2]
#ifdef P4_TO_P8
                        , &v[4], &v[5], &v[7], &v[6]
#endif
                );
                if (retval != CHARM_CHILDREN + 1) {
                    CHARM_LERROR ("Premature end of file");
                    CHARM_FREE (line);
                    return 1;
                }

                if (num_elements >= *num_trees) {
                    CHARM_LERROR ("Encountered element that will not fit into"
                                          " tree_to_vertex array. More elements than expected.\n");
                    CHARM_FREE (line);
                    return 1;
                }

                for (n = 0; n < CHARM_CHILDREN; ++n)
                    tree_to_vertex[CHARM_CHILDREN * num_elements + n] =
                            v[n] - 1;
            }

            ++num_elements;
        }

        ++lines_free;
        CHARM_FREE (line);
    }

    *num_vertices = num_nodes;
    *num_trees = num_elements;

    if (num_nodes == 0 || num_elements == 0) {
        CHARM_LERROR ("No elements or nodes found in mesh file.\n");
        return -1;
    }
    else {
        return 0;
    }
}


p4est_connectivity_t *
charm_conn_reader_inp (charm_ctx_t *ctx)
{
    int                 retval;
    p4est_topidx_t      num_vertices = 0, num_trees = 0, tree;
    int                 face;

    p4est_connectivity_t *conn = NULL;
    FILE               *fid = NULL;

    char *filename;

    CHARM_ASSERT(ctx->msh->type != CHARM_MESH_UNKNOWN);
    CHARM_ASSERT(ctx->msh->filename);

    filename = ctx->msh->filename;

    CHARM_GLOBAL_ESSENTIALF ("Reading connectivity from '%s'.\n", filename);

    fid = fopen (filename, "rb");
    if (fid == NULL) {
        CHARM_LERRORF ("Failed to open %s\n", filename);
        goto dead;
    }

    if (charm_connectivity_read_inp_stream
            (fid, &num_vertices, &num_trees, NULL, NULL)) {
        CHARM_LERRORF ("Failed to read %s: pass 1\n", filename);
        goto dead;
    }

    rewind (fid);

    conn = p4est_connectivity_new (num_vertices, num_trees,
#ifdef P4_TO_P8
                                   0, 0,
#endif
                                   0, 0);

    if (charm_connectivity_read_inp_stream (fid, &conn->num_vertices,
                                            &conn->num_trees, conn->vertices,
                                            conn->tree_to_vertex)) {
        CHARM_LERRORF ("Failed to read %s: pass 2\n", filename);
        goto dead;
    }

    /*
     * Fill tree_to_tree and tree_to_face to make sure we have a valid
     * connectivity.
     */
    for (tree = 0; tree < conn->num_trees; ++tree) {
        for (face = 0; face < CHARM_FACES; ++face) {
            conn->tree_to_tree[CHARM_FACES * tree + face] = tree;
            conn->tree_to_face[CHARM_FACES * tree + face] = face;
        }
    }
    CHARM_ASSERT (p4est_connectivity_is_valid (conn));

    /* Compute real tree_to_* fields and complete (edge and) corner fields. */
    p4est_connectivity_complete (conn);

    retval = fclose (fid);
    fid = NULL;
    if (retval) {
        CHARM_LERRORF ("Failed to close %s\n", filename);
        goto dead;
    }

    CHARM_GLOBAL_ESSENTIALF
    ("New connectivity with %lld trees and %lld vertices\n",
     (long long) conn->num_trees, (long long) conn->num_vertices);

    return conn;

    dead:
    /* clean up on error */
    if (fid != NULL) {
        fclose (fid);
    }
    if (conn != NULL) {
        p4est_connectivity_destroy (conn);
    }
    return NULL;
}


/**
 *
 */

charm_conn_reader_t charm_get_conn_reader(charm_ctx_t *ctx)
{
    switch (ctx->msh->type) {
        case CHARM_MESH_GMSH_MSH:
            return charm_conn_reader_msh;
        case CHARM_MESH_GMSH_INP:
            return charm_conn_reader_inp;
        case CHARM_MESH_UNKNOWN:
        default:
            return NULL;
    }
}


p4est_connectivity_t* charm_conn_create(charm_ctx_t *ctx)
{
    charm_conn_reader_t     reader_fn;
    p4est_connectivity_t*   conn;
    const p4est_topidx_t    num_vertices    = 44;
    const p4est_topidx_t    num_trees       = 10;

    reader_fn = charm_get_conn_reader(ctx);

    conn = reader_fn(ctx);


    return conn;
}

