//
// Created by zhrv on 19.10.17.
//

#include "charm_connectivity.h"

typedef p4est_connectivity_t* (*charm_conn_reader_t)(const char*);

const int charm_face_corners[6][4] =
        {{ 0, 2, 4, 6 },
         { 1, 3, 5, 7 },
         { 0, 1, 4, 5 },
         { 2, 3, 6, 7 },
         { 0, 1, 2, 3 },
         { 4, 5, 6, 7 }};


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
    char               *line = P4EST_ALLOC (char, 1024), *linep = line;
    size_t              lenmax = 1024, len = lenmax;
    int                 c;

    if (line == NULL)
        return NULL;

    for (;;) {
        c = fgetc (stream);
        c = toupper (c);
        if (c == EOF && linep == line) {
            P4EST_FREE (linep);
            return NULL;
        }

        if (--len == 0) {
            char               *linen;

            len = lenmax;
            lenmax *= 2;

            linen = P4EST_REALLOC (linep, char, lenmax);
            if (linen == NULL) {
                P4EST_FREE (linep);
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


int
charm_connectivity_read_inp_stream (FILE * stream,
                                    p4est_topidx_t * num_vertices,
                                    p4est_topidx_t * num_trees,
                                    double *vertices,
                                    p4est_topidx_t * tree_to_vertex)
{
    int                 reading_nodes = 0, reading_elements = 0;
    int                 lines_read = 0, lines_free = 0;
    char               *line;
    p4est_topidx_t      num_nodes = 0;
    p4est_topidx_t      num_elements = 0;
    int                 fill_trees_and_vertices = (vertices != NULL &&
                                                   tree_to_vertex != NULL);

    P4EST_ASSERT ((vertices == NULL && tree_to_vertex == NULL) ||
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
                P4EST_FREE (line);
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
                    P4EST_FREE (line);
                    continue;
                }
            }
        }

        if (reading_nodes) {
            if (fill_trees_and_vertices) {
                long long int       node;
                double              x, y, z;
                int                 retval;

                retval = sscanf (line, "%lld, %lf, %lf, %lf", &node, &x, &y, &z);
                if (retval != 4) {
                    P4EST_LERROR ("Premature end of file");
                    P4EST_FREE (line);
                    return 1;
                }

                if (node > *num_vertices) {
                    P4EST_LERRORF
                    ("Encountered vertex %lld that will not fit in vertices"
                             " array of length %lld.  Are the vertices contiguously"
                             " numbered?\n", node, (long long int) *num_vertices);
                    P4EST_FREE (line);
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
                long long int       v[P4EST_CHILDREN];
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
                if (retval != P4EST_CHILDREN + 1) {
                    P4EST_LERROR ("Premature end of file");
                    P4EST_FREE (line);
                    return 1;
                }

                if (num_elements >= *num_trees) {
                    P4EST_LERROR ("Encountered element that will not fit into"
                                          " tree_to_vertex array. More elements than expected.\n");
                    P4EST_FREE (line);
                    return 1;
                }

                for (n = 0; n < P4EST_CHILDREN; ++n)
                    tree_to_vertex[P4EST_CHILDREN * num_elements + n] =
                            v[n] - 1;
            }

            ++num_elements;
        }

        ++lines_free;
        P4EST_FREE (line);
    }

    *num_vertices = num_nodes;
    *num_trees = num_elements;

    if (num_nodes == 0 || num_elements == 0) {
        P4EST_LERROR ("No elements or nodes found in mesh file.\n");
        return -1;
    }
    else {
        return 0;
    }
}


p4est_connectivity_t *
charm_conn_reader_inp (const char *filename)
{
    int                 retval;
    p4est_topidx_t      num_vertices = 0, num_trees = 0, tree;
    int                 face;

    p4est_connectivity_t *conn = NULL;
    FILE               *fid = NULL;

    P4EST_GLOBAL_PRODUCTIONF ("Reading connectivity from %s\n", filename);

    fid = fopen (filename, "rb");
    if (fid == NULL) {
        P4EST_LERRORF ("Failed to open %s\n", filename);
        goto dead;
    }

    if (charm_connectivity_read_inp_stream
            (fid, &num_vertices, &num_trees, NULL, NULL)) {
        P4EST_LERRORF ("Failed to read %s: pass 1\n", filename);
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
        P4EST_LERRORF ("Failed to read %s: pass 2\n", filename);
        goto dead;
    }

    /*
     * Fill tree_to_tree and tree_to_face to make sure we have a valid
     * connectivity.
     */
    for (tree = 0; tree < conn->num_trees; ++tree) {
        for (face = 0; face < P4EST_FACES; ++face) {
            conn->tree_to_tree[P4EST_FACES * tree + face] = tree;
            conn->tree_to_face[P4EST_FACES * tree + face] = face;
        }
    }
    P4EST_ASSERT (p4est_connectivity_is_valid (conn));

    /* Compute real tree_to_* fields and complete (edge and) corner fields. */
    p4est_connectivity_complete (conn);

    retval = fclose (fid);
    fid = NULL;
    if (retval) {
        P4EST_LERRORF ("Failed to close %s\n", filename);
        goto dead;
    }

    P4EST_GLOBAL_PRODUCTIONF
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


int8_t charm_conn_faces_is_equal(p4est_topidx_t  *face_vert_1, p4est_topidx_t  *face_vert_2)
{
    if (    (face_vert_1[0] == face_vert_2[0]) && (face_vert_1[1] == face_vert_2[1]) &&
            (face_vert_1[2] == face_vert_2[2]) && (face_vert_1[3] == face_vert_2[3]) ) {
        return 1;
    }
    else {
        return 0;
    }
}

int charm_conn_find_tree_by_face(p4est_connectivity_t  *conn, p4est_topidx_t  *face_vert, int8_t bnd_type)
{
    int i, j, k;
    p4est_topidx_t  tree_faces[4];
    charm_tree_attr_t *attr;
    p4est_topidx_t *ttv;

    int find_count = 0;
    p4est_topidx_bsort(face_vert, 4);
    for (i = 0, attr = (charm_tree_attr_t*)(conn->tree_to_attr), ttv = conn->tree_to_vertex; i < conn->num_trees; i++, ++attr, ttv+=P4EST_CHILDREN) {
        for (j = 0; j < P4EST_FACES; j++) {
            for (k = 0; k < 4; k++) {
                tree_faces[k] = ttv[charm_face_corners[j][k]];
            }
            p4est_topidx_bsort (tree_faces, 4);
            if (charm_conn_faces_is_equal(tree_faces, face_vert)) {
                attr->bnd_type[j] = bnd_type;
                find_count++;
            }
        }
    }
    P4EST_ASSERT(find_count);
}

int8_t charm_conn_parse_face(char *line, p4est_topidx_t *fv)
{
    int retval;
    int       node, type;

    retval = sscanf (line, "%lld %lld", &node, &type);
    if (retval != 2) {
        P4EST_LERROR ("Premature end of file");
        return 0;
    }

    if (type == 3) {
        int tagc, tag1, tag2;
        retval = sscanf (line, "%d %d %d %d %d %d %d %d %d %d %d %d %d",
                         &node, &type, &tagc, &fv[4], &tag2,
                         &fv[0], &fv[1], &fv[3], &fv[2] );
        if (retval != 9) {
            P4EST_LERROR ("Premature end of file");
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

p4est_connectivity_t * charm_conn_reader_msh (const char *filename)
{
    int                 retval;
    p4est_topidx_t      num_vertices = 0, num_trees = 0, num_trees_real = 0, tree;
    int                 i;
    int8_t              face;

    char                *line;

    long long int       node, type;
    double              x, y, z;

    p4est_connectivity_t *conn = NULL;
    FILE                 *fid  = NULL;
    fpos_t                fpos;

    P4EST_GLOBAL_PRODUCTIONF ("Reading connectivity from %s\n", filename);

    fid = fopen (filename, "rb");
    if (fid == NULL) {
        P4EST_LERRORF ("Failed to open %s\n", filename);
        goto dead;
    }

    double          *vert;
    p4est_topidx_t  *trees;

    while(1) {
        line = charm_connectivity_getline_upper (fid);

        if (line == NULL) {
            break;
        }

        if (line[0] == '$') {
            if (strstr(line, "$NODES")) {
                P4EST_FREE(line);
                line = charm_connectivity_getline_upper (fid);
                sscanf(line, "%d", &num_vertices);
                vert = P4EST_ALLOC(double, num_vertices*3);
                for (i = 0; i < num_vertices; i++) {
                    P4EST_FREE(line);
                    line = charm_connectivity_getline_upper (fid);
                    retval = sscanf (line, "%lld %lf %lf %lf", &node, &x, &y, &z);
                    if (retval != 4) {
                        P4EST_LERROR ("Premature end of file");
                        P4EST_FREE (line);
                        P4EST_FREE (vert);
                        return NULL;
                    }
                    vert[3 * i + 0] = x;
                    vert[3 * i + 1] = y;
                    vert[3 * i + 2] = z;
                }
            }
            else if(strstr(line, "$ELEMENTS")) {
                fgetpos (fid, &fpos);
                P4EST_FREE(line);
                line = charm_connectivity_getline_upper (fid);
                sscanf(line, "%d", &num_trees);
                trees = P4EST_ALLOC(p4est_topidx_t, num_trees*8);
                for (i = 0; i < num_trees; i++) {
                    P4EST_FREE(line);
                    line = charm_connectivity_getline_upper (fid);
                    retval = sscanf (line, "%lld %lld", &node, &type);
                    if (retval != 2) {
                        P4EST_LERROR ("Premature end of file");
                        P4EST_FREE (line);
                        P4EST_FREE (vert);
                        P4EST_FREE (trees);
                        return NULL;
                    }

                    if (type == 5) {
                        int tagc, tag1, tag2;
                        retval = sscanf (line, "%lld %lld %d %d %d %d %d %d %d %d %d %d %d", &node, &type, &tagc, &tag1, &tag2,
                                         &trees[num_trees_real*P4EST_CHILDREN+0], &trees[num_trees_real*P4EST_CHILDREN+1], &trees[num_trees_real*P4EST_CHILDREN+3], &trees[num_trees_real*P4EST_CHILDREN+2],
                                         &trees[num_trees_real*P4EST_CHILDREN+4], &trees[num_trees_real*P4EST_CHILDREN+5], &trees[num_trees_real*P4EST_CHILDREN+7], &trees[num_trees_real*P4EST_CHILDREN+6] );
                        if (retval != 13) {
                            P4EST_LERROR ("Premature end of file");
                            P4EST_FREE (line);
                            P4EST_FREE (vert);
                            P4EST_FREE (trees);
                            return NULL;
                        }
                        for (int j = 0; j < P4EST_CHILDREN; j++) {
                            --trees[num_trees_real*P4EST_CHILDREN+j];
                        }
                        ++num_trees_real;
                    }
                }
            }
        }

        P4EST_FREE(line);
    }

    conn = p4est_connectivity_new(num_vertices, num_trees_real, 0, 0, 0, 0);
    memcpy(conn->vertices, vert, conn->num_vertices*3*sizeof(double));
    memcpy(conn->tree_to_vertex, trees, conn->num_trees*P4EST_CHILDREN*sizeof(p4est_topidx_t));

    P4EST_FREE(vert);
    P4EST_FREE(trees);

    /*
     * Fill tree_to_tree and tree_to_face to make sure we have a valid
     * connectivity.
     */
    for (tree = 0; tree < conn->num_trees; ++tree) {
        for (face = 0; face < P4EST_FACES; ++face) {
            conn->tree_to_tree[P4EST_FACES * tree + face] = tree;
            conn->tree_to_face[P4EST_FACES * tree + face] = face;
        }
    }
    P4EST_ASSERT (p4est_connectivity_is_valid (conn));

    /* Compute real tree_to_* fields and complete (edge and) corner fields. */
    p4est_connectivity_complete (conn);

    p4est_connectivity_set_attr(conn, sizeof(charm_tree_attr_t));
    charm_tree_attr_t * attr;
    for (tree = 0, attr = (charm_tree_attr_t*)(conn->tree_to_attr); tree < conn->num_trees; ++tree, ++attr) {
        for (face = 0; face < P4EST_FACES; ++face) {
            attr->bnd_type[face]    = CHARM_FACE_TYPE_INNER;
        }
        attr->region            = 0;
    }


    fsetpos (fid, &fpos);
    line = charm_connectivity_getline_upper(fid);
    int ele_count = 0;
    sscanf(line, "%d", &ele_count);
    P4EST_FREE(line);
    p4est_topidx_t fv[5];
    for (; ele_count > 0; ele_count--) {
        line = charm_connectivity_getline_upper(fid);
        if (charm_conn_parse_face(line, fv)) {
            charm_conn_find_tree_by_face(conn, fv, fv[4]);
        }
        P4EST_FREE(line);
    }

    retval = fclose (fid);
    fid = NULL;
    if (retval) {
        P4EST_LERRORF ("Failed to close %s\n", filename);
        goto dead;
    }

    P4EST_GLOBAL_PRODUCTIONF
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


charm_conn_reader_t charm_get_conn_reader(charm_ctx_t *ctx)
{
    return charm_conn_reader_msh;
}


p4est_connectivity_t* charm_conn_create(charm_ctx_t *ctx)
{
    charm_conn_reader_t     reader_fn;
    p4est_connectivity_t*   conn;
    const p4est_topidx_t    num_vertices    = 44;
    const p4est_topidx_t    num_trees       = 10;

    reader_fn = charm_get_conn_reader(ctx);

    conn = reader_fn("hole_3d_gmsh.msh");


    return conn;
}

