//
// Created by appmath on 29.10.17.
//

#include "charm_xml.h"
#include "charm_globals.h"


mxml_node_t* charm_xml_node_get_child_or_null(mxml_node_t * n, char* name)
{
    mxml_node_t *node = mxmlFindElement(n, n, name, NULL, NULL, MXML_DESCEND);
    return node;
}

mxml_node_t* charm_xml_node_get_child(mxml_node_t * n, char* name)
{
    mxml_node_t *node = mxmlFindElement(n, n, name, NULL, NULL, MXML_DESCEND);
    if (node == NULL) {
        CHARM_GLOBAL_LERRORF("Undefined parameter '%s'\n", name);
        charm_abort(NULL, 1);
    }
    return node;
}

mxml_node_t* charm_xml_node_get_next_child(mxml_node_t * current, mxml_node_t * n, char* name)
{
    return mxmlFindElement(current, n, name,
                           NULL, NULL,
                           MXML_DESCEND);
}

int charm_xml_node_value_dbl(mxml_node_t *n, double*d)
{
    int whitespace = 0;
    const char* str = mxmlGetText(n, &whitespace);
    *d = atof(str);
    return 0;
}


int charm_xml_node_value_int(mxml_node_t *n, int*val)
{
    int whitespace = 0;
    const char* str = mxmlGetText(n, &whitespace);
    *val = atoi(str);
    return 0;
}


int charm_xml_node_value_str(mxml_node_t *n, char*val)
{
    int whitespace = 0;
    const char* str = mxmlGetText(n, &whitespace);
    strcpy(val, str);
    return 0;
}



int charm_xml_node_attr_dbl(mxml_node_t *n, char*name, double*val)
{
    const char* sval = mxmlElementGetAttr(n, name);
    *val = atof(sval);
    return 1;
}


int charm_xml_node_attr_int(mxml_node_t *n, char*name, int*val)
{
    const char* sval = mxmlElementGetAttr(n, name);
    *val = atoi(sval);
    return 1;
}


int charm_xml_node_attr_str(mxml_node_t *n, char*name, char*val)
{
    const char* sval = mxmlElementGetAttr(n, name);
    strcpy(val, sval);
    return 1;
}



int charm_xml_node_child_dbl(mxml_node_t *n, char *name, double *val)
{
    mxml_node_t* node = charm_xml_node_get_child(n, name);
    return charm_xml_node_value_dbl(node, val);
}


int charm_xml_node_child_value_int(mxml_node_t *n, char*name, int*val)
{

}


int charm_xml_node_child_value_str(mxml_node_t *n, char*name, char*val)
{

}


int charm_xml_node_child_param_dbl(mxml_node_t *n, char* name, double* val)
{
    mxml_node_t *node = charm_xml_node_get_child(n, name);
    return charm_xml_node_attr_dbl(node, "value", val);
}


int charm_xml_node_child_param_int(mxml_node_t *n, char* name, int* val)
{
    mxml_node_t *node = charm_xml_node_get_child(n, name);
    return charm_xml_node_attr_int(node, "value", val);
}


int charm_xml_node_child_param_str(mxml_node_t *n, char* name, char* val)
{
    mxml_node_t *node = charm_xml_node_get_child(n, name);
    return charm_xml_node_attr_str(node, "value", val);
}


