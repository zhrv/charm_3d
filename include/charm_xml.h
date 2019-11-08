//
// Created by appmath on 29.10.17.
//

#ifndef CHAMR_3D_CHARM_XML_H
#define CHAMR_3D_CHARM_XML_H

#include "mxml.h"


mxml_node_t* charm_xml_node_get_child_or_null(mxml_node_t * n, char* name);
mxml_node_t* charm_xml_node_get_child(mxml_node_t *, char*);
mxml_node_t* charm_xml_node_get_next_child(mxml_node_t * current, mxml_node_t * n, char* name);

int charm_xml_node_value_dbl(mxml_node_t *, double*);
int charm_xml_node_value_int(mxml_node_t *, int*);
int charm_xml_node_value_str(mxml_node_t *, char*);

int charm_xml_node_attr_dbl(mxml_node_t *, char*, double*);
int charm_xml_node_attr_int(mxml_node_t *, char*, int*);
int charm_xml_node_attr_str(mxml_node_t *, char*, char*);

int charm_xml_node_child_value_dbl(mxml_node_t *, char*, double*);
int charm_xml_node_child_value_int(mxml_node_t *, char*, int*);
int charm_xml_node_child_value_str(mxml_node_t *, char*, char*);

int charm_xml_node_child_param_dbl(mxml_node_t *, char*, double*);
int charm_xml_node_child_param_int(mxml_node_t *n, char* name, int* val);
int charm_xml_node_child_param_str(mxml_node_t *n, char* name, char* val);


#endif //CHAMR_3D_CHARM_XML_H
