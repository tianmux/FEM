#include "Header.h"
#ifndef MY_ULT_FOR_FEA_H
#define MY_ULT_FOR_FEA_H

#define ERRCODE 2
#define ERR(e) {printf("Error: %s\n", nc_strerror(e)); /*exit(ERRCODE);*/}

std::string exitchar;

struct node{
    double x;
    double y;
    double z;
    std::vector<int> neighbors; // store the ids of the adjustant nodes
};

struct element{ // node based element, only support tet4 now
    // four ids of four nodes.
    int nd1;
    int nd2;
    int nd3;
    int nd4;
};

struct edge{
    int nd1;
    int nd2;
};

struct element_edge{ // edge based element, linear
/*
    int edge1; // store the global edge id, corresponding to local node id 1,2
    int edge2; // corresponding to local node id 1,3
    int edge3; // corresponding to local node id 1,4
    int edge4; // corresponding to local node id 2,3
    int edge5; // corresponding to local node id 4,2
    int edge6; // corresponding to local node id 3,4
    */
    std::vector<int> edges;
    std::vector<int> direction; // store the direction of each edges
};
// open the netCDF file.
void get_ncid(int &ncid,const char* FILE_NAME){
    int retval;//error code
	if ((retval = nc_open(FILE_NAME, NC_NOWRITE, &ncid)))
	{
		ERR(retval);
		std::cout << std::endl;
		std::cout << "Press Enter to Exit";
		std::getline(std::cin, exitchar);
		exit(ERRCODE);
	}
	printf("The File ID is: %d\n", ncid);
};

// get the general info of the file
void get_general_info(int &ncid,int &ndims, int &nvars, int &ngatts, int &unlimdimid){
    int retval;//error code
    std::cout<<ncid<<std::endl;
	if ((retval = nc_inq(ncid, &ndims, &nvars, &ngatts, &unlimdimid)))
	{
		ERR(retval);
		std::cout << std::endl;
		std::cout << "Press Enter to Exit";
		std::getline(std::cin, exitchar);
		exit(ERRCODE);
	}
};

// get the node info
void get_node_info(int &ncid,
                    std::string name_dim_numNodes,
                    std::string name_var_coordx,
                    std::string name_var_coordy,
                    std::string name_var_coordz,
                    int &coordx_var_id,
                    int &coordy_var_id,
                    int &coordz_var_id,
                    int &num_node_dim_id,
                    int &dimids,
                    size_t &num_node_len){
    int retval;//error code
    if ((retval = nc_inq_varid(ncid, name_var_coordx.c_str(),&coordx_var_id)))
    {
        ERR(retval);
        std::cout << std::endl;
        std::cout << "Press Enter to Exit";
        std::getline(std::cin, exitchar);
        exit(ERRCODE);
    }
    std::cout<<"The id of coordx : "<<coordx_var_id<<std::endl;
    if ((retval = nc_inq_varid(ncid, name_var_coordy.c_str(),&coordy_var_id)))
    {
        ERR(retval);
        std::cout << std::endl;
        std::cout << "Press Enter to Exit";
        std::getline(std::cin, exitchar);
        exit(ERRCODE);
    }
    std::cout<<"The id of coordy : "<<coordy_var_id<<std::endl;
    if ((retval = nc_inq_varid(ncid, name_var_coordz.c_str(),&coordz_var_id)))
    {
        ERR(retval);
        std::cout << std::endl;
        std::cout << "Press Enter to Exit";
        std::getline(std::cin, exitchar);
        exit(ERRCODE);
    }
    std::cout<<"The id of coordz : "<<coordz_var_id<<std::endl;

    if ((retval = nc_inq_var(ncid, coordx_var_id, NULL,NULL,NULL,&dimids,NULL)))
    {
        ERR(retval);
        std::cout << std::endl;
        std::cout << "Press Enter to Exit";
        std::getline(std::cin, exitchar);
        exit(ERRCODE);
    }
    std::cout<<"The id of dimension used by coords : "<<dimids<<std::endl;

    if ((retval = nc_inq_dimid(ncid, name_dim_numNodes.c_str(),&num_node_dim_id)))
    {
        ERR(retval);
        std::cout << std::endl;
        std::cout << "Press Enter to Exit";
        std::getline(std::cin, exitchar);
        exit(ERRCODE);
    }
    std::cout<<"The id of dimension Num_node : "<<num_node_dim_id<<std::endl;

    if(num_node_dim_id!=dimids){
        std::cout<<"Something is wrong, the id numbers did not check..."<<std::endl;
        exit(ERRCODE);
    }
    std::cout<<"Id number check, move on..."<<std::endl;

    if ((retval = nc_inq_dimlen(ncid, num_node_dim_id,&num_node_len)))
    {
        ERR(retval);
        std::cout << std::endl;
        std::cout << "Press Enter to Exit";
        std::getline(std::cin, exitchar);
        exit(ERRCODE);
    }
    std::cout << "Total number of nodes : " << num_node_len << std::endl;

};

void read_nodes(int &ncid, 
                int &coordx_var_id,
                int &coordy_var_id,
                int &coordz_var_id,
                size_t &num_node_len,
                std::vector<node> &nodes,
                std::vector<double> &tmpx,
                std::vector<double> &tmpy,
                std::vector<double> &tmpz){
    int retval;//error code
    if ((retval = nc_get_var_double(ncid, coordx_var_id, &tmpx[0])))//read in the value of node coordinates
		ERR(retval);
    if ((retval = nc_get_var_double(ncid, coordy_var_id, &tmpy[0])))//read in the value of node coordinates
		ERR(retval);
    if ((retval = nc_get_var_double(ncid, coordz_var_id, &tmpz[0])))//read in the value of node coordinates
		ERR(retval);
    for(int i = 0;i<num_node_len;++i){
        nodes[i].x = tmpx[i]; 
        nodes[i].y = tmpy[i];
        nodes[i].z = tmpz[i];
    }
#if 0
    for(int i = 0;i<num_node_len;++i){
        std::cout<<"node "<<i<<" : "<<std::setw(10)<<nodes[i].x<<','<<std::setw(10)<<nodes[i].y<<','<<std::setw(10)<<nodes[i].z<<std::endl;
    }
#endif
	std::cout << "Nodes info readin complete." << std::endl;
};

void get_elem_info(int &ncid, 
                    std::string name_dim_numElem,
                    std::string name_var_connect,
                    int &num_elem_dim_id,
                    int &connect_var_id,
                    size_t &num_elem_len){
    int retval;//error code
    if ((retval = nc_inq_varid(ncid, name_var_connect.c_str(),&connect_var_id)))
    {
        ERR(retval);
        std::cout << std::endl;
        std::cout << "Press Enter to Exit";
        std::getline(std::cin, exitchar);
        exit(ERRCODE);
    }
    std::cout<<"The id of connect : "<<connect_var_id<<std::endl;
    if ((retval = nc_inq_dimid(ncid, name_dim_numElem.c_str(),&num_elem_dim_id)))
    {
        ERR(retval);
        std::cout << std::endl;
        std::cout << "Press Enter to Exit";
        std::getline(std::cin, exitchar);
        exit(ERRCODE);
    }
    std::cout<<"The id of dimension num_elem : "<<num_elem_dim_id<<std::endl;

    if ((retval = nc_inq_dimlen(ncid, num_elem_dim_id,&num_elem_len)))
    {
        ERR(retval);
        std::cout << std::endl;
        std::cout << "Press Enter to Exit";
        std::getline(std::cin, exitchar);
        exit(ERRCODE);
    }
    std::cout << "Total number of elements : " << num_elem_len << std::endl;
};

void read_elem(int &ncid,
                int &connect_var_id,
                size_t &num_elem_len,
                std::vector<int> &tmpElem, 
                std::vector<element> &elements){
    int retval;//error code
    if ((retval = nc_get_var_int(ncid, connect_var_id, &tmpElem[0])))//read in the value of element node ids
    	ERR(retval);
    for(int i = 0; i<num_elem_len;++i){
        elements[i].nd1 = tmpElem[i*4]-1;// renumbering the node id from 1,2,3... to 0,1,2,...
        elements[i].nd2 = tmpElem[i*4+1]-1;
        elements[i].nd3 = tmpElem[i*4+2]-1;
        elements[i].nd4 = tmpElem[i*4+3]-1;
    }
    std::cout << "Elements info readin complete." << std::endl;
};

void get_side_set_info(int &ncid, 
                        std::string name_dim_num_side_ss1,
                        std::string name_var_elem_ss1,
                        std::string name_var_side_ss1,
                        int &num_side_ss1_dim_id,
                        int &elem_ss1_var_id,
                        int &side_ss1_var_id,
                        size_t &num_side_ss1_len){
    int retval;//error code
    if ((retval = nc_inq_varid(ncid, name_var_elem_ss1.c_str(),&elem_ss1_var_id)))
    {
        ERR(retval);
        std::cout << std::endl;
        std::cout << "Press Enter to Exit";
        std::getline(std::cin, exitchar);
        exit(ERRCODE);
    }
    if ((retval = nc_inq_varid(ncid, name_var_side_ss1.c_str(),&side_ss1_var_id)))
    {
        ERR(retval);
        std::cout << std::endl;
        std::cout << "Press Enter to Exit";
        std::getline(std::cin, exitchar);
        exit(ERRCODE);
    }    
    if ((retval = nc_inq_dimid(ncid, name_dim_num_side_ss1.c_str(),&num_side_ss1_dim_id)))
    {
        ERR(retval);
        std::cout << std::endl;
        std::cout << "Press Enter to Exit";
        std::getline(std::cin, exitchar);
        exit(ERRCODE);
    }

    if ((retval = nc_inq_dimlen(ncid, num_side_ss1_dim_id,&num_side_ss1_len)))
    {
        ERR(retval);
        std::cout << std::endl;
        std::cout << "Press Enter to Exit";
        std::getline(std::cin, exitchar);
        exit(ERRCODE);
    }
    std::cout << "Total number of sites in side set 1 : " << num_side_ss1_len << std::endl;
};
void read_side_set(int &ncid,
                int &elem_ss1_var_id,
                int &side_ss1_var_id,
                size_t &num_side_ss1_len,
                std::vector<int> &elem_ss1, 
                std::vector<int> &side_ss1){
    int retval;//error code
    if ((retval = nc_get_var_int(ncid, elem_ss1_var_id, &elem_ss1[0])))//read in the ids of elements in ss1
    	ERR(retval);
    if ((retval = nc_get_var_int(ncid, side_ss1_var_id, &side_ss1[0])))//read in the value of element node ids
    	ERR(retval);

    for(int i = 0; i<num_side_ss1_len;++i){
        elem_ss1[i] -=1;
        //std::cout<< elem_ss1[i]<<',';
    }
    std::cout<<std::endl;
    for(int i = 0; i<num_side_ss1_len;++i){
        side_ss1[i] -=1;
        //std::cout<< side_ss1[i]<<',';
    }
    std::cout<<std::endl;
    
    std::cout << "Side set info readin complete." << std::endl;
};

// find the global id of the edge in the array for all edges.
int find_edge_id(int &nd1, int &nd2, std::vector<edge> &edges){
    for(int i = 0; i<edges.size(); ++i){
        if((nd1 == edges[i].nd1 && nd2 == edges[i].nd2) || (nd1 == edges[i].nd2 && nd2 == edges[i].nd1)){
            return i; // not new
        }
    }
    edge tmpedge = {nd1,nd2};
    edges.push_back(tmpedge);
    return edges.size()-1; // is new
};
int find_edge_direction(int &nd1, int &nd2, std::vector<edge> &edges){
    for(int i = 0; i<edges.size(); ++i){
        if((nd1 == edges[i].nd1 && nd2 == edges[i].nd2)){
            return 1; 
        }
        else if(nd1 == edges[i].nd2 && nd2 == edges[i].nd1){
            return -1;
        }
    }
};
void convert_to_edge(std::vector<element> &elements,
                        std::vector<element_edge> &elements_edge,
                        std::vector<edge> &edges){
    for(int i = 0; i<elements.size(); ++i){
        elements_edge[i].edges.push_back(find_edge_id(elements[i].nd1,elements[i].nd2,edges));
        elements_edge[i].edges.push_back(find_edge_id(elements[i].nd1,elements[i].nd3,edges));
        elements_edge[i].edges.push_back(find_edge_id(elements[i].nd1,elements[i].nd4,edges));
        elements_edge[i].edges.push_back(find_edge_id(elements[i].nd2,elements[i].nd3,edges));
        elements_edge[i].edges.push_back(find_edge_id(elements[i].nd4,elements[i].nd2,edges));
        elements_edge[i].edges.push_back(find_edge_id(elements[i].nd3,elements[i].nd4,edges));
    }
    for(int i = 0; i<elements.size(); ++i){
        elements_edge[i].direction.push_back(find_edge_direction(elements[i].nd1,elements[i].nd2,edges));
        elements_edge[i].direction.push_back(find_edge_direction(elements[i].nd1,elements[i].nd3,edges));
        elements_edge[i].direction.push_back(find_edge_direction(elements[i].nd1,elements[i].nd4,edges));
        elements_edge[i].direction.push_back(find_edge_direction(elements[i].nd2,elements[i].nd3,edges));
        elements_edge[i].direction.push_back(find_edge_direction(elements[i].nd4,elements[i].nd2,edges));
        elements_edge[i].direction.push_back(find_edge_direction(elements[i].nd3,elements[i].nd4,edges));
    }
};

void convert_boundary(std::vector<int> &elem_ss1,
                    std::vector<int> &side_ss1,
                    std::vector<edge> &edges,
                    std::vector<element> &elements_node,
                    std::vector<int> &side_edge){
    // go through all the elements that has side on the boundary,
    // construct the edge that is on the boundary,
    // find the global id of the edge that is on the boundary and store it in the side_edge vector.
    std::vector<edge> tmpEdges = {{0,0},{0,0},{0,0}};
    for(int i = 0; i<elem_ss1.size(); ++i){
        switch(side_ss1[i]){
            case 0:
                tmpEdges[0] = {elements_node[elem_ss1[i]].nd1,elements_node[elem_ss1[i]].nd2};
                tmpEdges[1] = {elements_node[elem_ss1[i]].nd1,elements_node[elem_ss1[i]].nd4};
                tmpEdges[2] = {elements_node[elem_ss1[i]].nd2,elements_node[elem_ss1[i]].nd4};
                break;
            case 1:
                tmpEdges[0] = {elements_node[elem_ss1[i]].nd2,elements_node[elem_ss1[i]].nd3};
                tmpEdges[1] = {elements_node[elem_ss1[i]].nd2,elements_node[elem_ss1[i]].nd4};
                tmpEdges[2] = {elements_node[elem_ss1[i]].nd3,elements_node[elem_ss1[i]].nd4};
                break;
            case 2:
                tmpEdges[0] = {elements_node[elem_ss1[i]].nd1,elements_node[elem_ss1[i]].nd3};
                tmpEdges[1] = {elements_node[elem_ss1[i]].nd1,elements_node[elem_ss1[i]].nd4};
                tmpEdges[2] = {elements_node[elem_ss1[i]].nd3,elements_node[elem_ss1[i]].nd4};
                break;
            case 3:
                tmpEdges[0] = {elements_node[elem_ss1[i]].nd1,elements_node[elem_ss1[i]].nd2};
                tmpEdges[1] = {elements_node[elem_ss1[i]].nd1,elements_node[elem_ss1[i]].nd3};
                tmpEdges[2] = {elements_node[elem_ss1[i]].nd2,elements_node[elem_ss1[i]].nd3};
                break;
        }
        for(int j = 0; j < 3; ++j){
            side_edge.push_back(find_edge_id(tmpEdges[j].nd1,tmpEdges[j].nd2,edges));
        }
    }
    std::sort(side_edge.begin(),side_edge.end());
    side_edge.erase(std::unique(side_edge.begin(),side_edge.end()),side_edge.end());
}

void print_node_elements(std::vector<element> &elements){
    for(int i = 0; i<elements.size(); ++i){
                std::cout<< "elem "<<i<<" : "<<std::setw(10)<<
                elements[i].nd1<<','<< std::setw(10)<<
                elements[i].nd2<<','<< std::setw(10)<<
                elements[i].nd3<<','<< std::setw(10)<<
                elements[i].nd4<<std::endl;

    }
};
void print_edge(std::vector<edge> &edges){
    for(int i = 0;i<edges.size();++i){
        std::cout<<"Edge id "<<std::setw(5)<<std::right<<i<<" :"<<std::setw(5)<<edges[i].nd1<<"->"<<edges[i].nd2<<std::endl;
    }
}
void print_edge_elements(std::vector<element_edge> &elements_edge, std::vector<edge> &edges, std::vector<int> &side_edge){
    for(int i = 0; i<elements_edge.size(); ++i){
        std::cout<<std::setw(10)<<"Edge elem "<<std::setw(4)<<i<<" : ";
        for(int j = 0 ;j<elements_edge[0].edges.size();++j){
            std::cout<< std::setw(10)<<"Edge "<<std::setw(5)<<elements_edge[i].edges[j]<<": "<<
                    std::setw(6)<<
                    edges[elements_edge[i].edges[j]].nd1<<"->"<<
                    edges[elements_edge[i].edges[j]].nd2;
            if(elements_edge[i].direction[j]==-1){
                std::cout<<" (Ops)";
            }
            else{
                std::cout<<" (Pos)";
            }
            int k=0;
            while(k<side_edge.size() && elements_edge[i].edges[j]!=side_edge[k]){
                k++;
            }
            if(k==side_edge.size()){
                std::cout<<" (In)"<<" ;";
            }
            else{
                std::cout<<" (Ex)"<<" ;";
            }
        }  
        
        std::cout<<std::endl;
    }
}

void print_bundary(std::vector<element> &elements_node,
                    std::vector<node> &nodes,
                    std::vector<int> &elem_ss1,
                    std::vector<int> &side_ss1){
    for(int i = 0; i<elem_ss1.size(); ++i){
        std::cout<<"The id of element on ss1 : "<<elem_ss1[i]<<std::endl;
        std::cout<<"The node list of this element : "<<elements_node[elem_ss1[i]].nd1<<" : "
        <<nodes[elements_node[elem_ss1[i]].nd1].x<<" , "<<
        nodes[elements_node[elem_ss1[i]].nd1].y<<" , "<<
        nodes[elements_node[elem_ss1[i]].nd1].z<<std::endl;
        std::cout<<"The node list of this element : "<<elements_node[elem_ss1[i]].nd2<<" : "
        <<nodes[elements_node[elem_ss1[i]].nd2].x<<" , "<<
        nodes[elements_node[elem_ss1[i]].nd2].y<<" , "<<
        nodes[elements_node[elem_ss1[i]].nd2].z<<std::endl;
        std::cout<<"The node list of this element : "<<elements_node[elem_ss1[i]].nd3<<" : "
        <<nodes[elements_node[elem_ss1[i]].nd3].x<<" , "<<
        nodes[elements_node[elem_ss1[i]].nd3].y<<" , "<<
        nodes[elements_node[elem_ss1[i]].nd3].z<<std::endl;
        std::cout<<"The node list of this element : "<<elements_node[elem_ss1[i]].nd4<<" : "
        <<nodes[elements_node[elem_ss1[i]].nd4].x<<" , "<<
        nodes[elements_node[elem_ss1[i]].nd4].y<<" , "<<
        nodes[elements_node[elem_ss1[i]].nd4].z<<std::endl;
        std::cout<<"The local ID of the side that is on ss1 : "<<side_ss1[i]<<std::endl;
    }
}

void print_edges_on_boundary(std::vector<int> &side_edges){
    for(int i = 0; i < side_edges.size(); ++i){
        std::cout<<side_edges[i]<<std::endl;
    }
};
void out_put_maxtrix(std::vector<std::vector<double>> &M,std::string &fn){
    std::ofstream out_fn(fn, std::ofstream::out);

    for(int i = 0;i<M.size();++i){
        for(int j = 0;j<M[0].size();++j){
            out_fn<<std::setprecision(3)<<std::setw(9)<<std::right<<M[i][j]<<',';
        }
        out_fn<<std::endl;
    }
};

void print_bcd(std::vector<node> &nodes, 
                element &element_node_i, 
                std::vector<double> &b, 
                std::vector<double> &c, 
                std::vector<double> &d){
    std::cout<<"Node1 : "<<std::setw(8)<<nodes[element_node_i.nd1].x<<','
                         <<std::setw(8)<<nodes[element_node_i.nd1].y<<','
                         <<std::setw(8)<<nodes[element_node_i.nd1].z<<std::endl;
    std::cout<<"Node2 : "<<std::setw(8)<<nodes[element_node_i.nd2].x<<','
                         <<std::setw(8)<<nodes[element_node_i.nd2].z<<','
                         <<std::setw(8)<<nodes[element_node_i.nd2].y<<std::endl;
    std::cout<<"Node3 : "<<std::setw(8)<<nodes[element_node_i.nd3].x<<','
                         <<std::setw(8)<<nodes[element_node_i.nd3].y<<','
                         <<std::setw(8)<<nodes[element_node_i.nd3].z<<std::endl;
    std::cout<<"Node4 : "<<std::setw(8)<<nodes[element_node_i.nd4].x<<','
                         <<std::setw(8)<<nodes[element_node_i.nd4].y<<','
                         <<std::setw(8)<<nodes[element_node_i.nd4].z<<std::endl;
    for(int  i = 0;i<b.size(); ++i){
        std::cout<<"b_"<<i<<" : "<<b[i]<<std::endl;
        std::cout<<"c_"<<i<<" : "<<c[i]<<std::endl;
        std::cout<<"d_"<<i<<" : "<<d[i]<<std::endl;
    }
}
#endif 
