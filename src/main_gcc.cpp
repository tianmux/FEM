#include "Header.h"
#include "ult.h"
#include "matrix_cal.h"
const double pi = M_PI;
const double c_light = 299792458;
int main(int argc, char *argv[]){
    

/* Open the file. NC_NOWRITE tells netCDF we want read-only access
	* to the file.*/

    std::string modelfile = argv[1];//"box.gen";
    std::string fn = "matrix.txt";
    const char* FILE_NAME = modelfile.c_str();

    int status, ncid, ndims, nvars, ngatts, dimids,unlimdimid;
    char recname[NC_MAX_NAME+1];
    size_t lenDim;
    
    std::string name_dim_numNodes = "num_nodes"; // the name of the dimension used by the coordx,y,z
    std::string name_var_coordx = "coordx";// the name of the variable represents the coordinatex, inquary this to see if it's using the id of num_nodes to double check.
    std::string name_var_coordy = "coordy";// the name of the variable represents the coordinatex, inquary this to see if it's using the id of num_nodes to double check.
    std::string name_var_coordz = "coordz";// the name of the variable represents the coordinatex, inquary this to see if it's using the id of num_nodes to double check.
    std::string name_dim_numElem = "num_elem"; // the name of the dimension used 
    std::string name_var_connect = "connect1"; // the name of the variable that stores the connection info of elements
    std::string name_dim_num_side_ss1 = "num_side_ss1"; // the name of the dimension for side set 1
    std::string name_var_elem_ss1 = "elem_ss1";// the name of the variable for elements in ss1.
    std::string name_var_side_ss1 = "side_ss1"; // the name of the variable for sides in ss1.

    int num_node_dim_id; // the id of the dimension "num_nodes
    int num_elem_dim_id; // the id of the dimension "num_elem"
    int num_side_ss1_dim_id; // the id of the dimension "num_side_ss1"

    int coordx_var_id; // the id of the variable "coordx_var"
    int coordy_var_id; // the id of the variable "coordy_var"
    int coordz_var_id; // the id of the variable "coordz_var"
    int connect_var_id; // the id of the variable "connect1"
    int elem_ss1_var_id; // the id of the variable "elem_ss1"
    int side_ss1_var_id; // the id of the variable "side_ss1"

    size_t num_node_len; // the length of the dimension "num_node", namely the number of nodes in this file.
    size_t num_elem_len; // the length fo the dimension "num_elem", namely the nubmer of elements in this model. 
    size_t num_side_ss1_len; // the length of the dimension "num_side_ss1", namely the number of sides in the side set 1, can 

    std::vector<node> nodes; // store the nodes. size should be num_node_len.
    std::vector<double> tmpx;
    std::vector<double> tmpy;
    std::vector<double> tmpz;
    std::vector<element> elements; // to store the elements. size should be num_elem.
    std::vector<edge> edges; // store the array of all the edges.size unknown. 
    std::vector<element_edge> elements_edge;
    std::vector<int> tmpElem; 
    std::vector<int> elem_ss1; // to store the element id that has at least one side in side set 1
    std::vector<int> side_ss1; // to store the LOCAL side id in the element that in side set 1.
                               // need to figure out the rule of numbering. 
                               // seems like the numbering rule is the following:
                               // side id 1 in the model corresponding to node ids : 1,2,4
                               // side id 2 in the model corresponding to node ids : 2,3,4
                               // side id 3 in the model corresponding to node ids : 3,4,1
                               // side id 4 in the model corresponding to node ids : 1,2,3
                               
    std::vector<int> side_edge; // to store the global id of the edges that are in side set 1.

// prepare the model
    // first, read the model in, including the node, element, side set
    get_ncid(ncid,FILE_NAME);
    get_general_info(ncid, ndims,  nvars,  ngatts,  unlimdimid);
    // Node
    get_node_info(ncid, name_dim_numNodes, name_var_coordx,name_var_coordy, name_var_coordz,coordx_var_id, coordy_var_id,coordz_var_id,num_node_dim_id, dimids,num_node_len);
    nodes.resize(num_node_len);
    tmpx.resize(num_node_len);
    tmpy.resize(num_node_len);
    tmpz.resize(num_node_len);
    read_nodes(ncid, coordx_var_id,coordy_var_id, coordz_var_id, num_node_len, nodes,tmpx,tmpy,tmpz);
    // Element
    get_elem_info(ncid, name_dim_numElem,name_var_connect,num_node_dim_id,connect_var_id,num_elem_len);
    elements.resize(num_elem_len);
    elements_edge.resize(num_elem_len);
    tmpElem.resize(num_elem_len*4);
    read_elem(ncid, connect_var_id, num_elem_len, tmpElem, elements);

    // side set
    get_side_set_info(ncid,name_dim_num_side_ss1,name_var_elem_ss1, name_var_side_ss1,num_side_ss1_dim_id, elem_ss1_var_id,side_ss1_var_id,num_side_ss1_len);
    elem_ss1.resize(num_side_ss1_len);
    side_ss1.resize(num_side_ss1_len);
    read_side_set(ncid,elem_ss1_var_id,side_ss1_var_id,num_side_ss1_len, elem_ss1, side_ss1);
    
    // Now we need to convert the node based_elements into edge_based elements
    convert_to_edge(elements,elements_edge,edges);
    
    convert_boundary(elem_ss1,side_ss1,edges,elements,side_edge);
#if DEBUG
    print_node_elements(elements);
    print_edge(edges);
    print_edge_elements(elements_edge,edges,side_edge);
    print_edges_on_boundary(side_edge);
#endif    
    
    // Now we have the node and node based elements converted to edge and edge based elements, 
    // we can start building the matrix 
    // First, we reserve the full space of matrix assuming it's dense
    std::vector<std::vector<double>> A;
    std::vector<std::vector<double>> B;
    std::vector<std::vector<double>> LB;
    std::vector<std::vector<double>> LBT;
    
    A.resize(edges.size(),std::vector<double>(edges.size(),0)); // num_edge x num_edge matrix elements
    B.resize(edges.size(),std::vector<double>(edges.size(),0));
    
    // then build the elemental matrix for each element
    std::vector<std::vector<double>> E;
    std::vector<std::vector<double>> N;
    
    E.resize(6,std::vector<double>(6,-1e9)); // 6 x 6 elements.
    N.resize(6,std::vector<double>(6,0)); // 6 x 6 elements.
    // go through all edge based elements.
    std::cout<<"Generating matrices..."<<std::endl;
    for (int i = 0; i<elements_edge.size();++i){
        get_elemental_matrix(elements_edge[i],elements[i],edges,nodes,E,N);
        update_main_matrix(elements_edge[i],A,E,B,N);
    }
    std::cout<<"Generated A and B matrices."<<std::endl;
    // Apply boundary condition
    // for now, only support PEC

    //print_bundary(elements,nodes,elem_ss1,side_ss1);
    
    apply_boundary(A,side_edge);
    apply_boundary(B,side_edge);
    //print_matrix(A);
    //print_matrix(B);
    std::cout<<"Compressing matrices to get rid of boundary rows and columns..."<<std::endl;
    compress_matrix(A,side_edge);
    compress_matrix(B,side_edge);
#if 0
    fn = "A.txt";
    out_put_maxtrix(A,fn);
    fn = "B.txt";
    out_put_maxtrix(B,fn);
#endif
    std::cout<<"Finish compressing matrices."<<std::endl;
    LB.resize(A.size(),std::vector<double>(A.size(),0));
    LBT.resize(A.size(),std::vector<double>(A.size(),0));
    
    std::cout<<"New size of A : "<<A.size()<<"x"<<A[0].size()<<std::endl;
    //print_matrix(A);
    //print_matrix(B);
    decompose_B(B,LB);
    for(int i = 0;i<LB.size();++i){
        for(int j = 0;j<LB[0].size();++j){
            LBT[i][j] = LB[i][j];
        }
    }
    // inverse LB
    inverse(LB,'L','N');
    std::cout<<"Inversed LB: "<<std::endl;

    //check_unitary(LB,LBT);

    transpose(LBT);   
#if DEBUG    
    print_matrix(LB);
    print_matrix(LBT);
#endif
    
    // inverse the LBT
    inverse(LBT,'U','N');
    std::cout<<"Inversed LBT: "<<std::endl;    
#if DEBUG
    print_matrix(LB);
    print_matrix(LBT);
#endif
    get_final_matrix(LB,A,LBT);
#if 0
    fn = "matrixA.txt";
    out_put_maxtrix(A,fn);
#endif 
    // We are only interested in the eigenvalues within certain range
    std::vector<double> d(A.size(),0); // store the diagonal elemets of T after we turn A into a tridiagonal matrix
    std::vector<double> e(A.size(),0); // store the off diagonal elemts of T.
    std::vector<double> tau(A.size(),0); // for elementary reflectors. 
    
    // convert the matrix A into a tridiagonal matrix T in order to do the bisection eigenvalue search later.
    A_to_T(A,d,e,tau);
#if 0
    fn = "matrixT.txt";
    out_put_maxtrix(A,fn); 
#endif
    // Now we got all the matrices we need, try to solve the eigen problem.
    // solve for eigenvalues
    int *m, *isuppz;
    double *w,*z;
    double vl,vu;
    m = (int *)mkl_malloc( sizeof( int ), 64 );
    isuppz = (int *)mkl_malloc( A.size()*sizeof( int ), 64 );
    w = (double *)mkl_malloc( A.size()*sizeof( double ), 64 );
    z = (double *)mkl_malloc( A.size()*sizeof( double ), 64 );
    std::cout<<"Lower freq : "<<std::atof(argv[2])<<std::endl;
    std::cout<<"Upper freq : "<<std::atof(argv[3])<<std::endl;
    
    vl = 2*pi*std::atof(argv[2])*2*pi*std::atof(argv[2])/c_light/c_light;
    vu = 2*pi*std::atof(argv[3])*2*pi*std::atof(argv[3])/c_light/c_light;
    solve_eigen_T(A,d,e,m,isuppz,w,z,vl,vu);
    print_eigen_values(m,w);
    return 0;
}