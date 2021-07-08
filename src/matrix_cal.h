#include "Header.h"
#include "ult.h"
#define DEBUG 0

//calculate det with MKL
double det_MKL(std::vector<std::vector<double>> &M){
    double *tmpM = (double*) mkl_malloc(M.size()*M.size()*sizeof(double),64);
    int *ip = (int*) mkl_malloc(M.size()*sizeof(double),64);
    double det = 1;
    for(int i = 0;i<M.size();++i){
        for(int j = 0;j<M.size();++j){
            tmpM[i*M.size()+j] = M[i][j];
        }
    }
    LAPACKE_dgetrf(LAPACK_ROW_MAJOR,M.size(),M.size(),tmpM,M.size(),ip);
    for(int i = 0;i<M.size();++i){
        det*=tmpM[i*M.size()+i];
        if(i+1!=ip[i]){
            det*=-1;
        }
    }
#if 0
    std::cout<<"ip:"<<std::endl;
    for(int i = 0;i<M.size();++i){
        std::cout<<ip[i]<<std::endl;
    }
    std::cout<<"---"<<std::endl;
#endif
    return det;
};

// home made det finder
double determinant(std::vector<std::vector<double>> &M,int n){
    double det = 0;
#if DEBUG
    std::cout<<"-------------------------"<<std::endl;
    for(int i = 0;i<n;++i){
        for (int j = 0;j<n;++j){
            std::cout<<M[i][j]<<',';
        }
        std::cout<<std::endl;
    }
    std::cout<<"-------------------------"<<std::endl;
#endif
    if(n==2){
        det = M[0][0]*M[1][1]-M[1][0]*M[0][1];
#if DEBUG
        std::cout<<"det : "<<det<<std::endl;
#endif
        return det;
    }
    std::vector<std::vector<double>> subM(n,std::vector<double>(n,0)); // temp matrix to store the sub matrix
    for(int a = 0;a<n;++a){ // go through all the elemet in first row of M
        int factor = 1;
        // copy the correct elements into sub matrix
        int subi = 0;
        int subj = 0; // the index for sub matrix
        int i = 0;
        int j = 0; // the index for main matrix

        for(subi = 0;subi<n-1;++subi){
            if( i == 0){
                i++;
            }
            for(subj = 0;subj<n-1;++subj){
                if( j == a){
                    j++;
                }
                subM[subi][subj] = M[i][j];
                j++;
            }
            i++;
            j = 0;
        }
        det += pow(-1,a)*M[0][a]*determinant(subM,n-1);
#if DEBUG
        std::cout<<"det : "<<det<<std::endl;
#endif
    }
    return det;
}
// Calculate the length of the edge. 
double get_length(node nd1, node nd2){
    return sqrt((nd1.x-nd2.x)*(nd1.x-nd2.x)+(nd1.y-nd2.y)*(nd1.y-nd2.y)+(nd1.z-nd2.z)*(nd1.z-nd2.z));
}
// Calculate b
double get_b(node nd1, node nd2, node nd3, int factor){
    std::vector<std::vector<double>> tmpM(3,std::vector<double>(3,0));
    tmpM[0] = {1,1,1};
    tmpM[1][0] = nd1.y;
    tmpM[1][1] = nd2.y;
    tmpM[1][2] = nd3.y;
    tmpM[2][0] = nd1.z;
    tmpM[2][1] = nd2.z;
    tmpM[2][2] = nd3.z;
    return det_MKL(tmpM)*factor;
    return determinant(tmpM,3)*factor;
};
// Calculate c
double get_c(node nd1, node nd2, node nd3, int factor){
    std::vector<std::vector<double>> tmpM(3,std::vector<double>(3,0));
    tmpM[0] = {1,1,1};
    tmpM[1][0] = nd1.x;
    tmpM[1][1] = nd2.x;
    tmpM[1][2] = nd3.x;
    tmpM[2][0] = nd1.z;
    tmpM[2][1] = nd2.z;
    tmpM[2][2] = nd3.z;
    return det_MKL(tmpM)*factor;
    return determinant(tmpM,3)*factor;
};
// Calculate d
double get_d(node nd1, node nd2, node nd3, int factor){
    std::vector<std::vector<double>> tmpM(3,std::vector<double>(3,0));
    tmpM[0] = {1,1,1};
    tmpM[1][0] = nd1.x;
    tmpM[1][1] = nd2.x;
    tmpM[1][2] = nd3.x;
    tmpM[2][0] = nd1.y;
    tmpM[2][1] = nd2.y;
    tmpM[2][2] = nd3.y;
    return det_MKL(tmpM)*factor;
    return determinant(tmpM,3)*factor;
};

// Calculate f
void get_f(std::vector<std::vector<double>> &f, std::vector<double> &b,  std::vector<double> &c,std::vector<double> &d){
    for(int i = 0;i<4;++i){
        for(int j = 0;j<4; ++j){
            f[i][j] = b[i]*b[j]+c[i]*c[j]+d[i]*d[j];
        }
    }
}
// Calculate the elementary matrix for each element
// using linear shape function
// based on the formula from Jin JianMing Book page181
// we need the following values
// V--- the volume of element_edge_i
// b,c,d---for both nodes of each of all six edges
// f --- for all six edges
// l --- the length of each of all six edges.
void get_elemental_matrix(element_edge element_edge_i, // the element of interest.
                            element element_node_i, // same element in node base.
                            std::vector<edge> &edges, // node ids of all edges,
                            std::vector<node> &nodes, // coordinates of all nodes,
                            std::vector<std::vector<double>> &E, // 
                            std::vector<std::vector<double>> &N){
    
    node tmpNode = {0,0,0};
    std::vector<node> nd_cords = {4,tmpNode};
    // first, get the node coordinates
    nd_cords[0].x = nodes[element_node_i.nd1].x;
    nd_cords[0].y = nodes[element_node_i.nd1].y;
    nd_cords[0].z = nodes[element_node_i.nd1].z;
    nd_cords[1].x = nodes[element_node_i.nd2].x;
    nd_cords[1].y = nodes[element_node_i.nd2].y;
    nd_cords[1].z = nodes[element_node_i.nd2].z;
    nd_cords[2].x = nodes[element_node_i.nd3].x;
    nd_cords[2].y = nodes[element_node_i.nd3].y;
    nd_cords[2].z = nodes[element_node_i.nd3].z;
    nd_cords[3].x = nodes[element_node_i.nd4].x;
    nd_cords[3].y = nodes[element_node_i.nd4].y;
    nd_cords[3].z = nodes[element_node_i.nd4].z;

    // Then, calculate the volume of the element
    std::vector<std::vector<double>> tmpM(4,std::vector<double>(4,0));
    tmpM[0] = {1,1,1,1};
    for(int i = 0; i<4; ++i){
        tmpM[1][i] = nd_cords[i].x;
        tmpM[2][i] = nd_cords[i].y;
        tmpM[3][i] = nd_cords[i].z;
    }
#if 0
    for(int i = 0; i<4; ++i){
        std::cout<<tmpM[i][0]<<','<<tmpM[i][1]<<','<<tmpM[i][2]<<','<<tmpM[i][3]<<std::endl;
    }
#endif
    double V;
    //V = 1/6.0*determinant(tmpM,4);
    V = fabs(1/6.0*det_MKL(tmpM));
    
#if 0
    std::cout<<"Volume from home made det : "<<1/6.0*determinant(tmpM,4)<<std::endl;
    std::cout<<"Volume from det_MKL       : "<<1/6.0*det_MKL(tmpM)<<std::endl;
#endif
    // Then calculate the length of each edge
    std::vector<double> l(6,0);
    for(int i = 0;i<6;++i){
        l[i] = get_length(nodes[edges[element_edge_i.edges[i]].nd1],nodes[edges[element_edge_i.edges[i]].nd2])*element_edge_i.direction[i];
        //std::cout<<"length of edge : "<<l[i]<<std::endl;
    }

    // Then calculate the b,c,d's of each node
    std::vector<double> b(4,0);
    std::vector<double> c(4,0);
    std::vector<double> d(4,0);
    
    
    b[0] = get_b(nodes[element_node_i.nd2],nodes[element_node_i.nd3],nodes[element_node_i.nd4],-1);
    c[0] = get_c(nodes[element_node_i.nd2],nodes[element_node_i.nd3],nodes[element_node_i.nd4],1);
    d[0] = get_d(nodes[element_node_i.nd2],nodes[element_node_i.nd3],nodes[element_node_i.nd4],-1);
    b[1] = get_b(nodes[element_node_i.nd1],nodes[element_node_i.nd3],nodes[element_node_i.nd4],1);
    c[1] = get_c(nodes[element_node_i.nd1],nodes[element_node_i.nd3],nodes[element_node_i.nd4],-1);
    d[1] = get_d(nodes[element_node_i.nd1],nodes[element_node_i.nd3],nodes[element_node_i.nd4],1);
    b[2] = get_b(nodes[element_node_i.nd1],nodes[element_node_i.nd2],nodes[element_node_i.nd4],-1);
    c[2] = get_c(nodes[element_node_i.nd1],nodes[element_node_i.nd2],nodes[element_node_i.nd4],1);
    d[2] = get_d(nodes[element_node_i.nd1],nodes[element_node_i.nd2],nodes[element_node_i.nd4],-1);
    b[3] = get_b(nodes[element_node_i.nd1],nodes[element_node_i.nd2],nodes[element_node_i.nd3],1);
    c[3] = get_c(nodes[element_node_i.nd1],nodes[element_node_i.nd2],nodes[element_node_i.nd3],-1);
    d[3] = get_d(nodes[element_node_i.nd1],nodes[element_node_i.nd2],nodes[element_node_i.nd3],1);
    //print_bcd(nodes,element_node_i,b,c,d);

    // Now construct bij, cij, dij for each edge for easier coding later
    std::vector<std::vector<double>> bij(6,std::vector<double>(2,0));
    std::vector<std::vector<double>> cij(6,std::vector<double>(2,0));
    std::vector<std::vector<double>> dij(6,std::vector<double>(2,0));
    bij[0] = {b[0],b[1]};
    cij[0] = {c[0],c[1]};
    dij[0] = {d[0],d[1]};
    bij[1] = {b[0],b[2]};
    cij[1] = {c[0],c[2]};
    dij[1] = {d[0],d[2]};
    bij[2] = {b[0],b[3]};
    cij[2] = {c[0],c[3]};
    dij[2] = {d[0],d[3]};
    bij[3] = {b[1],b[2]};
    cij[3] = {c[1],c[2]};
    dij[3] = {d[1],d[2]};
    bij[4] = {b[3],b[1]};
    cij[4] = {c[3],c[1]};
    dij[4] = {d[3],d[1]};
    bij[5] = {b[2],b[3]};
    cij[5] = {c[2],c[3]};
    dij[5] = {d[2],d[3]}; 
    

    // Now construct elementary matrix E
    for(int i = 0; i<6; ++i){
        for(int j = 0; j<6; ++j){
            E[i][j] = 4*l[i]*l[j]*V/(pow(6*V,4))*((cij[i][0]*dij[i][1]-dij[i][0]*cij[i][1])*(cij[j][0]*dij[j][1]-dij[j][0]*cij[j][1])+
                                                (dij[i][0]*bij[i][1]-bij[i][0]*dij[i][1])*(dij[j][0]*bij[j][1]-bij[j][0]*dij[j][1])+
                                                (bij[i][0]*cij[i][1]-cij[i][0]*bij[i][1])*(bij[j][0]*cij[j][1]-cij[j][0]*bij[j][1]));

        }
    }

    // Now construct matrix N.
    // To calculate this, we need some helper variable f
    std::vector<std::vector<double>> f(4,std::vector<double>(4,0)); // 4x4 for each node
    get_f(f,b,c,d); 

    N[0][0] = l[0]*l[0]/(360*V)*(f[1][1]-f[0][1]+f[0][0]);
    N[0][1] = l[0]*l[1]/(720*V)*(2*f[1][2]-f[1][0]-f[0][2]+f[0][0]);
    N[1][0] = N[0][1];
    N[0][2] = l[0]*l[2]/(720*V)*(2*f[1][3]-f[1][0]-f[0][3]+f[0][0]);
    N[2][0] = N[0][2];
    N[0][3] = l[0]*l[3]/(720*V)*(f[1][2]-f[1][1]-2*f[0][2]+f[0][1]);
    N[3][0] = N[0][3];
    N[0][4] = l[0]*l[4]/(720*V)*(f[1][1]-f[1][3]-f[0][1]+2*f[0][3]);
    N[4][0] = N[0][4];
    N[0][5] = l[0]*l[5]/(720*V)*(f[1][3]-f[1][2]-f[0][3]+f[0][2]);
    N[5][0] = N[0][5];

    N[1][1] = l[1]*l[1]/(360*V)*(f[2][2]-f[0][2]+f[0][0]);
    N[1][2] = l[1]*l[2]/(720*V)*(2*f[2][3]-f[0][2]-f[0][3]+f[0][0]);
    N[2][1] = N[1][2];
    N[1][3] = l[1]*l[3]/(720*V)*(f[2][2]-f[1][2]-f[0][2]+2*f[0][1]);
    N[3][1] = N[1][3];
    N[1][4] = l[1]*l[4]/(720*V)*(f[1][2]-f[2][3]-f[0][1]+f[0][3]);
    N[4][1] = N[1][4];
    N[1][5] = l[1]*l[5]/(720*V)*(f[0][2]-f[2][2]-2*f[0][3]+f[2][3]);
    N[5][1] = N[1][5];

    N[2][2] = l[2]*l[2]/(360*V)*(f[3][3]-f[0][3]+f[0][0]);
    N[2][3] = l[2]*l[3]/(720*V)*(f[2][3]-f[1][3]-f[0][2]+f[0][1]);
    N[3][2] = N[2][3];
    N[2][4] = l[2]*l[4]/(720*V)*(f[1][3]-f[3][3]-2*f[0][1]+f[0][3]);
    N[4][2] = N[2][4];
    N[2][5] = l[2]*l[5]/(720*V)*(f[3][3]-f[2][3]-f[0][3]+2*f[0][2]);
    N[5][2] = N[2][5];

    N[3][3] = l[3]*l[3]/(360*V)*(f[2][2]-f[1][2]+f[1][1]);
    N[3][4] = l[3]*l[4]/(720*V)*(f[1][2]-2*f[2][3]-f[1][1]+f[1][3]);
    N[4][3] = N[3][4];
    N[3][5] = l[3]*l[5]/(720*V)*(f[2][3]-f[2][2]-2*f[1][3]+f[1][2]);
    N[5][3] = N[3][5];

    N[4][4] = l[4]*l[4]/(360*V)*(f[1][1]-f[1][3]+f[3][3]);
    N[4][5] = l[4]*l[5]/(720*V)*(f[1][3]-2*f[1][2]-f[3][3]+f[2][3]);
    N[5][4] = N[4][5];

    N[5][5] = l[5]*l[5]/(360*V)*(f[3][3]-f[2][3]+f[2][2]);
    //std::cout<<"Volume of the element : "<<V<<std::endl;

};

void update_main_matrix(element_edge &elements_edge_i,
                        std::vector<std::vector<double>> &A,
                        std::vector<std::vector<double>> &E,
                        std::vector<std::vector<double>> &B,
                        std::vector<std::vector<double>> &N){
    for(int i = 0; i<E.size();++i){
        for(int j = 0; j<E[0].size();++j){
            A[elements_edge_i.edges[i]][elements_edge_i.edges[j]] += E[i][j];
            B[elements_edge_i.edges[i]][elements_edge_i.edges[j]] += N[i][j];
        }
    }

};

void print_matrix(std::vector<std::vector<double>> &M){
    if(M.size()==0){
        return;
    }
    std::cout<<"++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"<<std::endl;
    for(int i = 0;i<M.size();++i){
        for(int j = 0;j<M[0].size();++j){
            std::cout<<std::setprecision(3)<<std::setw(9)<<std::right<<M[i][j]<<',';
        }
        std::cout<<std::endl;
    }
    std::cout<<"++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"<<std::endl;
};

void compress_matrix(std::vector<std::vector<double>> &M, std::vector<int> &side_edge){
    std::vector<std::vector<double>> tmpM;
    int newsize = M.size()-side_edge.size();
    std::cout<<"Old size : "<<M.size()<<std::endl;
    std::cout<<"New size : "<<newsize<<std::endl;
    tmpM.resize(newsize,std::vector<double>(newsize,0));
    int oldI=0, oldJ=0;
    int side_edge_idI = 0;
    int side_edge_idJ = 0;
    
    for(int i = 0; i<newsize; ++i){
        while(oldI == side_edge[side_edge_idI]){
            oldI++;
            side_edge_idI++;
        }
        
        for(int j = 0; j<newsize; ++j){
            while(oldJ == side_edge[side_edge_idJ]){
                oldJ++;
                side_edge_idJ++;
            }
            tmpM[i][j] = M[oldI][oldJ];
            oldJ++;
        }
        oldJ = 0;
        side_edge_idJ = 0;
        oldI++;
    }
    M.resize(newsize);
    for(int i = 0;i<newsize;++i){
        M[i].resize(newsize);
    }
    std::cout<<"Resized old matrix..."<<std::endl;
    for(int i = 0; i<newsize; ++i){
        for(int j = 0; j<newsize; ++j){
            M[i][j] = tmpM[i][j];
        }
    }
};
void apply_boundary(std::vector<std::vector<double>> &M,std::vector<int> &side_edge){
    for(int i = 0; i < side_edge.size(); ++i){
        for(int j = 0; j< M.size(); ++j){
            M[side_edge[i]][j] = 0;
            M[j][side_edge[i]] = 0;
        }
        M[side_edge[i]][side_edge[i]] = 1;
    }
};

void decompose_B(std::vector<std::vector<double>> &B,std::vector<std::vector<double>> &LB){
    double *L = new double[B.size()*B.size()];
    
    for(int i = 0 ;i<B.size();++i){
        for(int j = 0 ;j<B.size();++j){
            L[i*B.size()+j] = B[i][j];
        }
    }
    // Do the Cholesky factorization.
    std::cout<<"decomposing B : ";
    std::cout<<LAPACKE_dpotrf(LAPACK_ROW_MAJOR,'L',B.size(),L,B.size())<<std::endl;

    for(int i = 0 ;i<B.size();++i){
        for(int j = 0 ;j<i+1;++j){
            LB[i][j] =  L[i*B.size()+j];
        }
    }
    for(int i = 0 ;i<B.size();++i){
        for(int j = i+1 ;j<B.size();++j){
            LB[i][j] =  0;
        }
    }
    free(L);
};
void transpose(std::vector<std::vector<double>> &M){
    double *tmpM = new double[M.size()*M.size()];
    
    for(int i = 0 ;i<M.size();++i){
        for(int j = 0 ;j<M.size();++j){
            tmpM[i*M.size()+j] = M[i][j];
        }
    }
    mkl_dimatcopy('r','t',M.size(),M[0].size(),1,tmpM,M[0].size(),M[0].size());
    for(int i = 0 ;i<M.size();++i){
        for(int j = 0 ;j<M.size();++j){
            M[i][j] = tmpM[i*M.size()+j];
        }
    }
    free(tmpM);
};

void inverse(std::vector<std::vector<double>> &M,
            const char &UL,
            const char &NU){
    double *tmpM = new double[M.size()*M.size()];
    
    for(int i = 0 ;i<M.size();++i){
        for(int j = 0 ;j<M.size();++j){
            tmpM[i*M.size()+j] = M[i][j];
        }
    }
    int n = M.size();
    int info = 0;
    LAPACKE_dtrtri(LAPACK_ROW_MAJOR,UL,NU,M.size(),tmpM,M.size());
    //std::cout<<info<<std::endl;
    for(int i = 0 ;i<M.size();++i){
        for(int j = 0 ;j<M.size();++j){
            M[i][j] = tmpM[i*M.size()+j];
        }
    }
    free(tmpM);
}
void check_unitary(std::vector<std::vector<double>> &A,
                      std::vector<std::vector<double>> &B){
    double *M1, *M2, *M3;
    M1 = (double *)mkl_malloc( A.size()*A.size()*sizeof( double ), 64 );
    M2 = (double *)mkl_malloc( A.size()*A.size()*sizeof( double ), 64 );
    M3 = (double *)mkl_malloc( A.size()*A.size()*sizeof( double ), 64 );
   
    for(int i = 0; i<A.size(); ++i){
        for(int j = 0; j < A.size(); ++j){
            M1[i*A.size()+j] = A[i][j];
            M2[i*A.size()+j] = B[i][j];
        }
    }
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 
                A.size(), A.size(), A.size(), 1, M1, A.size(),
                M2, A.size(), 0, M3, A.size());
    for(int i = 0; i<A.size(); ++i){
        for(int j = 0; j < A.size(); ++j){
            std::cout<<std::setw(10)<<M3[i*A.size()+j]<<',';
        }
        std::cout<<std::endl;
    }             
    std::cout<<"Matrix update finished. "<<std::endl;
    mkl_free(M1);
    mkl_free(M2);
    mkl_free(M3);
}
void get_final_matrix(std::vector<std::vector<double>> &LB,
                      std::vector<std::vector<double>> &A,
                      std::vector<std::vector<double>> &LBT){
    double *M1, *M2, *M3, *M4;    
    int max_threads = mkl_get_max_threads();

    M1 = (double *)mkl_malloc( A.size()*A.size()*sizeof( double ), 64 );
    M2 = (double *)mkl_malloc( A.size()*A.size()*sizeof( double ), 64 );
    M3 = (double *)mkl_malloc( A.size()*A.size()*sizeof( double ), 64 );
    M4 = (double *)mkl_malloc( A.size()*A.size()*sizeof( double ), 64 );
    std::cout<<"Copying temp matrices..."<<std::endl;

    for(int i = 0; i<A.size(); ++i){
        for(int j = 0; j < A.size(); ++j){
            M1[i*LB[0].size()+j] = LB[i][j];
            M2[i*LB[0].size()+j] = A[i][j];
            M3[i*LB[0].size()+j] = 0;
            M4[i*LB[0].size()+j] = 0;
        }
    }
    std::cout<<"Trying to use "<<max_threads<<" threads."<<std::endl;
    mkl_set_num_threads(max_threads);
    std::cout<<"Multiplicating matrices LB-1 x A..."<<std::endl;
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,A.size(), A.size(), A.size(), 
                1, M1, A.size(),
                M2, A.size(), 0, M3, A.size());
    std::cout<<"Copying intermediate matrix..."<<std::endl;
    for(int i = 0; i<A.size(); ++i){
        for(int j = 0; j < A.size(); ++j){
            M1[i*LB[0].size()+j] = LBT[i][j];
        }
    }
    std::cout<<"Calculating matrix A..."<<std::endl;
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 
                A.size(), A.size(), A.size(), 1, M3, A.size(),
                M1, A.size(), 0, M4, A.size());
    std::cout<<"Updating old matrix LB-1 x A x LBT-1..."<<std::endl;
    for(int i = 0; i<A.size(); ++i){
        for(int j = 0; j < A.size(); ++j){
            A[i][j] = M4[i*LB[0].size()+j];
        }
    }             
    std::cout<<"Matrix update finished. "<<std::endl;
    mkl_free(M1);
    mkl_free(M2);
    mkl_free(M3);
    mkl_free(M4);
};

void A_to_T(std::vector<std::vector<double>> &A,
             std::vector<double> &d,
             std::vector<double> &e,
             std::vector<double> &tau){
    double *tmpA, *tmpd, *tmpe, *tmptau;
    tmpA = (double *)mkl_malloc( A.size()*A.size()*sizeof( double ), 64 );
    tmpd = (double *)mkl_malloc( A.size()*sizeof( double ), 64 );
    tmpe = (double *)mkl_malloc( A.size()*sizeof( double ), 64 );
    tmptau = (double *)mkl_malloc( A.size()*sizeof( double ), 64 );
    std::cout<<"Copying A to tempA..."<<std::endl;
    for(int i = 0; i<A.size(); ++i){
        for(int j = 0; j < A.size(); ++j){
            tmpA[i*A.size()+j] = A[i][j];
        }
    }
    std::cout<<"Finished copy;"<<std::endl;
    std::cout<<"Turning A to T... : "<<std::endl;
    std::cout<<"Finished A->T : "<<LAPACKE_dsytrd(LAPACK_ROW_MAJOR,'L',A.size(),tmpA,A.size(), tmpd,tmpe,tmptau)<<std::endl;
    for(int i = 0; i<A.size(); ++i){
        for(int j = 0; j < A.size(); ++j){
            A[i][j] = tmpA[i*A.size()+j];
        }
    }
    for(int i = 0; i<A.size(); ++i){
        d[i] = tmpd[i];
        e[i] = tmpe[i];
        tau[i] = tmptau[i];
    }
    mkl_free(tmpA);
    mkl_free(tmpd);
    mkl_free(tmpe);
    mkl_free(tmptau);
};

void solve_eigen_T(std::vector<std::vector<double>> &T,
                    std::vector<double> &d,
                    std::vector<double> &e,
                    int * m,
                    int * isuppz,
                    double *w,
                    double *z,
                    double vl,double vu){
    double  *tmpd, *tmpe;
    
    tmpd = (double *)mkl_malloc( T.size()*sizeof( double ), 64 );
    tmpe = (double *)mkl_malloc( T.size()*sizeof( double ), 64 );
    for(int i = 0; i<T.size(); ++i){
        tmpd[i] = d[i];
        tmpe[i] = e[i];
    }
    std::cout<<"Solving... : ";
    std::cout<<LAPACKE_dstegr(LAPACK_ROW_MAJOR,'N','V',T.size(),tmpd,tmpe,vl,vu,1,T.size(),1,m,w,z,T.size(),isuppz)<<std::endl;
    std::cout<<"Found "<<*m<<" eigen value(s)."<<std::endl;
    mkl_free(tmpd);
    mkl_free(tmpe);
};

void print_eigen_values(int *m,double *w){
    for(int i = 0;i<*m;++i){
        std::cout<<"Eigen value id "<<i<<" : "<<std::endl<<
                    std::setw(20)<<"lambda = "<<w[i]<<std::endl<<
                    std::setw(20)<<"k = "<<sqrt(w[i])<<std::endl<<
                    std::setw(20)<<std::setprecision(10)<<"f = "<<std::setw(10)<<sqrt(w[i])/2/3.1415926*3e8<<std::endl;
    }
}
