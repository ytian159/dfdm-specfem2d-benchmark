#pragma once

#include <mpi.h>
#include "matrix.hpp"

/*
the idea behind the tag is that each process will have unique list of elements
hence a combination of an element and its owner process ID would be unique. Sp
as long as the number of elements per process don't exceed 255, the tag should
be unique. The direction and type mask are used to identify which matrix and
boundary direction is being shared.
*/

const int TAG_DIR_MASK =  0x00300000;   // 2 bits for boundary direction
const int TAG_TYPE_MASK = 0x000F0000;  // 4 bits for boundary type
const int TAG_PROC_MASK = 0x0000FF00;  // 8 bits for process ID, each process can have at most 4 neighbor processes, so this should work.
const int TAG_ELE_MASK =  0x000000FF;  // 8 bits for element ID, this caps the number of elements per process to 255 for now
const int maxtagsizeofi = 33554431; // obtained from MPI docs, this number is for slingshot NICS (SS-11) when suing OFI


namespace DFDM{
    enum BdryType {
        U12,
        U21,
        SXX11,
        SXX22,
        SZZ11,
        SZZ22,
    };
    enum BdryDirection {
        MO,
        PO,
        OM,
        OP
    };
    
    void mpi_init();
    void mpi_barrier();
    void mpi_finalize();
    int get_rank_id();
    int get_total_ranks();
    void recv_mat(int src_rank, int32_t tag, DFDM::matrix<double>& mat);
    void recv_mat_nb(int src_rank, int my_el_id, DFDM::matrix<double>& mat, MPI_Request& request_1, MPI_Request& request_2);
    void recv_matnb_dim(int src_rank, int32_t tag, DFDM::matrix<double>& mat, MPI_Request& request);
    void send_mat(int target_rank, int32_t tag, DFDM::matrix<double>& mat);
    void mpi_wait_all(uint32_t ccount, std::vector<MPI_Request>& request, std::vector<MPI_Status>& status);
    void mpi_wait(MPI_Request& reqeust1, MPI_Request& request2);
    void send_mat_nb(int target_rank, int target_el, DFDM::matrix<double>& mat, MPI_Request& request_1, MPI_Request& request_2);
    void send_matnb_dim(int target_rank, int32_t tag, DFDM::matrix<double>& mat, MPI_Request& request);
    int32_t generateTag(BdryType boundaryType, BdryDirection boundaryDirection, uint32_t procId, uint32_t eleid);
}