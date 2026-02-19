#include "mpi_helper.hpp"
#include <cmath>
void DFDM::mpi_init(){
    MPI_Init(NULL, NULL);
}

void DFDM::mpi_barrier() {
  MPI_Barrier(MPI_COMM_WORLD);
}
void DFDM::mpi_finalize(){
    MPI_Finalize();
}

int DFDM::get_rank_id(){
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    return rank;
}

int DFDM::get_total_ranks(){
    int comm_size;
    MPI_Comm_size(MPI_COMM_WORLD, &comm_size);

    return comm_size;
}

void DFDM::recv_mat(int src_rank, int32_t tag, DFDM::matrix<double>& mat){ // if an element is receiving, it will provide its id as tag
  std::vector<uint32_t> dims(2);
  uint32_t rec_tag_dim = tag;
  uint32_t rec_tag_data = tag ^ (1<< 22);

  if (rec_tag_data > maxtagsizeofi){
    std::cout << "MPI::ERROR::TAG SIZE larger than allowable on OFI with SS-11" << std::endl;
    exit(-1);
  }

  MPI_Recv((void*)dims.data(), dims.size(), MPI_INT, src_rank, rec_tag_dim, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  mat.rows = dims[0];
  mat.cols = dims[1];
  mat.resize(mat.rows, mat.cols);
  MPI_Recv((void*)mat.data_vec.data(), mat.data_vec.size(), MPI_DOUBLE, src_rank, rec_tag_data, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
}
//non blocking receive
void DFDM::recv_mat_nb(int src_rank, int32_t tag, DFDM::matrix<double>& mat, MPI_Request& request_1, MPI_Request& request_2){ // if an element is receiving, it will provide its id as tag
  std::vector<uint32_t> dims(2);  
  
  uint32_t rec_tag_dim = tag;
  uint32_t rec_tag_data = tag ^ (1<< 22);

  if (rec_tag_data > maxtagsizeofi || rec_tag_dim > maxtagsizeofi){
    std::cout << "MPI::ERROR::TAG SIZE larger than allowable on OFI with SS-11" << std::endl;
    exit(-1);
  }

  MPI_Status status;


  MPI_Irecv((void*)dims.data(), dims.size(), MPI_UINT32_T, src_rank, rec_tag_dim, MPI_COMM_WORLD, &request_1);
  MPI_Wait(&request_1,&status);
  mat.rows = dims[0];
  mat.cols = dims[1];

  if (dims[0] == 0 || dims[1] == 0) {
    int received_count;
    MPI_Get_count(&status, MPI_INT, &received_count);
    printf("Received %d integers from rank %d with tag %d\n",received_count, status.MPI_SOURCE, status.MPI_TAG);
    std::cerr << "MPI::ERROR::Received matrix dimensions are zero. "
              << "Source Rank: " << src_rank 
              << ", Target Rank: " << get_rank_id() 
              << ", Tag: " << tag 
              << ", Dimensions: [" << dims[0] << ", " << dims[1] << "]" << std::endl;
    exit(-1);
  }
  mat.resize(mat.rows, mat.cols);
  MPI_Irecv((void*)mat.data_vec.data(), mat.data_vec.size(), MPI_DOUBLE, src_rank, rec_tag_data, MPI_COMM_WORLD, &request_2);
}

// non blocking receive with the assumption that the receiver knows the dimension of incoming matrix.
void DFDM::recv_matnb_dim(int src_rank, int32_t tag, DFDM::matrix<double>& mat, MPI_Request& request){ 
  uint32_t rec_tag_data = tag;

  if (rec_tag_data > maxtagsizeofi) {
    std::cerr << "MPI::ERROR::Tag size exceeds the maximum allowable size. "
              << "Source Rank: " << src_rank 
              << ", Target Rank: " << get_rank_id() 
              << ", Tag Data: " << rec_tag_data 
              << ", Maximum Allowed Tag Size: " << maxtagsizeofi << std::endl;
    exit(-1);
  }

  

  if (mat.rows == 0 || mat.cols == 0) {
    std::cerr << "MPI::ERROR::Matrix dimensions are zero. "
              << "Source Rank: " << src_rank 
              << ", Target Rank: " << get_rank_id() 
              << ", Tag: " << tag 
              << ", Dimensions: [" << mat.rows << ", " << mat.cols << "]" << std::endl;
    exit(-1);
  }

  size_t message_size = mat.data_vec.size() * sizeof(double);
  if (message_size > 16364) {
    std::cerr << "MPI::WARNING::Message size exceeds 16364 bytes. "
              << "Rendezvous protocol may be used, which could lead to a potential deadlock. "
              << "Message Size: " << message_size << " bytes." << std::endl;
  }

  MPI_Irecv((void*)mat.data_vec.data(), mat.data_vec.size(), MPI_DOUBLE, src_rank, rec_tag_data, MPI_COMM_WORLD, &request);
}


void DFDM::send_mat(int target_rank, int32_t tag, DFDM::matrix<double>& mat){ // if an element is sending, it will provide the destination element id as tag
  std::vector<uint32_t> dims(2);
  dims[0] = mat.rows;
  dims[1] = mat.cols;
  int32_t rec_tag_dim = tag;
  int32_t rec_tag_data = tag ^ (1<< 22);

  if (rec_tag_data > maxtagsizeofi || rec_tag_dim > maxtagsizeofi){
    std::cout << "MPI::ERROR::TAG SIZE larger than allowable on OFI with SS-11" << std::endl;
    exit(-1);
  }

  MPI_Send((void*)dims.data(), dims.size(), MPI_INT, target_rank, rec_tag_dim, MPI_COMM_WORLD); // using my element_id as a tag, receiver will use target id
  MPI_Send((void*)mat.data_vec.data(), mat.data_vec.size(), MPI_DOUBLE, target_rank, rec_tag_data, MPI_COMM_WORLD);
}

//non blocking send
void DFDM::send_mat_nb(int target_rank, int32_t tag, DFDM::matrix<double>& mat, MPI_Request& request_1, MPI_Request& request_2){ // if an element is sending, it will provide the destination element id as tag
  std::vector<uint32_t> dims(2);
  dims[0] = mat.rows;
  dims[1] = mat.cols;
  int32_t rec_tag_dim = tag;
  int32_t rec_tag_data = tag ^ (1<< 22);

  if (rec_tag_data > maxtagsizeofi || rec_tag_dim > maxtagsizeofi) {
    std::cerr << "MPI::ERROR::Tag size exceeds the maximum allowable size. "
              << "Target Rank: " << target_rank 
              << ", Source Rank: " << get_rank_id() 
              << ", Tag Dim: " << rec_tag_dim 
              << ", Tag Data: " << rec_tag_data 
              << ", Maximum Allowed Tag Size: " << maxtagsizeofi << std::endl;
    exit(-1);
  }

  if (dims[0] == 0 || dims[1] == 0) {
    std::cerr << "MPI::ERROR::Matrix dimensions being sent are zero. "
              << "Target Rank: " << target_rank 
              << ", Source Rank: " << get_rank_id() 
              << ", Tag: " << tag 
              << ", Dimensions: [" << dims[0] << ", " << dims[1] << "]" << std::endl;
    exit(-1);
  }

  // std::cout << "Sending matrix dimensions to rank: " << target_rank
  //           << ", Current Rank: " << get_rank_id()
  //           << ", Tag Dim: " << rec_tag_dim
  //           << ", Tag Data: " << rec_tag_data
  //           << ", Dimensions: [" << dims[0] << ", " << dims[1] << "]"
  //           << ", Maximum Allowed Tag Size: " << maxtagsizeofi << std::endl;
            
  size_t message_size = dims.size() * sizeof(uint32_t) + mat.data_vec.size() * sizeof(double);
  if (message_size > 16364) {
    std::cerr << "MPI::WARNING::Message size exceeds 16364 bytes. "
              << "Rendezvous protocol may be used, which could lead to a potential deadlock. "
              << "Message Size: " << message_size << " bytes." << std::endl;
  }
  MPI_Isend((void*)dims.data(), dims.size(), MPI_UINT32_T, target_rank, rec_tag_dim, MPI_COMM_WORLD, &request_1); // when sending, using target element id as tag, receiver will use my element id as tag
  MPI_Wait(&request_1, MPI_STATUS_IGNORE);
  MPI_Isend((void*)mat.data_vec.data(), mat.data_vec.size(), MPI_DOUBLE, target_rank, rec_tag_data, MPI_COMM_WORLD, &request_2);
  MPI_Wait(&request_2, MPI_STATUS_IGNORE);
}

//non blocking send with the assumption that the receiver knows the dimension of incoming matrix.
void DFDM::send_matnb_dim(int target_rank, int32_t tag, DFDM::matrix<double>& mat, MPI_Request& request){ // if an element is sending, it will provide the destination element id as tag
  int32_t rec_tag_data = tag;

  if (rec_tag_data > maxtagsizeofi) {
    std::cerr << "MPI::ERROR::Tag size exceeds the maximum allowable size. "
              << "Target Rank: " << target_rank 
              << ", Source Rank: " << get_rank_id() 
              << ", Tag Data: " << rec_tag_data 
              << ", Maximum Allowed Tag Size: " << maxtagsizeofi << std::endl;
    exit(-1);
  }

  if (mat.rows == 0 || mat.cols == 0) {
    std::cerr << "MPI::ERROR::Matrix dimensions being sent are zero. "
              << "Target Rank: " << target_rank 
              << ", Source Rank: " << get_rank_id() 
              << ", Tag: " << tag 
              << ", Dimensions: [" << mat.rows << ", " << mat.cols << "]" << std::endl;
    exit(-1);
  }

  size_t message_size = mat.data_vec.size() * sizeof(double);
  if (message_size > 16364) {
    std::cerr << "MPI::WARNING::Message size exceeds 16364 bytes. "
              << "Rendezvous protocol may be used, which could lead to a potential deadlock. "
              << "Message Size: " << message_size << " bytes." << std::endl;
  }

  MPI_Isend((void*)mat.data_vec.data(), mat.data_vec.size(), MPI_DOUBLE, target_rank, rec_tag_data, MPI_COMM_WORLD, &request);
  // MPI_Wait(&request, MPI_STATUS_IGNORE);
}


void DFDM::mpi_wait_all(uint32_t ccount, std::vector<MPI_Request>& request, std::vector<MPI_Status>& status){
  switch(ccount){
    case 1:
      MPI_Wait(&request[1], MPI_STATUS_IGNORE);
      break;
    case 2:
      MPI_Wait(&request[3], MPI_STATUS_IGNORE);
      break;
    case 3:
      MPI_Wait(&request[1], MPI_STATUS_IGNORE);
      MPI_Wait(&request[3], MPI_STATUS_IGNORE);
      break;
    default:
      break;
  }  
}

// Function to wrap MPI_Waitall with MPI_STATUS_IGNORE
void DFDM::mpi_wait(MPI_Request& request1, MPI_Request& request2){
    int err1 = MPI_Wait(&request1, MPI_STATUS_IGNORE);   
    if (err1!= MPI_SUCCESS) {
        std::cerr << "MPI_Waitall failed with error code " << err1 << std::endl;
        exit(-1);
    }

    int err2 = MPI_Wait(&request2, MPI_STATUS_IGNORE);
    
    if (err2!= MPI_SUCCESS) {
        std::cerr << "MPI_Waitall failed with error code " << err2 << std::endl;
        exit(-1);
    }
}

  /**
   * Constructs a unique tag by combining multiple components into a single 32-bit integer.
   * The tag is used for MPI communication to uniquely identify messages based on boundary
   * type, boundary direction, process ID, and element ID. Each component is masked and shifted
   * to ensure no overlap between their bits.
   *
   * The masks and their purposes are as follows:
   * - TAG_DIR_MASK (0x00300000): Allocates 2 bits for the boundary direction.
   * - TAG_TYPE_MASK (0x000F0000): Allocates 4 bits for the boundary type.
   * - TAG_PROC_MASK (0x0000FF00): Allocates 8 bits for the process ID. This allows up to 256 unique process IDs.
   * - TAG_ELE_MASK (0x000000FF): Allocates 8 bits for the element ID. This caps the number of elements per process to 255.
   *
   * The bitwise operations ensure that:
   * - `boundaryType` is shifted left by 20 bits and masked with TAG_DIR_MASK.
   * - `boundaryDirection` is shifted left by 16 bits and masked with TAG_TYPE_MASK.
   * - `procId` is shifted left by 8 bits and masked with TAG_PROC_MASK.
   * - `eleid` is masked with TAG_ELE_MASK.
   *
   * The resulting tag is a combination of these components, ensuring uniqueness for MPI communication.
   */
int32_t DFDM::generateTag(BdryType boundaryType, BdryDirection boundaryDirection, uint32_t procId, uint32_t eleid){
  if (boundaryType > 15 || boundaryDirection > 4 || procId > 255 || eleid > 255){
    std::cerr << "generateTag: Invalid boundary type, direction, procId, or eleid. "
              << std::endl;
    return -1;
  }


  int32_t tag = ((boundaryDirection << 20) & TAG_DIR_MASK) | 
         ((boundaryType << 16) & TAG_TYPE_MASK) | 
         ((procId << 8) & TAG_PROC_MASK) | 
         (eleid & TAG_ELE_MASK);

  if (tag > maxtagsizeofi){
    std::cout << "MPI::ERROR::TAG SIZE larger than allowable on OFI with SS-11" << std::endl;
    exit(-1);
  }

  return tag;
}
