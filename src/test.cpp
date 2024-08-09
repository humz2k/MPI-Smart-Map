#include "smart-map.hpp"
#include <mpi.h>

class AllToAll : public SmartMap {
  private:
    int m_size;

  public:
    AllToAll(int size, MPI_Comm comm) : m_size(size) {
        int comm_size;
        MPI_Comm_size(comm, &comm_size);
        compile(m_size * comm_size, comm);
    }

    smart_map_index map(smart_map_index in) {
        int rank = in.idx / m_size;
        int idx = in.rank * (m_size) + (in.idx % m_size);
        return smart_map_index(rank, idx);
    }
};

int main() {
    MPI_Init(NULL, NULL);

    int comm_size;
    int comm_rank;
    MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &comm_rank);
    int alltoall_size = 10;

    AllToAll map(alltoall_size, MPI_COMM_WORLD);

    std::vector<int> test_in;
    for (int i = 0; i < alltoall_size * comm_size; i++) {
        test_in.push_back(i);
    }

    std::vector<int> correct_out;
    correct_out.resize(test_in.size());

    MPI_Alltoall(test_in.data(), alltoall_size * sizeof(int), MPI_BYTE,
                 correct_out.data(), alltoall_size * sizeof(int), MPI_BYTE,
                 MPI_COMM_WORLD);

    std::vector<int> test_out;
    test_out.resize(test_in.size());

    map.exec(test_in.data(), test_out.data());

    int test_rank = 1;

    if (comm_rank == test_rank) {
        for (auto& i : test_out) {
            printf("%d ", i);
        }
        printf("\n");
    }

    if (comm_rank == test_rank) {
        for (auto& i : correct_out) {
            printf("%d ", i);
        }
        printf("\n");
    }

    MPI_Finalize();
    return 0;
}