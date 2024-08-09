#ifndef _SMART_MAP_HPP_
#define _SMART_MAP_HPP_

#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <map>
#include <mpi.h>
#include <unordered_map>
#include <vector>

struct smart_map_index {
    int rank = -1;
    int idx = -1;
    smart_map_index() {}
    smart_map_index(int rank_, int idx_) : rank(rank_), idx(idx_) {}
    void print() { printf("(%d %d)\n", rank, idx); }
};

struct Send {
    int source;
    int dest;

    Send() {}

    Send(int source_, int dest_) : source(source_), dest(dest_) {}

    void print() { printf("%d -> %d\n", source, dest); }

    bool is_adjacent_to(const Send& other) {
        return (source == other.source) || (dest == other.dest);
    }
};

/**
 * @class Node
 * @brief Represents a node in a graph of sends.
 */
class Node {
  private:
    /** @brief The send operation associated with this node. */
    Send m_send;
    /** @brief Neighboring nodes */
    std::vector<Node*> m_neighbors;
    /** @brief The 'color' assigned to this node (practically, the ordering of
     * the sends/recvs). */
    int m_color = -1;

  public:
    Node(Send send) : m_send(send) {}

    /**
     * @brief Get the send operation associated with this node.
     * @return The send operation.
     */
    Send send() { return m_send; }

    /**
     * @brief Add a neighboring node.
     * @param n Pointer to the neighboring node.
     */
    void add_neighbor(Node* n) { m_neighbors.push_back(n); }

    /**
     * @brief Get the color of this node (practically, the ordering of the
     * sends/recvs).
     * @return The color of this node.
     */
    int color() { return m_color; }

    /**
     * @brief Calculate and set the color for this node.
     *
     * This does greedy graph coloring. We set the color of this node to the
     * first color (starting at 0) that has not already been used by neighboring
     * nodes.
     *
     * @return The calculated color.
     */
    int calc_color() {
        int out_color = 0;
        while (true) {
            bool found = true;
            for (auto i : m_neighbors) {
                if (i->color() == out_color) {
                    found = false;
                    break;
                }
            }
            if (found) {
                break;
            }
            out_color++;
        }
        m_color = out_color;
        return m_color;
    }

    /**
     * @brief Print the send operation and color of this node.
     */
    void print() {
        printf("%d -> %d : %d", m_send.source, m_send.dest, m_color);
    }
};

/**
 * @struct Action
 * @brief Represents an action involving send and receive operations.
 */
struct Action {
    /** @brief The rank from which to receive. `-1` means that there is no recv
     * operation in this action. */
    int recv_from = -1;
    /** @brief The rank to send to. `-1` means there is no send operation in
     * this action. */
    int send_to = -1;
    /** @brief The color associated with this action. Practically, the ordering
     * of this operation. */
    int color = -1;

    /**
     * @brief Print the details of the action.
     */
    void print() {
        std::string out = "";
        if (send_to != -1) {
            out += " send to " + std::to_string(send_to);
        }
        if (recv_from != -1) {
            out += " recv from " + std::to_string(recv_from);
        }
        printf("    %d. %s\n", color + 1, out.c_str());
    }

    /**
     * @brief Check if there is a receive operation.
     * @return True if there is a receive operation, false otherwise.
     */
    bool has_recv() { return recv_from != -1; }

    /**
     * @brief Check if there is a send operation.
     * @return True if there is a send operation, false otherwise.
     */
    bool has_send() { return send_to != -1; }

    /**
     * @brief Get the rank from which to receive.
     * @return The receive rank or `-1` if there is no recv in this action.
     */
    int get_recv() {
        assert(has_recv());
        return recv_from;
    }

    /**
     * @brief Get the rank to which to send.
     * @return The send rank or `-1` if there is no send in this action.
     */
    int get_send() {
        assert(has_send());
        return send_to;
    }
};

/**
 * @class Decider
 * @brief Decides the ordering of send/recv operations to ensure each rank only
 * sends to/recvs from one other rank at a time.
 */
class Decider {
  private:
    /** @brief List of nodes. */
    std::vector<Node> m_nodes;

  public:
    /**
     * @brief Constructor to initialize the Decider with a list of sends.
     *
     * We take in a list of sends, and then create from that a list of nodes.
     * We then construct adjacency lists for the nodes.
     * Finally, we do the greedy graph coloring on the nodes.
     *
     * @param sends The list of send operations.
     */
    Decider(std::vector<Send> sends) {
        for (auto& i : sends) {
            m_nodes.push_back(Node(i));
        }
        for (size_t i = 0; i < m_nodes.size(); i++) {
            for (size_t j = 0; j < m_nodes.size(); j++) {
                Node* node_i = &m_nodes[i];
                Node* node_j = &m_nodes[j];
                if (node_i->send().is_adjacent_to(node_j->send())) {
                    node_i->add_neighbor(node_j);
                }
            }
        }
        int loops = 0;
        for (auto& i : m_nodes) {
            int tmp = i.calc_color();
            if (tmp > loops) {
                loops = tmp;
            }
        }
    }

    /**
     * @brief Print the details of all nodes.
     */
    void print() {
        for (auto& i : m_nodes) {
            i.print();
        }
    }

    /**
     * @brief Get the actions for a specific rank.
     * @param rank The rank for which to get actions.
     * @return A vector of actions.
     */
    std::vector<Action> get_actions(int rank) {
        // ordered map of color->actions (we want it to be ordered so we can
        // convert it to a vector afterwards, sorted by color/order)
        std::map<int, Action> actions_by_color;
        for (auto& i : m_nodes) {
            if (i.send().source == rank) {
                actions_by_color[i.color()].send_to = i.send().dest;
                actions_by_color[i.color()].color = i.color();
            }
            if (i.send().dest == rank) {
                actions_by_color[i.color()].recv_from = i.send().source;
                actions_by_color[i.color()].color = i.color();
            }
        }
        // convert the map to vector for simplicity later.
        std::vector<Action> out;
        for (auto& i : actions_by_color) {
            out.push_back(i.second);
        }
        return out;
    }
};

class SmartMap {
  private:
    int m_buffsz;
    MPI_Comm m_comm;
    int m_comm_size;
    int m_comm_rank;
    std::vector<int> m_send_counts;
    std::vector<int> m_recv_counts;
    std::vector<Action> m_actions;
    std::vector<int> m_send_reorder_map;
    std::vector<int> m_recv_reorder_map;
    std::vector<int> m_send_offsets;
    std::vector<int> m_recv_offsets;

    void calc_actions_and_counts() {
        m_send_counts.resize(m_comm_size);
        m_recv_counts.resize(m_comm_size);

        std::unordered_map<int, std::unordered_map<int, int>> send_map;
        std::vector<Send> sends;

        for (int rank = 0; rank < m_comm_size; rank++) {
            for (int idx = 0; idx < m_buffsz; idx++) {
                int out_rank = this->map(smart_map_index(rank, idx)).rank;
                send_map[rank][out_rank]++;
                if (rank == m_comm_rank) {
                    m_send_counts[out_rank]++;
                }
                if (out_rank == m_comm_rank) {
                    m_recv_counts[rank]++;
                }
            }
        }
        for (int source = 0; source < m_comm_size; source++) {
            for (int dest = 0; dest < m_comm_size; dest++) {
                if (send_map[source][dest]) {
                    sends.push_back(Send(source, dest));
                }
            }
        }

        m_actions = Decider(sends).get_actions(m_comm_rank);
    }

    void calc_offsets() {
        m_send_offsets.resize(m_comm_size);
        m_recv_offsets.resize(m_comm_size);

        int tot_send = 0;
        int tot_recv = 0;
        for (int i = 0; i < m_comm_size; i++) {
            m_send_offsets[i] = tot_send;
            tot_send += m_send_counts[i];

            m_recv_offsets[i] = tot_recv;
            tot_recv += m_recv_counts[i];
        }

        assert(tot_send == m_buffsz);
        assert(tot_recv == m_buffsz);
    }

    void calc_send_reorder_map() {
        std::vector<int> offsets = m_send_offsets;

        m_send_reorder_map.resize(m_buffsz);

        for (int i = 0; i < m_buffsz; i++) {
            auto out_loc = this->map(smart_map_index(m_comm_rank, i));
            m_send_reorder_map[i] = offsets[out_loc.rank];
            offsets[out_loc.rank]++;
        }
    }

    void calc_recv_reorder_map() {
        std::vector<int> offsets = m_recv_offsets;

        m_recv_reorder_map.resize(m_buffsz);

        for (int rank = 0; rank < m_comm_size; rank++) {
            if (m_recv_counts[rank]) {
                for (int i = 0; i < m_buffsz; i++) {
                    auto recv_loc = this->map(smart_map_index(rank, i));
                    if (recv_loc.rank == m_comm_rank) {
                        m_recv_reorder_map[offsets[rank]] = recv_loc.idx;
                        offsets[rank]++;
                    }
                }
            }
        }
    }

    int send_idx(int rank) { return m_send_offsets[rank]; }

    int send_count(int rank) { return m_send_counts[rank]; }

    int recv_idx(int rank) { return m_recv_offsets[rank]; }

    int recv_count(int rank) { return m_recv_counts[rank]; }

  public:
    void compile(int buffsz, MPI_Comm comm) {
        m_buffsz = buffsz;
        m_comm = comm;

        MPI_Comm_size(m_comm, &m_comm_size);
        MPI_Comm_rank(m_comm, &m_comm_rank);

        calc_actions_and_counts();
        calc_offsets();
        calc_send_reorder_map();
        calc_recv_reorder_map();
    }

    SmartMap() {}

    virtual smart_map_index map(smart_map_index idx) = 0;

    template <class T> void exec(T* in, T* out) {
        for (int i = 0; i < m_buffsz; i++) {
            out[m_send_reorder_map[i]] = in[i];
        }
        for (auto& action : m_actions) {
            if (action.has_send() && action.has_recv()) {
                MPI_Sendrecv(
                    &out[send_idx(action.get_send())],
                    send_count(action.get_send()) * sizeof(T), MPI_BYTE,
                    action.get_send(), 0, &in[recv_idx(action.get_recv())],
                    recv_count(action.get_recv()) * sizeof(T), MPI_BYTE,
                    action.get_recv(), 0, m_comm, MPI_STATUS_IGNORE);
            } else if (action.has_send()) {
                MPI_Send(&out[send_idx(action.get_send())],
                         send_count(action.get_send()) * sizeof(T), MPI_BYTE,
                         action.get_send(), 0, m_comm);
            } else if (action.has_recv()) {
                MPI_Recv(&in[recv_idx(action.get_recv())],
                         recv_count(action.get_recv()) * sizeof(T), MPI_BYTE,
                         action.get_recv(), 0, m_comm, MPI_STATUS_IGNORE);
            }
        }
        for (int i = 0; i < m_buffsz; i++) {
            out[m_recv_reorder_map[i]] = in[i];
        }
    }
};

#endif // _SMART_MAP_HPP_