#ifndef PLOTTER_HPP
#define PLOTTER_HPP

#include <iostream>
#include <vector>
#include "penguin.hpp"
#include "buffer.hpp"
#include "park.hpp"
#include "plotter.hpp"
#include "packed_array.hpp"
#include "thread_mgr.hpp"
#include "penguin_entries.hpp"

//uint8_t getK(std::string plot_fname);

template <PlotConf conf>
class Plotter
{
public:
    uint32_t num_threads;

    class Context
    {
    public:
        Plotter<conf>* plotter;
        std::vector<uint32_t> cpu_ids;
        uint32_t numa_node;
        std::vector<void*> forward_pass_gid_penguins;
        Penguin<FwdYCEntry<conf, 5>, true>* forward_pass_final_yc_penguin;

        inline Context(Plotter<conf>* plotter_in, std::vector<uint32_t> cpu_ids_in, uint32_t numa_node_in):
            plotter(plotter_in),
            cpu_ids(cpu_ids_in),
            numa_node(numa_node_in)
        {
            forward_pass_gid_penguins.resize(6);
        }

        Penguin<FwdYCEntry<conf, -1>, true>* createFirstTable();
        template <uint8_t table_index>
        Penguin<FwdYCEntry<conf, table_index>, true>* createTable(
                std::vector<Penguin<FwdYCEntry<conf, table_index-1>, true>*> prev_penguins,
                std::atomic<uint64_t>& num_matches_out
                );
    };

    inline Plotter(std::vector<uint32_t> cpu_ids_in)
    {
        num_threads = cpu_ids_in.size();
        for (auto & numa_node : GetNUMANodesFromCpuIds(cpu_ids_in))
        {
            std::vector<uint32_t> cpu_ids_for_context;
            for (auto & cpu_id : cpu_ids_in)
            {
                if (numa_node_of_cpu(cpu_id) == (int)numa_node)
                {
                    cpu_ids_for_context.push_back(cpu_id);
                }
            }
            contexts.push_back(new Context(this, cpu_ids_for_context, numa_node));
        }
    }

    void create(std::array<uint8_t, 32> id_in);
    //void write(std::string filename, vector<uint8_t> memo);

    //void read(std::string plot_fname, bool decompress=true);
    int32_t find_proofs(uint128_t challenge);
    void find_many_proofs(uint32_t n);

    void checkFullPlot();

    std::vector<Context*> contexts;
    std::atomic<uint64_t> coordinator;

private:
    Buffer* output_buffer;
    std::vector<uint8_t> memo;
    std::array<uint8_t, 32> id;
    std::string filename;



    uint64_t final_table_begin_pointers[10];
    uint64_t pointer_table_offset;


    template <int8_t table_index>
    void readGraphTable(std::ifstream& file_stream, bool decompress);

    void readFinalTable(std::ifstream& file_stream);

    template<uint8_t table_index>
    uint8_t check_match_and_return_values(uint128_t gid_x, uint64_t& y, uint128_t& c);
};

#endif
