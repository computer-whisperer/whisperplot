//
// Created by christian on 5/31/21.
//

#ifndef WHISPERPLOT_EXPLICIT_TEMPLATES_HPP
#define WHISPERPLOT_EXPLICIT_TEMPLATES_HPP

constexpr uint32_t global_interlace_factor = 1ULL << 6;
constexpr uint32_t global_num_rows = global_interlace_factor*(1ULL << 8);


template class Plotter<PlotConf{18, global_num_rows, global_interlace_factor}>;
template class Plotter<PlotConf{22, global_num_rows, global_interlace_factor}>;
template class Plotter<PlotConf{26, global_num_rows, global_interlace_factor}>;
template class Plotter<PlotConf{32, global_num_rows, global_interlace_factor}>;

#endif //WHISPERPLOT_EXPLICIT_TEMPLATES_HPP
