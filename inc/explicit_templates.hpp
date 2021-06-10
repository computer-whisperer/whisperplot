//
// Created by christian on 5/31/21.
//

#ifndef WHISPERPLOT_EXPLICIT_TEMPLATES_HPP
#define WHISPERPLOT_EXPLICIT_TEMPLATES_HPP

constexpr uint32_t test_num_rows = 1ULL << 8;

template class Plotter<18, test_num_rows>;
template class Plotter<22, test_num_rows>;
template class Plotter<26, test_num_rows>;
template class Plotter<32, test_num_rows>;

#endif //WHISPERPLOT_EXPLICIT_TEMPLATES_HPP
