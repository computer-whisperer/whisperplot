//
// Created by christian on 6/14/21.
//

#ifndef WHISPERPLOT_STATUS_UPDATE_HPP
#define WHISPERPLOT_STATUS_UPDATE_HPP

#include <string>


class StatusUpdate {

public:
    static void Init();
    static void StartSeg(std::string seg_name);
    static void EndSeg();
};

#endif //WHISPERPLOT_STATUS_UPDATE_HPP
