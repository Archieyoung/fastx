#include <climits>
#include <cstdlib>
#include <iostream>
#include <cmath>
#include "utils.hpp"


int64_t KmgStrToInt(const std::string &str) {
    if (str.empty())
    {
        std::cerr << "Error! Input bases string is empty!" << std::endl;
        std::exit(1);
    }

    char *tmp;
    int64_t bases = 0;
    switch (str.back()) {
        case 'G':
        case 'g':
            bases = static_cast<int64_t>(std::round(1000000000 *
                strtod(
                str.substr(0, str.size() - 1).c_str(), &tmp)));
            break;
        case 'M':
        case 'm':
            bases = static_cast<int64_t>(std::round(1000000 *
                strtod(
                str.substr(0, str.size() - 1).c_str(), &tmp)));
            break;
        case 'K':
        case 'k':
            bases = static_cast<int64_t>(std::round(1000 *
                strtod(
                str.substr(0, str.size() - 1).c_str(), &tmp)));
            break;
        default:
            bases = strtol(
                str.substr(0, str.size()).c_str(), &tmp, 10);
            break;
    }

    return bases;
}



long SafeStrtol(const char *str, int base)
{
    char *tmp;
    long res = strtol(str, &tmp, base);
    if (res == 0 && tmp == str) {
        std::cerr << "[SafeStrtol] Error! Can not convert " << str
            << " into long with base " << base << std::endl;
        std::exit(1);
    }

    if (res == LONG_MAX || res == LONG_MIN) {
        std::cerr << "[SafeStrtol] Error! Input value out of range! str: "
            << str << " base: " << base << ". ";
        perror("");
        std::exit(1);
    }

    return res;
}


double SafeStrtod(const char *str)
{
    char *tmp;
    double res = strtod(str, &tmp);
    if (res == 0.0 && tmp == str) {
        std::cerr << "[SafeStrtod] Error! Can not convert " << str
            << " into double" << std::endl;
        std::exit(1);
    }

    if (res == HUGE_VAL || res == -HUGE_VAL) {
        std::cerr << "[SafeStrtol] Error! Input value out of range! str: "
            << str << ". ";
        perror("");
        std::exit(1);
    }

    return res;
}
