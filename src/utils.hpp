#ifndef FASTX_COMMON_HPP
#define FASTX_COMMON_HPP

#include <cstdint>
#include <string>


/**
 * @brief convert string representation of integral number(with K/M/G suffix)
 * into integer
 * 
 * @param bases_str 
 * @return int64_t 
 */
int64_t KmgStrToInt(const std::string &str);


/**
 * @brief safe wrapper of strtol.
 * 
 * @param str C-string beginning with the representation of an integral number.
 * @param base Numerical base (radix) that determines the valid characters and their interpretation.
 * @return long 
 */
long SafeStrtol(const char *str, int base);


/**
 * @brief safe wrapper of strtod.
 * 
 * @param str C-string beginning with the representation of a floating-point number.
 * @return double 
 */
double SafeStrtod(const char *str);


#endif  // FASTX_COMMON_HPP
