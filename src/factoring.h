#pragma once

#include <memory>
#include <vector>
#include "arithmetic.h"
#include "edwards.h"
#include "inputnum.h"
#include "logging.h"
#include "file.h"

#define FACTORING_APPID 3

int factoring_main(int argc, char *argv[]);

class Factoring
{
public:
    Factoring(InputNum& input, arithmetic::GWState& gwstate, Logging& logging) : _input(input), _gwstate(gwstate), _logging(logging), _gw(gwstate), _ed(_gw)
    {
    }

    void write_points(File& file);
    bool read_points(File& file);
    bool read_state(File& file, uint64_t B1);

    std::string generate(int seed, int count);
    std::string verify(bool verify_curve);
    void modulus(int curve, File& file_result);
    bool split(int offset, int count, Factoring& result);
    bool merge(Factoring& other);
    void stage1(uint64_t B1, File& file_state, File& file_result);

    std::vector<std::unique_ptr<arithmetic::EdPoint>>& points() { return _points; }
    uint64_t B0() { return _B0; }

private:
    void write_file(File& file, char type, uint64_t B1, std::vector<std::unique_ptr<arithmetic::EdPoint>>& points);
    bool read_file(File& file, char type, int& seed, uint64_t& B0, std::vector<std::unique_ptr<arithmetic::EdPoint>>& points);

private:
    InputNum& _input;
    arithmetic::GWState& _gwstate;
    Logging& _logging;
    arithmetic::GWArithmetic _gw;
    arithmetic::EdwardsArithmetic _ed;

    int _seed;
    uint64_t _B0;
    std::vector<std::unique_ptr<arithmetic::EdPoint>> _points;
    std::vector<std::unique_ptr<arithmetic::EdPoint>> _state;
};