#pragma once

#include <stdexcept>

#include "giant.h"

namespace arithmetic
{
    class ArithmeticException : public std::runtime_error
    {
    public:
        ArithmeticException() : std::runtime_error("Arithmetic exception.") { }
        ArithmeticException(const std::string& message) : std::runtime_error(message) { }
    };

    class InvalidFFTDataException : public ArithmeticException
    {
    public:
        InvalidFFTDataException() : ArithmeticException("Invalid FFT data.") { }
    };

    class NoInverseException : public ArithmeticException
    {
    public:
        NoInverseException(Giant& divisor) : ArithmeticException("The inverse does not exist."), divisor(GiantsArithmetic::default_arithmetic())
        {
            this->divisor = divisor;
            if (this->divisor < 0)
                this->divisor.arithmetic().neg(this->divisor, this->divisor);
        }
    public:
        Giant divisor;
    };
}
