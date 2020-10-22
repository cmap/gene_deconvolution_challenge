//
// Created by Wojciech 'Ardavel' Szałapski
//

#ifndef DPEAK_FAST_INPUT_H
#define DPEAK_FAST_INPUT_H

int readUnsignedInt(
    const char*& ptr,
    const char* const endPtr)
{
    unsigned int result = 0;

    int c;

    do
    {
        c = *ptr;
        ++ptr;

        if (ptr >= endPtr)
        {
            return -1;
        }
    } while (c < '0' || c > '9');

    do
    {
        result *= 10;
        result += c - '0';
        c = *ptr;
        ++ptr;
    } while (c >= '0' && c <= '9');

    return static_cast<int>(result);
}

#endif //DPEAK_FAST_INPUT_H
