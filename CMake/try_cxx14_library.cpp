
// Test for C++14 library support
#include <tuple>

int main(int argc, char **argv)
{
    std::tuple<int,double> t(1, 3.5);
    // Accessing a tuple member by type is a C++14 feature
    int j = std::get<int>(t);
    return 0;
}

