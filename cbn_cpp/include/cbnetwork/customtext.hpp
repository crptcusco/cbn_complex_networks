#ifndef CUSTOMTEXT_HPP
#define CUSTOMTEXT_HPP

#include <iostream>
#include <string>

namespace cbnetwork {

class CustomText {
public:
    static void print_duplex_line() {
        std::cout << "================================================================================" << std::endl;
    }

    static void print_simple_line() {
        std::cout << "--------------------------------------------------------------------------------" << std::endl;
    }

    static void make_title(const std::string& title) {
        print_duplex_line();
        std::cout << title << std::endl;
        print_duplex_line();
    }

    static void make_sub_title(const std::string& subtitle) {
        print_simple_line();
        std::cout << subtitle << std::endl;
        print_simple_line();
    }

    static void make_sub_sub_title(const std::string& subsubtitle) {
        std::cout << ">>> " << subsubtitle << " <<<" << std::endl;
    }
};

} // namespace cbnetwork

#endif // CUSTOMTEXT_HPP
